#!/usr/bin/bash


# DATA PRE-PROCESSING
# download dataset
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz

# download and unzip reference genome
wget https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz
gunzip hg19.chr5_12_17.fa.gz

# quality control with fastqc and multiqc
fastqc *.fastq.gz -o fastqc_reports/
multiqc fastqc_reports -o fastqc_reports/

# trimming with trimmomatic and quality control of trimmed reads
mkdir -p trimmed_reads
for sample in `cat list.txt`
do
       trimmomatic PE -threads 8 raw_data/${sample}_r1_chr5_12_17.fastq.gz raw_data/${sample}_r2_chr5_12_17.fastq.gz \
               trimmed_reads/${sample}_r1_paired.fq.gz trimmed_reads/${sample}_r1_unpaired.fq.gz \
               trimmed_reads/${sample}_r2_paired.fq.gz trimmed_reads/${sample}_r2_unpaired.fq.gz \
               ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:keepBothReads \
               LEADING:3 TRAILING:10 MINLEN:25
       
       fastqc  trimmed_reads/${sample}_r1_paired.fq.gz  trimmed_reads/${sample}_r2_paired.fq.gz \
                 -o trimmed_reads/Fastqc_results
done 
multiqc  trimmed_reads/Fastqc_results  -o trimmed_reads/Fastqc_results


# READ MAPPING
# index reference
bwa index refgenome/hg19.chr5_12_17.fa

mkdir mapping

# perform alignment with BWA
bwa mem -R '@RG\tID:231335\tSM:Normal' refgenome/hg19.chr5_12_17.fa trimmed_reads/SLGFSK-N_231335_r1_paired.fq.gz \
      trimmed_reads/SLGFSK-N_231335_r2_paired.fq.gz > mapping/SLGFSK-N_231335.sam

bwa mem -R '@RG\tID:231336\tSM:Tumor' refgenome/hg19.chr5_12_17.fa trimmed_reads/SLGFSK-T_231336_r1_paired.fq.gz \
       trimmed_reads/SLGFSK-T_231336_r2_paired.fq.gz > mapping/SLGFSK-T_231336.sam

# convert SAM to BAM, sort and index
for sample in `cat list.txt`
do
        samtools view -@ 20 -S -b mapping/${sample}.sam | samtools sort -@ 32 > mapping/${sample}.sorted.bam
        samtools index mapping/${sample}.sorted.bam
done

# DATA POST-PROCESSING
# filter mapped reads
for sample in `cat list.txt`
do
        samtools view -q 1 -f 0x2 -F 0x8 -b mapping/${sample}.sorted.bam > mapping/${sample}.filtered.bam
done

# remove duplicates with markdup
for sample in `cat list.txt`
do
	samtools sort -n -o ${sample}.namesort.bam ${sample}.filtered.bam
        samtools fixmate -m ${sample}.namesort.bam ${sample}.fixmate.bam
        samtools sort -@ 32 -o ${sample}.positionsort.bam ${sample}.fixmate.bam
        samtools markdup -@32 -r ${sample}.positionsort.bam ${sample}.clean.bam
done

# Left-align reads around indels
for sample in `cat list.txt`
do      
        cat mapping/${sample}.clean.bam  | bamleftalign -f refgenome/hg19.chr5_12_17.fa -m 5 -c > mapping/${sample}.leftAlign.bam
done

# recalibrate read mapping qualities
for sample in `cat list.txt`
do
        samtools calmd -@ 32 -b mapping/${sample}.leftAlign.bam refgenome/hg19.chr5_12_17.fa > mapping/${sample}.recalibrate.bam
done

# refilter based on read mapping qualities
for sample in `cat list.txt`
do
        bamtools filter -in mapping/${sample}.recalibrate.bam -mapQuality "<=254" > mapping/${sample}.refilter.bam
done

# VARIANT CALLING AND CLASSIFICATION
# Convert BAM files to pileup
mkdir Variants
for sample in `cat list.txt`
do
        samtools mpileup -f refgenome/hg19.chr5_12_17.fa mapping/${sample}.refilter.bam --min-MQ 1 --min-BQ 28 \
                > Variants/${sample}.pileup
done

# call variants with VarScan Somatic
varscan somatic Variants/SLGFSK-N_231335.pileup \         Variants/SLGFSK-T_231336.pileup Variants/SLGFSK \         --normal-purity 1 --tumor-purity 0.5 --output-vcf 1

# merge vcf files with bcftools
bgzip Variants/SLGFSK.snp.vcf > Variants/SLGFSK.snp.vcf.gz
bgzip Variants/SLGFSK.indel.vcf > Variants/SLGFSK.indel.vcf.gz
tabix Variants/SLGFSK.snp.vcf.gz
tabix Variants/SLGFSK.indel.vcf.gz
bcftools merge --force-samples Variants/SLGFSK.snp.vcf.gz Variants/SLGFSK.indel.vcf.gz > Variants/SLGFSK.vcf


# VARIANT ANNOTATION
# annotate variants with snpEff database
snpEff hg19 Variants/SLGFSK.vcf > Variants/SLGFSK.ann.vcf


# The above steps were replicated and further annotation and data processing was carried out on Galaxy server
