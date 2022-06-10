#!/usr/bin/bash
# Count the number of sequences in DNA.fa
wget https://raw.githubusercontent.com/HackBio-Internship/wale-home-tasks/main/DNA.fa
grep -v "^>" DNA.fa | wc -c

#Write a one-line command in Bash to get the total A, T, G & C counts for all the sequences in DNA.fa
grep -v "^>" DNA.fa | grep -E -o "A|T|G|C" | wc -l

#Set up a conda (anaconda, miniconda or mini forge) environment on your terminal

sudo apt update
sudo apt upgrade
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
bash Miniconda3-py39_4.12.0-Linux-x86_64.sh
conda create --name HackBio-env
conda activate HackBio-env

#Install your three softwares
conda install -c bioconda fastqc fastp multiqc

#Downloads some (>2) sample datasets from here: https://github.com/josoga2/yt-dataset/tree/main/dataset

wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/ACBarrie_R1.fastq.gz?raw=true -O ACBarrie_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/ACBarrie_R2.fastq.gz?raw=true -O ACBarrie_R2.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Baxter_R1.fastq.gz?raw=true -O Baxter_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Baxter_R2.fastq.gz?raw=true -O Baxter_R2.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Drysdale_R1.fastq.gz?raw=true -O Drysdale_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Drysdale_R2.fastq.gz?raw=true -O Drysdale_R2.fastq.gz

#Create a folder called output
mkdir output

#Implement the three software on the downloaded files (sample datasets) and send all outputs to the output folder.

#USING FASTQC FOR QUALITY CONTROL (AND SAVING REPORTS TO OUTPUT FOLDER)
fastqc ACBarrie_R1.fastq.gz -O output/
fastqc ACBarrie_R2.fastq.gz -O output/
fastqc Baxter_R1.fastq.gz -O output/
fastqc Baxter_R2.fastq.gz -O output/
fastqc Drysdale_R1.fastq.gz -O output/
fastqc Drysdale_R2.fastq.gz -O output/

#USING FASTP TO TRIM ADAPTERS
touch trim.sh
which bash
nano trim.sh

#in trim.sh
#!/usr/bin/bash
mkdir qc_reads
SAMPLES=(
  "ACBarrie"
  "Alsen"
  "Baxter"
  "Chara"
  "Drysdale"
)

for SAMPLE in "${SAMPLES[@]}"; do

  fastp \
    -i "$PWD/${SAMPLE}_R1.fastq.gz" \
    -I "$PWD/${SAMPLE}_R2.fastq.gz" \
    -o "qc_reads/${SAMPLE}_R1.fastq.gz" \
    -O "qc_reads/${SAMPLE}_R2.fastq.gz" \
    --html "qc_reads/${SAMPLE}_fastp.html" 
done

bash trim.sh
mv qc_reads trimmed_datasets/
cd trimmed_datasets/
cp ACBarrie_fastp.html Baxter_fastp.html Drysdale_fastp.html ../output/

#USING MULTIQC TO ASSEMBLE QUALITY CONTROL REPORTS
multiqc output/
cp multiqc_report.html output/
