# Bulk-RNA-seq-project
Overview: 

This project aims to show a complete Bulk RNA-Seq workflow- from raw sequencing data up to differential expression analysis using DESeq2. We tried to analyse the effect of hypoxia on gene expression in two cancer cell lines LNCAP and PC3.

For each cell lines, Hypoxia vs Normoxia was compared to identify genes that responded to low O2 levels. 

# I. Data Sourcing : 
# 1. Download raw sequencing data (.sra files) using SRA toolkit.
   Raw RNA-seq datasets submitted to the NCBI Sequence Read Archive are stored in a compressed binary format called .sra.
   To begin any bulk RNA-seq workflow, the first step is obtaining these files from SRA using tools such as ```prefetch```.
# 2. Convert .sra → .fastq.gz using fastq-dump (or fasterq-dump)
   The .sra format is not directly usable by downstream RNA-seq tools such as FastQC or HISAT2.These tools require reads in FASTQ format, which contains: nucleotide sequences, per-base quality scores,
   identifiers.
   Therefore, each .sra file must be converted into .**fastq** (or better, .**fastq.gz** for compression).

  ```bash
sudo apt install sra-toolkit
prefetch SRR7179504

# Converting to fastq
fastq-dump --outdir fastq --gzip --skip-technical --readids \
--read-filter pass --dumpbase --split-3 --clip SRR7179504.sra
```
**_Output: FASTQ files (raw sequencing reads_)**
```bash
tiasha@Tiasha:~/bulk_rnaseq_work$ ls sra_files

SRR7179504.sra   SRR7179505.sra   SRR7179506.sra   SRR7179507.sra
SRR7179508.sra   SRR7179509.sra   SRR7179510.sra   SRR7179511.sra
SRR7179520.sra   SRR7179521.sra   SRR7179522.sra   SRR7179523.sra
SRR7179524.sra   SRR7179525.sra   SRR7179526.sra
SRR7179527.sra   SRR7179536.sra   SRR7179537.sra
SRR7179540.sra   SRR7179541.sra
```
Also, FASTQ files are extremely large, so compressing them with --gzip reduces storage space without affecting read quality or downstream analysis.

# II. Automated SRA downloading and FASTQ generation:
Processing many RNA-seq samples manually can be slow, repetitive, and error-prone.
To make the workflow reproducible and efficient, a Python automation script was used to:

1. Download all SRA files automatically (prefetch)
2. Convert each .sra file into compressed FASTQ format (fastq-dump --gzip)
3. Record how long each step took
This automation ensures consistent processing, prevents typos in SRR IDs and saves a significant amount of manual effort when working with many samples.
```bash
import subprocess
import time

sra_numbers = [
    "SRR7179504", "SRR7179505", "SRR7179506", "SRR7179507",
    "SRR7179508", "SRR7179509", "SRR7179510", "SRR7179511",
    "SRR7179520", "SRR7179521", "SRR7179522", "SRR7179523",
    "SRR7179524", "SRR7179525", "SRR7179526", "SRR7179527",
    "SRR7179536", "SRR7179537", "SRR7179540", "SRR7179541"
]

# Download all SRA files
for sra_id in sra_numbers:
    print("\n=== Downloading:", sra_id, "===")
    prefetch_cmd = f"prefetch {sra_id}"
    print("Command:", prefetch_cmd)

    start_time = time.time()
    subprocess.call(prefetch_cmd, shell=True)
    end_time = time.time()

    elapsed_min = (end_time - start_time) / 60
    print(f"⏱ Download time for {sra_id}: {elapsed_min:.2f} minutes")

# Convert all SRA files to FASTQ
for sra_id in sra_numbers:
    sra_path = f"/{sra_id}/{sra_id}.sra"
    print("\n=== Generating FASTQ for:", sra_id, "===")
    fastq_dump_cmd = (
        f"fastq-dump --outdir fastq --gzip --skip-technical "
        f"--readids --read-filter pass --dumpbase --split-3 --clip {sra_path}"
    )
    print("Command:", fastq_dump_cmd)

    start_time = time.time()
    subprocess.call(fastq_dump_cmd, shell=True)
    end_time = time.time()

    elapsed_min = (end_time - start_time) / 60
    print(f"⏱ FASTQ generation time for {sra_id}: {elapsed_min:.2f} minutes")
```

# III. Quality Control Using FastQC

Once FASTQ files were generated, quality assessment was performed using FastQC.
Quality control is essential before alignment because it allows detection of:
1.Poor-quality reads
2.Adapter contamination
3.Over-represented sequences
4.GC content deviations
5.Per-base quality drops

Running FastQC ensures that the data is suitable for mapping and downstream differential expression analysis.
```bash
 # Run FastQC
 fastqc fastq/*.fastq.gz -o fastqc_results/ --threads 8
```
**_Output: HTML reports showing sequencing quality _**
<img width="1919" height="250" alt="Screenshot 2025-11-22 220153" src="https://github.com/user-attachments/assets/320f8bee-9141-43dd-98f7-f2814ee458fb" />

 **_Quality checking FastQC for a single sequence_**
<img width="1919" height="900" alt="image" src="https://github.com/user-attachments/assets/483e03b0-968f-4321-bc5f-8046e824cd67" />

**Running MultiQC:**
MultiQC aggregates all individual FastQC reports into a single unified report, making it easier to compare quality metrics across multiple samples in one place.
This helps quickly identify sample-specific issues, batch effects, or global sequencing problems without opening dozens of separate HTML files.
```bash
# Aggregate with MultiQC
multiqc fastqc_results/ -o multiqc_report/
```
_**Output: MultiQC report**_
<img width="1911" height="894" alt="image" src="https://github.com/user-attachments/assets/590461b1-9833-4aca-a04f-e92f46fe89c0" />

# IV. Read Trimming (Optional)
Helps in removing low-quality bases, sequencing adapters, or technical artifacts from raw reads, improving alignment accuracy and downstream gene-expression quantification. Although not always required, trimming can help when FastQC shows adapter contamination or poor quality at read ends.
```# Trimming low-quality bases using Trimmomatic (Single-End)
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar SE \
  -threads 4 \
  fastq/SRR7179504.fastq.gz \
  fastq/SRR7179504_trimmed.fastq.gz \
  TRAILING:10 \
  -phred33
```
# V. Combine Technical Replicates: concatenating FastQ files 
 Each LNCaP sample had four separate SRA runs, so running the file_name.py script produced four FASTQ files per sample. To get one final FASTQ file for each sample, we simply concatenate the four files using the cat command. 
```bash
cat SRR7179504_pass.fastq.gz SRR7179505_pass.fastq.gz SRR7179506_pass.fastq.gz SRR7179507_pass.fastq.gz  > LNCAP_Normoxia_S1.fastq.gz
cat SRR7179508_pass.fastq.gz SRR7179509_pass.fastq.gz SRR7179510_pass.fastq.gz SRR7179511_pass.fastq.gz  > LNCAP_Normoxia_S2.fastq.gz
cat SRR7179520_pass.fastq.gz SRR7179521_pass.fastq.gz SRR7179522_pass.fastq.gz SRR7179523_pass.fastq.gz  > LNCAP_Hypoxia_S1.fastq.gz
cat SRR7179524_pass.fastq.gz SRR7179525_pass.fastq.gz SRR7179526_pass.fastq.gz SRR7179527_pass.fastq.gz  > LNCAP_Hypoxia_S2.fastq.gz
```
For the PC3 samples, each one has only a single FASTQ file. So instead of merging, we rename the files from their SRR IDs to their actual sample names using the `mv` command.
```bash
mv SRR7179536_pass.fastq.gz PC3_Normoxia_S1.fastq.gz
mv SRR7179537_pass.fastq.gz PC3_Normoxia_S2.fastq.gz
mv SRR7179540_pass.fastq.gz PC3_Hypoxia_S1.fastq.gz
mv SRR7179541_pass.fastq.gz PC3_Hypoxia_S2.fastq.gz
```
The individual SRA runs are'nt required anymore, so we can remove them all using the command ```rm SRR*```, which removes all the files in the folder that begin with “SRR". Now, the folder contains a total of 8 FASTQ files: 4 for LNCaP and 4 for PC3.

_**Output:8 final FASTQ files**_
<img width="1897" height="119" alt="image" src="https://github.com/user-attachments/assets/157eec6d-a92d-4b87-9415-322ae414c3f5" />

# VI. Genome Reference and Annotation :
Retrieve the HISAT2-ready GRCh38 index:
```bash
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar -xvzf grch38_genome.tar.gz
```
Download Ensembl GTF annotation:
```bash
wget https://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.gtf.gz
gunzip Homo_sapiens.GRCh38.113.gtf.gz
```

# VII. Installation of HISAT2 and Alignment (HISAT2 → Samtools):
 Align the reads and convert to sorted & indexed BAM:
 ```bash
sudo apt install hisat2
sudo apt install samtools
```
```hisat2 -q -x grch38/genome -U fastq/LNCAP_Hypoxia_S1.fastq.gz | \
  samtools sort -o alignedreads/LNCAP_Hypoxia_S1.bam

samtools index alignedreads/LNCAP_Hypoxia_S1.bam
```
_**Output: Aligned, sorted, and indexed bam files**_
<img width="1912" height="489" alt="image" src="https://github.com/user-attachments/assets/6905edda-2d51-43cc-99b5-2a7f4c9e7044" />

 


















   
   













