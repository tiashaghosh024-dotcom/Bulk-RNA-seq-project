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
**Output: FASTQ files (raw sequencing reads)**
```bash
tiasha@Tiasha:~/bulk_rnaseq_work$ ls sra_files

SRR7179504.sra   SRR7179505.sra   SRR7179506.sra   SRR7179507.sra
SRR7179508.sra   SRR7179509.sra   SRR7179510.sra   SRR7179511.sra
SRR7179520.sra   SRR7179521.sra   SRR7179522.sra   SRR7179523.sra
SRR7179524.sra   SRR7179525.sra   SRR7179526.sra
SRR7179527.sra   SRR7179536.sra   SRR7179537.sra
SRR7179540.sra   SRR7179541.sra
```




   
   











