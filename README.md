# Bulk-RNA-seq-project
Overview: 

This project aims to show a complete Bulk RNA-Seq workflow- from raw sequencing data up to differential expression analysis using DESeq2. We tried to analyse the effect of hypoxia on gene expression in two cancer cell lines LNCAP and PC3.

For each cell lines, Hypoxia vs Normoxia was compared to identify genes that responded to low O2 levels. 

# Data Sourcing : 
# 1. Download raw sequencing data (.sra files) using SRA toolkit.
   Raw RNA-seq datasets submitted to the NCBI Sequence Read Archive are stored in a compressed binary format called .sra.
   To begin any bulk RNA-seq workflow, the first step is obtaining these files from SRA using tools such as ```prefetch```.
# 2. Convert .sra → .fastq.gz using fastq-dump (or fasterq-dump)
The .sra format is not directly usable by downstream RNA-seq tools such as FastQC or HISAT2.These tools require reads in FASTQ format, which contains: nucleotide sequences, per-base quality scores,
   identifiers.
   Therefore, each .sra file must be converted into .**fastq** (or better, .**fastq.gz** for compression).
```Downloading SRA
  sudo apt install sra-toolkit
  prefetch SRR7179504
```Converting to fastq
fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip SRR7179504.sra











