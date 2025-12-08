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
# Output: HTML reports showing sequencing quality 
<img width="1919" height="250" alt="Screenshot 2025-11-22 220153" src="https://github.com/user-attachments/assets/320f8bee-9141-43dd-98f7-f2814ee458fb" />

# Quality checking FastQC for a single sequence
<img width="1919" height="900" alt="image" src="https://github.com/user-attachments/assets/483e03b0-968f-4321-bc5f-8046e824cd67" />

# IV. 













   
   











