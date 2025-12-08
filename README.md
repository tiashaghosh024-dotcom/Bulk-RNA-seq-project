# Bulk-RNA-seq-project

```bash
# Download SRA
prefetch SRR7179504

# Convert to FASTQ
fastq-dump --outdir fastq --gzip --skip-technical --readids \
--read-filter pass --dumpbase --split-3 --clip SRR7179504/SRR7179504.sra
```

