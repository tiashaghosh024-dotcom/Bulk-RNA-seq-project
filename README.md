# Bulk-RNA-seq Analysis 
Overview: 

This project aims to show a complete Bulk RNA-Seq workflow- from raw sequencing data up to differential expression analysis using DESeq2. We tried to analyse the effect of hypoxia on gene expression in two cancer cell lines LNCAP and PC3.

For each cell lines, Hypoxia vs Normoxia was compared to identify genes that responded to low O2 levels. 

# I. Data Sourcing : 
**1. Download raw sequencing data (.sra files) using SRA toolkit.**
   Raw RNA-seq datasets submitted to the NCBI Sequence Read Archive are stored in a compressed binary format called .sra.
   To begin any bulk RNA-seq workflow, the first step is obtaining these files from SRA using tools such as ```prefetch```.
   
**2. Convert .sra → .fastq.gz using fastq-dump (or fasterq-dump)**
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

 hisat2 -q -x grch38/genome -U fastq/LNCAP_Hypoxia_S1.fastq.gz | \
 samtools sort -o alignedreads/LNCAP_Hypoxia_S1.bam

samtools index alignedreads/LNCAP_Hypoxia_S1.bam
```

For all the alignment process, create a new script file in the working directory : nano run_alignment.sh.
```bash
#!/bin/bash

# Directory containing FASTQ files
FASTQ_DIR="fastq"
GENOME_INDEX="genome_index/grch38/genome"
LOGFILE="alignment_log.txt"

# Clear or create logfile
> $LOGFILE

# List of FASTQ files
FILES=(
    "LNCAP_Hypoxia_S1.fastq.gz"
    "LNCAP_Hypoxia_S2.fastq.gz"
    "LNCAP_Normoxia_S1.fastq.gz"
    "LNCAP_Normoxia_S2.fastq.gz"
    "PC3_Hypoxia_S1.fastq.gz"
    "PC3_Hypoxia_S2.fastq.gz"
    "PC3_Normoxia_S1.fastq.gz"
    "PC3_Normoxia_S2.fastq.gz"
)

# Loop through each file
for f in "${FILES[@]}"; do
    SAMPLE_NAME=$(basename "$f" .fastq.gz)
    echo "Processing $SAMPLE_NAME at $(date)" | tee -a $LOGFILE
    START_TIME=$(date +%s)

    # Run HISAT2 alignment and Samtools sorting/indexing
    hisat2 -q -x $GENOME_INDEX -U $FASTQ_DIR/$f | \
    samtools sort -o ${SAMPLE_NAME}.bam
    samtools index ${SAMPLE_NAME}.bam

    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    echo "Finished $SAMPLE_NAME in $ELAPSED seconds at $(date)" | tee -a $LOGFILE
    echo "--------------------------------------" | tee -a $LOGFILE

done

echo "All files processed successfully at $(date)" | tee -a $LOGFILE
```
Then save and exit nano, make the script executable and run he script.


_**Output: Aligned, sorted, and indexed bam and bai files**_
<img width="1912" height="489" alt="image" src="https://github.com/user-attachments/assets/6905edda-2d51-43cc-99b5-2a7f4c9e7044" />


# VIII. Quantification:
This is for creating the gene expression count matrix.
```bash
featureCounts -S 2 \
  -a Homo_sapiens.GRCh38.113.gtf \
  -o quants/LNCAP_Normoxia_S1_featurecounts.txt \
  LNCAP_Normoxia_S1.bam
  ## to get counts for from .bam file
```
This will generete: quants/LNCAP_Normoxia_S1_featurecounts.txt
<img width="1508" height="755" alt="image" src="https://github.com/user-attachments/assets/4d1d6711-c2de-406f-8ae4-b8f2bd84f2d7" />

For all samples to be quantified:
```bash
#!/bin/bash

cd /home/tiasha/bulk_rnaseq_work

# Loop over all BAM files except temporary ones
for bam in *.bam; do

    if [[ "$bam" == *.tmp.*.bam ]]; then
        continue
    fi

    start=$(date +%s)

    echo "Processing $bam ..."

    featureCounts -s 2 \
        -a /home/tiasha/bulk_rnaseq_work/Homo_sapiens.GRCh38.113.gtf \
        -o /home/tiasha/bulk_rnaseq_work/quants/${bam%.bam}_featurecounts.txt \
        "$bam"

    end=$(date +%s)
    runtime=$(( (end - start) / 60 ))

    echo "Completed $bam in $runtime minutes."
    echo "------------------------------------"

done
```
Saved this in nano count_all.sh, made it executable and ran the script.

_**Output: Each ~50 MB, and a small .summary file for each.**_
<img width="1916" height="493" alt="image" src="https://github.com/user-attachments/assets/159a3c5e-abcd-4e54-93eb-4058fe0e1f3e" />

# IX. Creating count matrix:
It merges all 8 individual featureCounts result files into ONE big count matrix (genes × samples).
This final matrix will be used for:
1) Differential expression analysis (DESeq2)
2) Heatmaps
3) PCA plots
4) Volcano plots

This will merge files like the ones in the output of the picture above, into one .csv file and this matrix is what DESeq2 will see.
```bash
#!/usr/bin/env python
# coding: utf-8

import os
import glob
import pandas as pd
import time

path = "/home/tiasha/bulk_rnaseq_work/quants"

files = glob.glob(os.path.join(path, "*.txt"))

print("Files found:", files)

all_counts = []

for file in files:
    start_time = time.time()
    df = pd.read_csv(file, sep="\t", comment="#")
    
    sample_name = os.path.basename(file).replace("_featurecounts.txt", "")
    df = df[["Geneid", df.columns[-1]]]
    df.rename(columns={df.columns[-1]: sample_name}, inplace=True)
    
    all_counts.append(df)
    
    elapsed = (time.time() - start_time) / 60  # minutes
    print(f"Completed {sample_name} | Rows: {df.shape[0]} | Time: {elapsed:.2f} min")

counts_matrix = all_counts[0]
for df in all_counts[1:]:
    counts_matrix = counts_matrix.merge(df, on="Geneid", how="outer")

output_file = os.path.join(path, "GSE106305_counts_matrix_3011.csv")
counts_matrix.to_csv(output_file, index=False)

print("\nAll files processed!")
print("Merged matrix shape:", counts_matrix.shape)
print("Saved to:", output_file)
```

# X. Bulk RNA-seq Differential Expression Analysis- DESeq2:
**1. Installation and loading data and packages:**
   Created gse106305_project in working directory and the moved the count matrix in this folder.
   Installation:

 ```r
install.packages(c(
  "data.table",
  "dplyr",
  "tibble",
  "ggplot2",
  "forcats",
  "RColorBrewer",
  "pheatmap",
  "ggrepel",
  "stringr"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "DESeq2",
  "biomaRt",
  "msigdbr",
  "clusterProfiler",
  "ReactomePA",
  "org.Hs.eg.db",
  "fgsea",
  "AnnotationDbi",
  "DOSE",
  "enrichplot"
))
```
Set working directory inside R.
Load the data and packages:

```r
setwd("~/projects/GSE106305")   #<- changed to my path

library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(clusterProfiler)
library(fgsea)
library(DOSE)
library(enrichplot)
library(org.Hs.eg.db)
```

**2. Verification and reading counts file:**
```r
# Read counts
raw_counts <- read.csv("GSE106305_counts_matrix_3011.csv",
                       header = TRUE,
                       stringsAsFactors = FALSE)

# Check first columns & rownames
head(raw_counts)
colnames(raw_counts)[1:20]

# If Gene IDs are in first column (not row names), ensure correct shape:
# If the file has a "Geneid" column and then sample columns, set it as rownames:
if ("Geneid" %in% colnames(raw_counts)) {
  rownames(raw_counts) <- raw_counts$Geneid
  raw_counts$Geneid <- NULL
}

# Convert to matrix of integers if needed:
count_matrix <- as.matrix(raw_counts)
storage.mode(count_matrix) <- "integer"
dim(count_matrix)
```
```r
raw <- read.csv("GSE106305_counts_matrix_3011.csv")
colnames(raw)
```
_**Output:**_
<img width="1918" height="163" alt="image" src="https://github.com/user-attachments/assets/6817a709-fbd4-42c0-b30a-8a131b5565c9" />

**3. reloading CSV with proper row names and sorting columnns:**
```r
raw_counts <- read.csv("GSE106305_counts_matrix_3011.csv",
                       header = TRUE,
                       row.names = "Geneid",
                       stringsAsFactors = FALSE)

head(raw_counts)
dim(raw_counts)
```
<img width="1912" height="1031" alt="image" src="https://github.com/user-attachments/assets/0cc5c157-4e48-4f8c-b437-fe01cff9a139" />

For sorting :
```r
raw_counts <- raw_counts[, sort(colnames(raw_counts))]
```
Thw sorted order will be:
LNCAP_Hypoxia_S1
LNCAP_Hypoxia_S2
LNCAP_Normoxia_S1
LNCAP_Normoxia_S2
PC3_Hypoxia_S1
PC3_Hypoxia_S2
PC3_Normoxia_S1
PC3_Normoxia_S2

__**4. Check column sums:**__
   This checks that counts look reasonable.
   ```colSums(raw_counts)```
   <img width="1919" height="187" alt="image" src="https://github.com/user-attachments/assets/4d411ab8-90c0-4e40-a98e-8cd0d202f6c7" />

__**5. Creating metadata (colData)**__
```r
condition <- c(
  rep("LNCAP_Hypoxia", 2),
  rep("LNCAP_Normoxia", 2),
  rep("PC3_Hypoxia", 2),
  rep("PC3_Normoxia", 2)
)

my_colData <- as.data.frame(condition)
rownames(my_colData) <- colnames(raw_counts)

head(my_colData)
```
<img width="1582" height="548" alt="image" src="https://github.com/user-attachments/assets/849f43d7-4f3b-4656-a2a5-c6565c07b286" />

**6. Creating the DESeq2 dataset:**
   ```r
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData   = my_colData,
  design    = ~ condition
)

dds
```

Ran differetial expression analysis:
```dds <- DESeq(dds)```

_**Output:**_
<img width="1919" height="719" alt="image" src="https://github.com/user-attachments/assets/c0fc9bf5-d5b8-4380-b719-9524e77e022a" />
These lines mean that mean the normalization, dispersion estimation, and statistical testing are all completed.

Extracting the first differential expression results:
```r
dds$condition <- relevel(dds$condition, ref = "LNCAP_Normoxia")
```
_i.Extracting results for LNCAP Normoxia vs LNCAP Hypoxia:_

```r
res_LNCAP <- results(
  dds,
  contrast = c("condition", "LNCAP_Normoxia", "LNCAP_Hypoxia")
)

head(res_LNCAP)
summary(res_LNCAP)
```
This will give: log2FoldChange, p-values, adjusted p-values, number of up/downregulated genes.
A positive log2FoldChange → higher in Normoxia
A negative log2FoldChange → higher in Hypoxia

_ii. Extracting results for PC3 Hypoxia vs LNCAP Hypoxia:_

```r
res_PC3_Hypoxia <- results(dds,
                           name = "condition_PC3_Hypoxia_vs_LNCAP_Hypoxia")

summary(res_PC3_Hypoxia)
```

_iii. Extracting results for PC3 Normoxia vs LNCAP Hypoxia:_

```r
res_PC3_Normoxia <- results(dds, name = "condition_PC3_Normoxia_vs_LNCAP_Hypoxia")
summary(res_PC3_Normoxia)
```

_**Output:**_
```out of 44394 with nonzero total read count adjusted p-value < 0.1 LFC > 0 (up) : 2825, 6.4% LFC < 0 (down) : 3502, 7.9% outliers [1] : 0, 0% low counts [2] : 27240, 61% (mean count < 16) [1] see 'cooksCutoff' argument of ?results [2] see 'independentFiltering' argument of ?results```

**7. Ordering results and saving DE results to a CSV file:**
```r
res_LNCAP_ordered <- res_LNCAP[order(res_LNCAP$padj), ]
head(res_LNCAP_ordered)
```
<img width="1919" height="530" alt="image" src="https://github.com/user-attachments/assets/d4bb4e18-4974-4f3d-ace3-c6761c922456" />
This shows the top most significant genes.

Saving DE results to a csv file:
```r
write.csv(
  as.data.frame(res_LNCAP_ordered),
  file = "LNCAP_Normoxia_vs_Hypoxia_DEGs.csv"
)
```

**8. Volcano plot:**
   Used to quickly identify genes that show both large expression changes and strong statistical significance in differential expression analysis.

  ```r
res_df <- as.data.frame(res_LNCAP_ordered)
res_df <- na.omit(res_df)
res_df$gene <- rownames(res_df)

res_df$regulation <- "Not Significant"
res_df$regulation[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Upregulated"
res_df$regulation[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Downregulated"

qp <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = regulation)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c(
    "Upregulated" = "#FEA405",
    "Downregulated" = "purple",
    "Not Significant" = "gray"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value")

ggsave("vp_lncap.png", plot = qp, width = 8, height = 6, dpi = 300)
```

**9. PCA Analysis:**
   Used in bulk RNA-seq to quickly visualize sample-to-sample similarity and identify major sources of variation such as batch effects or clear group separation.

  ```r
vsd <- vst(dds, blind = TRUE)
```
   This created the vst object.
  ```r
plot_PCA = function(vsd.obj) {
  pcaData <- plotPCA(vsd.obj, intgroup = c("condition"), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  ggplot(pcaData, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    labs(
      x = paste0("PC1: ", percentVar[1], "% variance"),
      y = paste0("PC2: ", percentVar[2], "% variance"),
      title = "PCA Plot colored by condition"
    ) +
    ggrepel::geom_text_repel(aes(label = name), color = "black")
}

png(filename = "pcab.png", width = 2000, height = 2000, res = 300)
plot_PCA(vsd)
dev.off()
```

**10. Sample to sample distance heatmap:**
    Shows how similar or different the samples are from one another, helping detect outliers, batch effects and whether biological replicates cluster together as expected.
  ```r
plotDists = function (vsd.obj) {
  sampleDists <- dist(t(assay(vsd.obj)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd.obj$condition)
  
  colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(55)
  
  pheatmap::pheatmap(
    sampleDistMatrix,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists,
    col = colors,
    fontsize_row = 4,
    fontsize_col = 4,
    fontsize_legend = 4,
    fontsize = 4
  )
}

png(filename = "sampleheatmap1.png", width = 1000, height = 900, res = 300)
plotDists(vsd)
dev.off()
```

11. **Gene Set Enrichment Analysis (GSEA) of LNCaP Hypoxia RNA-seq:**
    i. Prepare pathways for fgsea (installed msigdbr)
  ```r
# Use new msigdbr syntax
m_df <- msigdbr(
  species    = "Homo sapiens",
  collection = "H"
)

# Convert to fgsea list format
pathways_hallmark <- m_df %>%
  split(x = .$ensembl_gene, f = .$gs_name)

length(pathways_hallmark)
head(names(pathways_hallmark))
```
 <img width="1919" height="784" alt="image" src="https://github.com/user-attachments/assets/a978f390-a4c5-48c5-909e-21a6cb7fc649" />

   ii. Run fgsea
   ```r
    fgsea_res <- fgsea(
    pathways = pathways_hallmark,
     stats = gene_ranks,
     minSize = 15,
     maxSize = 500,
      nperm = 1000
       )

    # Sort results
    fgsea_res <- fgsea_res[order(fgsea_res$pval), ]
    head(fgsea_res)
```
iii. Save GSEA table
   ```r
    write.csv(fgsea_res, "fgsea_hallmark_results.csv", row.names = FALSE)
```

iv. Plot NES-ranked pathways (barplot)
```r
library(ggplot2)

p_gsea <- ggplot(fgsea_res[1:20, ], 
                 aes(x = reorder(pathway, NES), 
                     y = NES, 
                     fill = padj)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_continuous(type = "viridis") +
  theme_bw() +
  labs(
    title = "Top 20 Hallmark Pathways (fgsea)"
  )

ggsave("gsea_barplot.png", 
       p_gsea, 
       width = 8, 
       height = 6, 
       dpi = 300)
```
v. Enrichment plot for the top pathway
```r
top_pathway <- fgsea_res$pathway[1]

png("fgsea_top_pathway.png", width = 2000, height = 1200, res = 300)
plotEnrichment(pathways_hallmark[[top_pathway]], gene_ranks) +
  ggtitle(paste("Top Enriched Pathway:", top_pathway))
dev.off()
```
vi. Barplot of top enriched pathways
```r
library(ggplot2)

top20 <- fgsea_res[1:20, ]

png("fgsea_barplot.png", width = 2000, height = 1400, res = 300)
ggplot(top20, aes(x = reorder(pathway, NES), y = NES, fill = padj)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(
    title = "Top 20 Enriched Hallmark Pathways",
    x = "Pathway",
    y = "Normalized Enrichment Score (NES)",
    fill = "Adj P-value"
  ) +
  theme_bw(base_size = 14)
dev.off()
```
vii. Leading edge enrichment plot
 ```r
top_pathway <- fgsea_res$pathway[1]
top_pathway
```
```r
png("fgsea_leading_edge.png", width = 2000, height = 1400, res = 300)
plotEnrichment(
  pathways_hallmark[[top_pathway]],
  gene_ranks
) + 
  ggtitle(paste("Leading Edge Plot for:", top_pathway))
dev.off()
```

**12. MA Plot for LNCAP contrast:**
```r
png("MAplot_LNCAP.png", width = 2000, height = 1600, res = 300)
plotMA(res_LNCAP, ylim = c(-5, 5), cex = 0.6)
dev.off()
```
**13. Creating boxplots:**
```r
# extract normalized counts
df_counts <- plotCounts(dds, gene = gene_id, returnData = TRUE)

library(ggplot2)

png("boxplot_IGFBP1_ggplot.png", width = 1600, height = 1200, res = 300)

ggplot(df_counts, aes(x = condition, y = count, fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.7) +
  theme_bw(base_size = 18) +
  ggtitle(paste("Expression of", "IGFBP1")) +
  ylab("Normalized counts")

dev.off()
```
**14. Heatmap — Top 30 Most Variable Genes**
    -used on mostly all RNA-seq QC pipelines
```r
# variance for each gene
gene_variances <- apply(assay(vsd), 1, var)

# Select top 30 most variable genes
top30 <- names(sort(gene_variances, decreasing = TRUE))[1:30]

mat_topvar <- assay(vsd)[top30, ]
mat_topvar <- mat_topvar - rowMeans(mat_topvar)  # center rows

annotation_col <- data.frame(Condition = vsd$condition)
rownames(annotation_col) <- colnames(mat_topvar)

png("top30_variable_genes_heatmap.png", width = 2000, height = 1600, res = 300)
pheatmap(
  mat_topvar,
  annotation_col = annotation_col,
  scale = "row",
  show_rownames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  main = "Top 30 Most Variable Genes"
)
dev.off()
```
Heatmap for Top 30 DE Genes (LFC-based):
```r
# Clean DE results
res_clean <- res_lncap[!is.na(res_lncap$log2FoldChange), ]
res_clean <- res_clean[order(abs(res_clean$log2FoldChange), decreasing = TRUE), ]

# Top 30 DE genes
top_genes_LFC <- rownames(res_clean)[1:30]

# Extract only genes present in VSD matrix
top_genes_LFC <- top_genes_LFC[top_genes_LFC %in% rownames(vsd)]

mat_lfc <- assay(vsd)[top_genes_LFC, ]
mat_lfc <- mat_lfc - rowMeans(mat_lfc)

annotation_col <- data.frame(Condition = vsd$condition)
rownames(annotation_col) <- colnames(mat_lfc)

png("top30_DE_genes_heatmap.png", width = 2000, height = 1600, res = 300)
pheatmap(
  mat_lfc,
  annotation_col = annotation_col,
  scale = "row",
  fontsize_row = 8,
  fontsize_col = 10,
  main = "Top 30 DE Genes by |LFC|"
)
dev.off()
```
**15. Zero count distribution table:**
Counts per gene = how many samples have zero counts. Saves a small summary table.
```r
# Zero-count distribution
counts_mat <- as.matrix(assay(dds))   # or raw count matrix you used
zero_counts <- rowSums(counts_mat == 0)

zc_table <- as.data.frame(table(zero_counts))
colnames(zc_table) <- c("zeros_in_samples", "gene_count")

write.csv(
  zc_table,
  "zero_count_distribution.csv",
  row.names = FALSE
)

summary(zero_counts)
```
**16. QC figure: Density plots raw vst**
```r
# Raw vs VST Density Plots (Grid View) 

library(ggplot2)
library(gridExtra)

# Convert raw counts to data frame for density plotting
raw_df <- as.data.frame(assay(dds)[, ])
vst_df <- as.data.frame(assay(vsd)[, ])

# Function to create density plots for raw + VST
make_density_plots <- function(sample_name) {
  df_raw <- data.frame(value = raw_df[[sample_name]])
  df_vst <- data.frame(value = vst_df[[sample_name]])
  
  p_raw <- ggplot(df_raw, aes(x = value)) +
    geom_density(color = "red") +
    ggtitle(paste("Raw - Sample", sample_name)) +
    theme_bw()
  
  p_vst <- ggplot(df_vst, aes(x = value)) +
    geom_density(color = "blue") +
    ggtitle(paste("VST - Sample", sample_name)) +
    theme_bw()
  
  list(p_raw, p_vst)
}

plots <- unlist(lapply(colnames(raw_df), make_density_plots), recursive = FALSE)

png("density_plots_raw_vst_grid.png", width = 4000, height = 4000, res = 300)
grid.arrange(grobs = plots, ncol = 4)
dev.off()
```
**17. Second PCA plot (After DESeq normalization):**
```r
vsd2 <- vst(dds, blind = TRUE)

plot_PCA2 <- function(vsd.obj) {
  pcaData <- plotPCA(vsd.obj, intgroup = c("condition"), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(aes(label = name), color = "black") +
    labs(
      x = paste0("PC1: ", percentVar[1], "% variance"),
      y = paste0("PC2: ", percentVar[2], "% variance"),
      title = "PCA Plot (After DESeq Normalization)"
    ) +
    theme_minimal()
}

png("PCA_after_DESeq.png", width = 2000, height = 2000, res = 300)
plot_PCA2(vsd2)
dev.off()
```
**18. Volcano plot on the DEGs (ordered DE results):**
```r
res_df <- as.data.frame(reslncapOrdered)
res_df <- na.omit(res_df)
res_df$gene <- rownames(res_df)

res_df$regulation <- "Not Significant"
res_df$regulation[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Upregulated"
res_df$regulation[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Downregulated"

library(ggplot2)

volcano2 <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = regulation)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c(
    "Upregulated" = "#FEA405",
    "Downregulated" = "purple",
    "Not Significant" = "gray"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Volcano Plot (DESeq2 Ordered Results)",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value"
  )

png("Volcano_LNCAP_Ordered.png", width = 2000, height = 1600, res = 300)
print(volcano2)
dev.off()
```
**19. Waterfall Plot:**
```r
library(ggplot2)
library(stringr)

# Wrap long labels to avoid clipping
top20$pathway_wrap <- str_wrap(top20$pathway_short, width = 25)

# Save high-resolution, very wide PNG
png("hallmark_waterfall_final_full.png", width = 3800, height = 2400, res = 300)

ggplot(top20, aes(x = reorder(pathway_wrap, NES), y = NES, fill = NES > 0)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#2E86AB", "FALSE" = "#E74C3C"), guide = FALSE) +
  labs(
    title = "Hallmark Pathways Altered in LNCaP Hypoxia",
    x = "Pathway",
    y = "Normalized Enrichment Score (NES)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 30, face = "bold"),
    axis.text.y = element_text(size = 14, lineheight = 0.9),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    plot.margin = margin(40, 60, 40, 60)
  )

dev.off()
```

**ALL .PNG FILES (RESULTS) OF DESEQ2 ARE UPLOADED IN THE FILES SECTION.**

# Conclusion:
This project successfully demonstrates a complete bulk RNA-seq analysis pipeline, starting from raw sequencing data in the SRA database and ending with biologically meaningful interpretations using DESeq2 and pathway enrichment tools. Beginning with FASTQ downloads and quantification steps, I generated a high-quality count matrix and performed QC to ensure sample consistency. Normalization and differential expression analysis using DESeq2 revealed key hypoxia-responsive genes in LNCaP prostate cancer cells, supported by clear visual outputs such as PCA plots, heatmaps, and volcano graphs. Integration of gene annotation and pathway enrichment methods like FGSEA Hallmark analysis—provided deeper insights into molecular pathways affected under hypoxic conditions, such as metabolic adaptation, stress signaling and cell-state transitions. Overall, the workflow illustrates how raw sequencing data can be transformed into a coherent biological story through careful computational analysis.

# Acknowledgement:
To Smriti Arora, her workshop helped me understand the basics of Bulk RNA seq and gave me the confidence to perform hands on pipeline for this project.
Reference of the tutorial, used for analysis: https://github.com/erilu/bulk-rnaseq-analysis

# Author details:
Tiasha Ghosh
tiashaghosh024@gmail.com

    









    





   
    


   





   







   











































    

 


















   
   













