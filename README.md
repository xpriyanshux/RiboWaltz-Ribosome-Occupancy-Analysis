# Ribosome Profiling Analysis with riboWaltz (RiboLace Experiments)

**Purpose:** Process BAM files from RiboLace experiments to calculate **P-site offsets**, summarize **gene-specific ribosome occupancy**, generate plots, and perform basic statistics.

---

## 📌 Overview

This workflow leverages **riboWaltz** to analyze ribosome profiling data from RiboLace experiments. The pipeline includes:

- Loading and preprocessing BAM files  
- Filtering duplicate reads  
- Calculating P-site offsets  
- Annotating reads to transcripts/genes  
- Summarizing P-sites per gene  
- Plotting ribosome occupancy per sample and across groups  
- Performing basic statistical tests (t-tests)  

Key genes analyzed in this example: `LDHA`, `VEGFA`, `SLC2A1` (GLUT1), `AK4`.

---

## ⚙️ Requirements

- **R (≥ 4.0)**  
- **riboWaltz**  
- **GenomicFeatures**, **TxDb.Hsapiens.UCSC.hg38.knownGene**  
- **org.Hs.eg.db**, **AnnotationDbi**  
- **tidyverse packages**: `ggplot2`, `dplyr`, `tidyr`  
- **fs** (for file operations)  

---

---

## 🚀 Workflow

### Step 0: Load libraries
Import `riboWaltz`, `GenomicFeatures`, annotation packages, and tidyverse dependencies.

### Step 1: Define BAM files
- Specify the BAM directory and list files.  
- Assign sample names and save a `bam_info.csv` table for reference.

### Step 2: Load annotation
- Create transcript annotation using `TxDb.Hsapiens.UCSC.hg38.knownGene`.  
- Required for mapping reads to P-sites.

### Step 3–5: Process BAM files
- Convert BAM to reads list using `bamtolist()`.  
- Filter duplicate reads (protocol-dependent).  
- Calculate P-site offsets (`psite()`) and annotate reads (`psite_info()`).  
- Save P-site info per sample (`psite_<sample>.rds`).

### Step 6: Summarize ribosome occupancy
- Map transcripts to genes using Entrez IDs.  
- Aggregate P-sites per gene.  
- Save per-sample gene counts (`gene_counts_<sample>.csv`).

### Step 7: Plot per-sample barplots
- Barplots of P-site counts per gene per sample.  
- Output saved to `plots/`.

### Step 8–9: Merge samples & plot group statistics
- Merge all per-sample gene counts.  
- Reshape for plotting across conditions (Control vs Treatment).  
- Plot mean ± SD barplots across groups (`grouped_gene_ribosome_occupancy.png`).

### Step 10: Statistical tests
- Perform per-gene t-tests between Control and Treatment.  
- Assign significance symbols: `***` (<0.001), `**` (<0.01), `*` (<0.05), `ns` (not significant).  
- Save results to `gene_ribosome_occupancy_stats.csv`.

---

## 📊 Outputs

- `bam_info.csv` → metadata for BAM files  
- `psite_<sample>.rds` → P-site information per sample  
- `gene_counts_<sample>.csv` → gene-level P-site counts  
- `merged_gene_data.csv` → combined counts for all samples  
- `plots/` → per-sample and grouped barplots  
- `gene_ribosome_occupancy_stats.csv` → t-test results  

---

## 🔍 Notes

- Change the `bam_dir` path to point to your local BAM files.  
- Adjust sample indices (`i`) when processing individual files.  
- Gene list can be modified as needed for your study.  
- Duplicate filtering parameters (`extremity`, `keep`) depend on your experimental protocol.  

---

## 📚 References

- **riboWaltz**:
  
  [https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/master/ReferenceManual.pdf](https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/master/ReferenceManual.pdf)


