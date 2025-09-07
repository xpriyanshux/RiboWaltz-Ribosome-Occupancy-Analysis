############################################################
# Purpose: Process BAM files from RiboLace experiments, calculate P-site offsets,
#          summarise gene-specific ribosome occupancy, plot results, and perform basic statistics.
############################################################

# ----------------------------------------
# STEP 0: Load required libraries
# ----------------------------------------
library(fs)                    # File system operations (copy, delete, list files)
library(riboWaltz)             # Core package for ribosome profiling analysis
library(GenomicFeatures)       # Handling genomic features, transcripts, etc.
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # Human transcript annotation
library(org.Hs.eg.db)          # Gene symbol <-> Entrez ID mapping
library(AnnotationDbi)         # Interface for annotation databases
library(ggplot2)               # Plotting
library(dplyr)                 # Data manipulation
library(tidyr)                 # Reshaping data

# ----------------------------------------
# STEP 1: Define BAM directory and list files
# ----------------------------------------
bam_dir <- "ribolace/bam/file/path"  # Change to your local BAM directory
bam_files <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE)

# Define sample names (should match BAM files)
sample_names <- c("_01", "_02", "_03", "_04", "_05", "_06")

# Save info as a data.frame for easy loop processing
bam_info <- data.frame(
  sample = sample_names,
  path = bam_files,
  stringsAsFactors = FALSE
)

# Save to disk for later reference
write.csv(bam_info, "bam_info.csv", row.names = FALSE)
cat("✅ BAM info saved to bam_info.csv\n")

# ----------------------------------------
# STEP 2: Load annotation
# ----------------------------------------
# riboWaltz requires transcript annotation to map reads to P-sites
annot <- create_annotation(
  txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
  dataSource = NA,
  organism = NA
)

# ----------------------------------------
# STEP 3: Process individual BAM files
# ----------------------------------------
# Load saved BAM info table
bam_info <- read.csv("bam_info.csv", stringsAsFactors = FALSE)

# Loop or process one sample manually by setting index `i`
i <- 7  # Change this to 1,2,...,6 depending on sample
cat("Processing:", bam_info$sample[i], "\n")

# Clear temp folder to avoid contamination from previous runs
dir.create("temp_bam", showWarnings = FALSE)
fs::file_delete(dir("temp_bam", full.names = TRUE))

# Copy BAM file to temp folder for processing
fs::file_copy(bam_info$path[i], "temp_bam")

# ----------------------------------------
# STEP 4: Convert BAM to reads list and filter duplicates
# ----------------------------------------
reads_list <- bamtolist(bamfolder = "temp_bam", annotation = annot)

# Remove duplicate reads (can vary depending on sequencing protocol)
filtered <- duplicates_filter(
  reads_list,
  extremity = "both",       # Consider both 5' and 3' ends for duplicates
  keep = "shortest",        # Keep shortest read in case of duplicates
  output_class = "datatable"
)

# ----------------------------------------
# STEP 5: Calculate P-site offsets and annotate reads
# ----------------------------------------
p_offset <- psite(reads_list, flanking = 6)  # Flanking region for P-site estimation
reads_psite <- psite_info(reads_list, p_offset)

# Save P-site info for this sample
saveRDS(reads_psite[[1]], file = paste0("psite_", bam_info$sample[i], ".rds"))
cat("✅ Saved:", paste0("psite_", bam_info$sample[i], ".rds"), "\n")

# Clean memory
rm(reads_list, reads_psite, filtered)
gc(); gc()

# ----------------------------------------
# STEP 6: Summarize ribosome occupancy per gene
# ----------------------------------------
# Load sample P-site data
i <- 6  # Change per sample
sample_id <- bam_info$sample[i]
cat("🔍 Processing:", sample_id, "\n")
psite_data <- readRDS(paste0("psite_", sample_id, ".rds"))

# Define genes of interest
genes <- c("LDHA", "VEGFA", "SLC2A1", "AK4")  # GLUT1 = SLC2A1

# Map gene symbols to Entrez IDs
entrez_ids <- mapIds(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", column = "ENTREZID")

# Map transcripts to genes
tx_map <- select(
  TxDb.Hsapiens.UCSC.hg38.knownGene,
  keys = entrez_ids,
  keytype = "GENEID",
  columns = c("TXNAME", "GENEID")
)

# Merge P-site data with transcript-gene mapping
merged <- merge(psite_data, tx_map, by.x = "transcript", by.y = "TXNAME")

# Count P-sites per gene
gene_counts <- aggregate(psite ~ GENEID, data = merged, FUN = length)

# Map Entrez IDs back to gene symbols
gene_counts$SYMBOL <- mapIds(org.Hs.eg.db, keys = gene_counts$GENEID,
                             keytype = "ENTREZID", column = "SYMBOL")

# Clean and rename columns
gene_counts <- gene_counts[, c("SYMBOL", "psite")]
colnames(gene_counts)[2] <- sample_id

# Save gene counts
out_csv <- paste0("gene_counts_", sample_id, ".csv")
write.csv(gene_counts, out_csv, row.names = FALSE)
cat("✅ Saved counts to:", out_csv, "\n")

# ----------------------------------------
# STEP 7: Plot per-sample barplot
# ----------------------------------------
dir.create("plots", showWarnings = FALSE)

p <- ggplot(gene_counts, aes(x = SYMBOL, y = .data[[sample_id]])) +
  geom_bar(stat = "identity", fill = "#4f81bd") +
  theme_minimal(base_size = 14) +
  labs(title = paste("Ribosome Occupancy in", sample_id),
       x = "Gene", y = "P-site Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

out_plot <- paste0("plots/gene_counts_", sample_id, ".png")
ggsave(out_plot, plot = p, width = 6, height = 4)
cat("✅ Saved plot to:", out_plot, "\n")

# ----------------------------------------
# STEP 8: Merge all gene count CSVs for group analysis
# ----------------------------------------
samples <- c("_01", "_02", "_03", "_04", "_05", "_06")

gene_data_list <- lapply(samples, function(s) {
  df <- read.csv(paste0("gene_counts_", s, ".csv"))
  colnames(df)[2] <- s
  return(df)
})

merged_gene_data <- Reduce(function(x, y) merge(x, y, by = "SYMBOL"), gene_data_list)

# Reshape for plotting and stats
long_df <- pivot_longer(merged_gene_data,
                        cols = -SYMBOL,
                        names_to = "sample",
                        values_to = "psite_count")

# Assign experimental conditions
long_df$condition <- ifelse(long_df$sample %in% c("_01", "_02", "_03"), "Control", "Treatment")

# ----------------------------------------
# STEP 9: Plot mean ± SD barplot across groups
# ----------------------------------------
summary_stats <- long_df %>%
  group_by(SYMBOL, condition) %>%
  summarise(mean = mean(psite_count), sd = sd(psite_count), .groups = "drop")

ggplot(summary_stats, aes(x = SYMBOL, y = mean, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2,
                position = position_dodge(0.8)) +
  theme_minimal(base_size = 14) +
  labs(title = "Average Ribosome Occupancy (P-sites)",
       y = "P-site Count (Mean ± SD)", x = "Gene") +
  scale_fill_manual(values = c("Control" = "#4f81bd", "Treatment" = "#c0504d")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("plots/grouped_gene_ribosome_occupancy.png", width = 7, height = 5)
cat("✅ Saved grouped plot\n")

# ----------------------------------------
# STEP 10: Perform per-gene t-tests
# ----------------------------------------
results <- long_df %>%
  group_by(SYMBOL) %>%
  summarise(p_value = t.test(psite_count ~ condition)$p.value, .groups = "drop")

# Assign significance symbols
results$significance <- cut(results$p_value,
                            breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                            labels = c("***", "**", "*", "ns"))

print(results)

# Save statistics
write.csv(results, "gene_ribosome_occupancy_stats.csv", row.names = FALSE)
cat("✅ Saved statistical test results\n")
