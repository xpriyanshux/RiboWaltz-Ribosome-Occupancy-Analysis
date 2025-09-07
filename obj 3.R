library(fs)

# Define BAM directory
bam_dir <- "D:/New folder/rome/rrrr/ROME ASSIGNMENT/bam 2"
bam_files <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE)
sample_names <- c("IMGN_B11_IB_01", "IMGN_B11_IB_02", "IMGN_B11_IB_03",
                  "IMGN_B11_IB_04", "IMGN_B11_IB_05", "IMGN_B11_IB_06")

# Save info as data.frame
bam_info <- data.frame(
  sample = sample_names,
  path = bam_files,
  stringsAsFactors = FALSE
)

# Save to disk for loop reference
write.csv(bam_info, "bam_info.csv", row.names = FALSE)


# Load required libraries
library(riboWaltz)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Load annotation
annot <- create_annotation(txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene", dataSource = NA, organism = NA)

# Load saved BAM info
bam_info <- read.csv("bam_info.csv", stringsAsFactors = FALSE)

# Set which index to process (set manually!)
i <- 7  # CHANGE THIS to 1, 2, 3, ..., 6

cat("Processing:", bam_info$sample[i], "\n")

# Step 1: Clear temp folder
dir.create("temp_bam", showWarnings = FALSE)
fs::file_delete(dir("temp_bam", full.names = TRUE))

# Step 2: Copy BAM
fs::file_copy(bam_info$path[i], "temp_bam")

# Step 3: Process
reads_list <- bamtolist(bamfolder = "temp_bam", annotation = annot)

filtered <- duplicates_filter(
  reads_list,
  extremity = "both",
  keep = "shortest",
  output_class = "datatable"
)

p_offset <- psite(reads_list, flanking = 6)
reads_psite <- psite_info(reads_list, p_offset)

# Step 4: Save P-site result
saveRDS(reads_psite[[1]], file = paste0("psite_", bam_info$sample[i], ".rds"))
cat("✅ Saved:", paste0("psite_", bam_info$sample[i], ".rds"), "\n")

rm(reads_list)
rm(reads_psite)
rm(filtered)
gc()
gc()

# Load libraries
library(riboWaltz)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)

# Load annotation
annot <- create_annotation(txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene", dataSource = NA, organism = NA)

# Load BAM info table
bam_info <- read.csv("bam_info.csv", stringsAsFactors = FALSE)

# Set index manually (change from 1 to 6)
i <- 6
sample_id <- bam_info$sample[i]
cat("🔍 Processing:", sample_id, "\n")

# Load P-site table
psite_data <- readRDS(paste0("psite_", sample_id, ".rds"))

# ---------------------------
# Gene symbol to transcript
# ---------------------------

genes <- c("LDHA", "VEGFA", "SLC2A1", "AK4")  # GLUT1 = SLC2A1
entrez_ids <- mapIds(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", column = "ENTREZID")

tx_map <- select(TxDb.Hsapiens.UCSC.hg38.knownGene,
                 keys = entrez_ids,
                 keytype = "GENEID",
                 columns = c("TXNAME", "GENEID"))

# ---------------------------
# Merge and count P-sites
# ---------------------------

merged <- merge(psite_data, tx_map, by.x = "transcript", by.y = "TXNAME")
gene_counts <- aggregate(psite ~ GENEID, data = merged, FUN = length)

# Add gene symbols
gene_counts$SYMBOL <- mapIds(org.Hs.eg.db,
                             keys = gene_counts$GENEID,
                             keytype = "ENTREZID",
                             column = "SYMBOL")

# Clean and rename
gene_counts <- gene_counts[, c("SYMBOL", "psite")]
colnames(gene_counts)[2] <- sample_id

# ---------------------------
# Save as CSV
# ---------------------------

out_csv <- paste0("gene_counts_", sample_id, ".csv")
write.csv(gene_counts, out_csv, row.names = FALSE)
cat("✅ Saved counts to:", out_csv, "\n")

# ---------------------------
# PLOT Barplot of gene counts
# ---------------------------

# Create plot folder
dir.create("plots", showWarnings = FALSE)

# Plot
p <- ggplot(gene_counts, aes(x = SYMBOL, y = .data[[sample_id]])) +
  geom_bar(stat = "identity", fill = "#4f81bd") +
  theme_minimal(base_size = 14) +
  labs(title = paste("Ribosome Occupancy in", sample_id),
       x = "Gene", y = "P-site Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save plot
out_plot <- paste0("plots/gene_counts_", sample_id, ".png")
ggsave(out_plot, plot = p, width = 6, height = 4)
cat("✅ Saved plot to:", out_plot, "\n")





# ----------------------------------------
# Load libraries
# ----------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)

# ----------------------------------------
# Step 1: Read all CSVs
# ----------------------------------------
samples <- c("IMGN_B11_IB_01", "IMGN_B11_IB_02", "IMGN_B11_IB_03",
             "IMGN_B11_IB_04", "IMGN_B11_IB_05", "IMGN_B11_IB_06")

gene_data_list <- lapply(samples, function(s) {
  df <- read.csv(paste0("gene_counts_", s, ".csv"))
  colnames(df)[2] <- s
  return(df)
})

# ----------------------------------------
# Step 2: Merge into single table
# ----------------------------------------
merged_gene_data <- Reduce(function(x, y) merge(x, y, by = "SYMBOL"), gene_data_list)

# ----------------------------------------
# Step 3: Reshape and annotate
# ----------------------------------------
long_df <- pivot_longer(merged_gene_data,
                        cols = -SYMBOL,
                        names_to = "sample",
                        values_to = "psite_count")

long_df$condition <- ifelse(long_df$sample %in% c("IMGN_B11_IB_01", "IMGN_B11_IB_02", "IMGN_B11_IB_03"), "Control", "Treatment")

# ----------------------------------------
# Step 4: Barplot (Mean ± SD per gene)
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
# Step 5: Statistical tests (t-test per gene)
# ----------------------------------------
results <- long_df %>%
  group_by(SYMBOL) %>%
  summarise(p_value = t.test(psite_count ~ condition)$p.value,
            .groups = "drop")

results$significance <- cut(results$p_value,
                            breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                            labels = c("***", "**", "*", "ns"))

print(results)
write.csv(results, "gene_ribosome_occupancy_stats.csv", row.names = FALSE)
cat("✅ Saved statistical test results\n")
