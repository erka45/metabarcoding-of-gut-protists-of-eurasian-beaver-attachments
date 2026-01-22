setwd("~/Documents/Documents/Magister/Thesis/Analysis")

# === Libraries ===
library(readxl)
library(dplyr)
library(tibble)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(vegan)

# === 1. Load the data ===

metadata <- read_excel("metadata.xlsx")
taxa <- read_excel("tax_counts_fasta.xlsx")

# === 2. Prepare the data ===

# Convert Name column to character (string)
metadata$Name <- as.character(metadata$Name)

# Identify all columns that contain sample abundances
abundance_cols <- grep("^EKK-", colnames(taxa), value = TRUE)

# Create a matching table for sample IDs
sample_mapping <- data.frame(
  SampleID = abundance_cols,
  Name = gsub("EKK-", "", abundance_cols),
  stringsAsFactors = FALSE
)

# Fix known mismatch: Z2 should be Ž2
sample_mapping$Name[sample_mapping$Name == "Z2"] <- "Ž2"

# Merge with metadata to match sample names and numbers
sample_mapping <- merge(sample_mapping, metadata[, c("Name", "Number", "Sex", "Area")], by = "Name")

# Transpose abundance data (samples as rows)
abundance_matrix <- taxa[, sample_mapping$SampleID] %>% t()
colnames(abundance_matrix) <- paste0("ASV", 1:ncol(abundance_matrix))
abundance_matrix <- as.data.frame(abundance_matrix)
abundance_matrix$SampleID <- sample_mapping$Number

# Merge abundance matrix with metadata
merged <- merge(abundance_matrix, metadata, by.x = "SampleID", by.y = "Number")
asv_data <- merged %>% select(starts_with("ASV"))
metadata_info <- merged %>% select(SampleID, Sex, Area)

# === 3. Check structure ===

cat("Number of males and females:\n")
print(table(metadata_info$Sex))

cat("\nNumber of samples per population:\n")
print(table(metadata_info$Area))

cat("\nCross-tab: population × sex:\n")
print(table(Population = metadata_info$Area, Sex = metadata_info$Sex))

# === 4. Normalize abundance data ===

asv_log <- log1p(asv_data)  # log(x + 1)

# === 5. PCA ===

res.pca <- PCA(asv_log, graph = FALSE)

# Simplify population names for plotting
metadata_info$Area <- recode(metadata_info$Area,
                             "Šluknovsko" = "Sluknovsko",
                             "Šumava" = "Sumava")

# PCA by sex
p1 <- fviz_pca_ind(res.pca,
                   geom.ind = "point",
                   col.ind = metadata_info$Sex,
                   palette = c("blue", "red"),
                   addEllipses = TRUE,
                   legend.title = "Sex")

ggsave("pca_by_sex.pdf", plot = p1, width = 6, height = 5)

# PCA by population
p2 <- fviz_pca_ind(res.pca,
                   geom.ind = "point",
                   col.ind = metadata_info$Area,
                   palette = "jco",
                   addEllipses = TRUE,
                   legend.title = "Population")

ggsave("pca_by_population.pdf", plot = p2, width = 6, height = 5)

# === 6. PERMANOVA ===

bray <- vegdist(asv_log, method = "bray")

sink("permanova_results.txt")
cat("PERMANOVA by Sex\n")
print(adonis2(bray ~ Sex, data = metadata_info))

cat("\nPERMANOVA by Population\n")
print(adonis2(bray ~ Area, data = metadata_info))
sink()