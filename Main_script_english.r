
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


# ===# === NMDS analysis for beaver ASV data ===

library(vegan)
library(ggplot2)
library(ggrepel)
library(dplyr)

# 1. Normalize ASV matrix (log-transformed already)
# If needed, use Hellinger transformation:
asv_hell <- decostand(asv_log, method = "hellinger")

# 2. Create NMDS (non-metric multidimensional scaling)
set.seed(123)  # reproducibility
nmds_result <- metaMDS(asv_hell, k = 2, trymax = 100, autotransform = FALSE)

# 3. Stressplot to assess model fit
pdf("nmds_stressplot.pdf")
stressplot(nmds_result)
dev.off()

# 4. Extract NMDS coordinates
nmds_coords <- as.data.frame(scores(nmds_result, display = "sites"))
nmds_coords$SampleID <- metadata_info$SampleID
nmds_coords$Population <- metadata_info$Area
nmds_coords$Sex <- metadata_info$Sex

# 5. Calculate distance from origin (optional, not used here)
nmds_coords$Distance <- sqrt(nmds_coords$NMDS1^2 + nmds_coords$NMDS2^2)

# 6. NMDS plot with labels for all samples
p_nmds <- ggplot(nmds_coords, aes(x = NMDS1, y = NMDS2, color = Sex, shape = Population)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text_repel(aes(label = SampleID),
                  size = 4, box.padding = 0.5, point.padding = 0.3, segment.alpha = 0.5) +
  scale_color_manual(values = c("F" = "#e41a1c", "M" = "#377eb8")) +
  scale_shape_manual(values = c("Sumava" = 16, "Sluknovsko" = 17)) +
  labs(title = "NMDS of Microbiota by sex and population",
       x = "NMDS1", y = "NMDS2") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("nmds_by_sex_population_labeled.pdf", plot = p_nmds, width = 8, height = 6)
