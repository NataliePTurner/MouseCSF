# ============================================
# MOUSE CSF PROTEOMICS ANALYSIS
# ============================================

# PART 1: LOAD PACKAGES ----
library(diann)
library(dplyr)
library(tidyverse)
library(data.table)
library(readxl)
library(arrow)
library(mixOmics)
library(DEqMS)
library(ggrepel)
library(VennDiagram)
library(grid)
library(enrichR)

# PART 2: LOAD DATA ----
df <- read_parquet("report.parquet",
                   col_select = NULL,
                   as_data_frame = TRUE,
                   mmap = TRUE)

# Load annotations
df_2_ann <- read.csv("annotations.csv")

# Merge with annotations
raw <- merge(df, df_2_ann, by = 'Run')

# Remove contaminants and outlier
raw_clean <- raw %>%
  filter(!str_detect(Protein.Group, "cRAP"),
         !str_detect(Protein.Ids, "Biognosys\\|iRT-Kit_WR_fusion|P00761\\|TRYP_PIG"),
         Run != "CSF_M15_run2.raw")

setDT(raw_clean)

# PART 3: PROTEIN QUANTIFICATION ----
formatted_data <- diann_maxlfq(
  raw_clean %>% filter(Q.Value < 0.01 & PG.Q.Value < 0.01 & Global.Q.Value < 0.01),
  sample.header = "Run",
  group.header = "Protein.Names",
  id.header = "Precursor.Id",
  quantity.header = "Precursor.Normalised"
)

# Log2 transform
log_data <- log2(formatted_data)

# PART 4: ALIGN SAMPLE METADATA TO DATA ORDER ----
# Clean annotations - remove empty rows and M15
ann_clean <- df_2_ann %>%
  filter(Run != "" & Run != "CSF_M15_run2.raw")

# Get the column order from log_data and match annotations
data_order <- colnames(log_data)
ann_matched <- ann_clean[match(data_order, ann_clean$Run), ]

# Create properly ordered factors
condition <- as.factor(ann_matched$Condition)
sex <- as.factor(tolower(ann_matched$Other))

# Create animal labels (numbered within each condition)
ann_matched$AnimalNum <- ave(seq_len(nrow(ann_matched)), 
                             ann_matched$Condition, 
                             FUN = seq_along)
animal_labels <- paste0(ifelse(ann_matched$Condition == "Dep", "D", "ND"), 
                        ann_matched$AnimalNum)

# Verify alignment
cat("=== SAMPLE MAPPING ===\n")
print(data.frame(
  Run = ann_matched$Run,
  Condition = condition,
  Sex = sex,
  Label = animal_labels
))

# PART 5: PREPARE DATA FOR ANALYSIS ----
final_data <- as.data.frame(log_data)
colnames(final_data) <- animal_labels

# Separate by condition
dep_labels <- animal_labels[condition == "Dep"]
nondep_labels <- animal_labels[condition == "Non-dep"]

Dep <- final_data[, dep_labels]
Non_dep <- final_data[, nondep_labels]

# Filter: proteins in at least 2 replicates per group
Dep_filtered <- Dep[rowSums(!is.na(Dep)) >= 2, ]
Non_dep_filtered <- Non_dep[rowSums(!is.na(Non_dep)) >= 2, ]

# write.csv(Dep_filtered, file = "Dep_filtered_PG.csv")
# write.csv(Non_dep_filtered, file = "Non_dep_filtered_PG.csv")

cat("\nDep proteins (>=2 reps):", nrow(Dep_filtered), "\n")
cat("Non-dep proteins (>=2 reps):", nrow(Non_dep_filtered), "\n")

# PART 6: NORMALIZATION AND QC ----
log_data_complete <- na.omit(log_data)
norm_data <- normalizeMedianValues(log_data_complete)

# Boxplot QC
par(las = 2, lwd = 1, mar = c(5, 10, 4, 2))
boxplot(as.data.frame(norm_data), main = "Boxplot of Normalized Data", horizontal = TRUE)

# PART 7: PCA ANALYSIS ----
X <- t(log_data_complete)

# Match labels to complete data
complete_samples <- colnames(log_data_complete)
idx <- match(complete_samples, ann_matched$Run)

rownames(X) <- animal_labels[idx]
pca_condition <- condition[idx]
pca_sex <- sex[idx]

# Run PCA
pca.protein <- pca(X, ncomp = 5, center = TRUE, scale = TRUE)

# Extract scores
pca_scores <- as.data.frame(pca.protein$variates$X)
pca_scores$Condition <- pca_condition
pca_scores$Sex <- pca_sex
pca_scores$Animal <- rownames(X)

var_explained <- pca.protein$prop_expl_var$X * 100

# Verify before plotting
cat("\n=== PCA MAPPING ===\n")
print(pca_scores[, c("Animal", "Condition", "Sex")])

# PCA Plot
ggplot(pca_scores, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(aes(shape = Sex), size = 4) +
  geom_text_repel(aes(label = Animal), size = 5, color = "black") +
  stat_ellipse(aes(group = Condition), level = 0.95, linewidth = 1.2) +
  scale_color_manual(values = c("Dep" = "#388ECC", "Non-dep" = "#F68B33")) +
  scale_shape_manual(values = c("female" = 16, "male" = 17),
                     labels = c("Female", "Male")) +
  labs(title = "PCA on Mouse CSF",
       x = paste0("PC1 (", round(var_explained[1]), "%)"),
       y = paste0("PC2 (", round(var_explained[2]), "%)")) +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

# PART 8: VENN DIAGRAM ----
filtered_df_list <- list(
  Non_dep = Non_dep_filtered,
  Dep = Dep_filtered
)

# Venn diagram function
plot_venn <- function(filtered_list, set_names = NULL) {
  sets <- lapply(filtered_list, rownames)
  if (is.null(set_names)) set_names <- names(filtered_list)
  
  venn_colors <- c("#F68B33", "#388ECC")
  
  venn.plot <- venn.diagram(
    x = sets,
    category.names = set_names,
    filename = NULL,
    col = "transparent",
    fill = venn_colors[1:length(sets)],
    alpha = 0.3,
    cat.col = venn_colors[1:length(sets)],
    cat.cex = 1.5,
    cex = 1.5,
    fontface = "bold",
    cat.fontface = "bold",
    cat.default.pos = "outer",
    margin = 0.1
  )
  
  grid.newpage()
  grid.draw(venn.plot)
}

plot_venn(filtered_df_list, c("Non-Dep", "Dep"))

# PART 9: EXTRACT VENN ELEMENTS ----
protein_to_gene <- raw_clean %>%
  as.data.frame() %>%
  distinct(Protein.Names, Genes) %>%
  rename(protein = Protein.Names, gene = Genes)

extract_venn_elements <- function(filtered_list, lookup_table) {
  sets <- lapply(filtered_list, rownames)
  set_names <- names(filtered_list)
  
  venn_elements <- list()
  
  for (i in seq_along(sets)) {
    unique_proteins <- setdiff(sets[[i]], unlist(sets[-i]))
    unique_genes <- lookup_table$gene[match(unique_proteins, lookup_table$protein)]
    venn_elements[[paste0(set_names[i], "_unique")]] <- unique_genes
  }
  
  overlap_proteins <- Reduce(intersect, sets)
  venn_elements$overlap <- lookup_table$gene[match(overlap_proteins, lookup_table$protein)]
  
  return(venn_elements)
}

venn_elements <- extract_venn_elements(filtered_df_list, protein_to_gene)

cat("\n=== VENN ELEMENT COUNTS ===\n")
cat("Non-dep unique:", length(venn_elements$Non_dep_unique), "\n")
cat("Dep unique:", length(venn_elements$Dep_unique), "\n")
cat("Overlap:", length(venn_elements$overlap), "\n")

# Save venn elements
write.csv(venn_elements$Non_dep_unique, file = "Venn_elements_Non_dep_unique.csv")
write.csv(venn_elements$Dep_unique, file = "Venn_elements_Dep_unique.csv")
write.csv(venn_elements$overlap, file = "Venn_elements_overlap.csv")

# PART 10: ENRICHMENT ANALYSIS ----
setEnrichrSite("Enrichr")
dbs <- c("Mouse_Gene_Atlas", "GO_Biological_Process_2023", "GO_Cellular_Component_2023", 
         "GO_Molecular_Function_2023", "Jensen_COMPARTMENTS", "Jensen_DISEASES")

# Non-dep enrichment
enriched_nondep <- enrichr(venn_elements$Non_dep_unique, dbs)
plotEnrich(enriched_nondep[["Jensen_COMPARTMENTS"]], showTerms = 10,
           numChar = 80, y = "Count", orderBy = "P.value",
           title = "Proteins enriched in Mouse CSF (Non-dep)")

# Dep enrichment
enriched_dep <- enrichr(venn_elements$Dep_unique, dbs)
plotEnrich(enriched_dep[["Jensen_COMPARTMENTS"]], showTerms = 10,
           numChar = 80, y = "Count", orderBy = "P.value",
           title = "Proteins enriched in Mouse CSF (Dep)")

# Overlap enrichment
enriched_overlap <- enrichr(venn_elements$overlap, dbs)
plotEnrich(enriched_overlap[["Jensen_COMPARTMENTS"]], showTerms = 10,
           numChar = 80, y = "Count", orderBy = "P.value",
           title = "Proteins identified in Mouse CSF (Both groups)")

