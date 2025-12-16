# Group-Specific Protein Detection Analysis in R
# Identifies proteins detected in ≥2 animals from one group and ≤1 animals from the opposite group

# Load required libraries
library(readr)
library(dplyr)
library(readxl)
library(tibble)
library(extrafont)

# Read the data
protein_data <- read_csv("final_data_outlier_rm.csv")
annotations <- read_excel("annotations_heatmap.xlsx", sheet = "annotations")

# Examine data structure
print("Data dimensions:")
print(dim(protein_data))
print("Column names:")
print(colnames(protein_data))

# Column mapping (after removing outlier at original position 6)
# Original positions were: Dep at 2,7,8,11 and Non-dep at 3,4,5,6,9,10
# After removing position 6, positions shift:
Dep_cols <- c(2, 6, 7, 10)  # 4 Dep samples
Non_dep_cols <- c(3, 4, 5, 8, 9)  # 5 Non-dep samples (outlier removed)

print(paste("Dep columns:", paste(Dep_cols, collapse = ", ")))
print(paste("Non-dep columns:", paste(Non_dep_cols, collapse = ", ")))

# Function to count non-NA detections
count_detections <- function(row, col_indices) {
  sum(!is.na(row[col_indices]))
}

# Create analysis dataframe
protein_analysis <- protein_data %>%
  mutate(
    protein = `...1`,
    dep_detections = apply(., 1, count_detections, Dep_cols),
    non_dep_detections = apply(., 1, count_detections, Non_dep_cols),
    total_detections = dep_detections + non_dep_detections
  ) %>%
  dplyr::select(protein, dep_detections, non_dep_detections, total_detections, everything())

# Identify group-specific proteins with improved ranking
dep_specific <- protein_analysis %>%
  dplyr::filter(dep_detections >= 2 & non_dep_detections <= 1) %>%
  mutate(
    dep_frequency = dep_detections / length(Dep_cols),
    non_dep_frequency = non_dep_detections / length(Non_dep_cols), 
    specificity_score = dep_frequency - non_dep_frequency,
    perfect_specificity = ifelse(non_dep_detections == 0, 1, 0)
  ) %>%
  arrange(desc(perfect_specificity), desc(specificity_score), desc(dep_detections)) %>%
  mutate(specificity_type = "dep_specific")

non_dep_specific <- protein_analysis %>%
  dplyr::filter(non_dep_detections >= 2 & dep_detections <= 1) %>%
  mutate(
    dep_frequency = dep_detections / length(Dep_cols),
    non_dep_frequency = non_dep_detections / length(Non_dep_cols),
    specificity_score = non_dep_frequency - dep_frequency,
    perfect_specificity = ifelse(dep_detections == 0, 1, 0)
  ) %>%
  arrange(desc(perfect_specificity), desc(specificity_score), desc(non_dep_detections)) %>%
  mutate(specificity_type = "non_dep_specific")

# Combine results
all_specific <- bind_rows(dep_specific, non_dep_specific)

# Print summary
cat("=== ANALYSIS SUMMARY ===\n")
cat("Total proteins analyzed:", nrow(protein_data), "\n")
cat("Dep-specific proteins:", nrow(dep_specific), "\n")
cat("Non-dep-specific proteins:", nrow(non_dep_specific), "\n")
cat("Total group-specific proteins:", nrow(all_specific), "\n\n")

# Display top results
cat("=== TOP 10 DEP-SPECIFIC PROTEINS ===\n")
dep_top10 <- head(dep_specific, 10)
for(i in 1:nrow(dep_top10)) {
  cat(sprintf("%2d. %-15s (%d/4 dep, %d/5 non-dep)\n", 
              i, dep_top10$protein[i], 
              dep_top10$dep_detections[i], 
              dep_top10$non_dep_detections[i]))
}

cat("\n=== TOP 10 NON-DEP-SPECIFIC PROTEINS ===\n")
non_dep_top10 <- head(non_dep_specific, 10)
for(i in 1:nrow(non_dep_top10)) {
  cat(sprintf("%2d. %-15s (%d/4 dep, %d/5 non-dep)\n", 
              i, non_dep_top10$protein[i], 
              non_dep_top10$dep_detections[i], 
              non_dep_top10$non_dep_detections[i]))
}

# Perfect group-specific proteins
perfect_dep <- dep_specific %>%
  dplyr::filter(dep_detections == 4 & non_dep_detections == 0)

perfect_non_dep <- non_dep_specific %>%
  dplyr::filter(dep_detections == 0 & non_dep_detections == 5)

cat("\n=== PERFECT DEP-SPECIFIC PROTEINS ===\n")
cat("(Detected in all 4 dep animals, 0 non-dep animals)\n")
cat("Count:", nrow(perfect_dep), "\n")
if(nrow(perfect_dep) > 0) {
  for(i in 1:min(10, nrow(perfect_dep))) {
    cat(sprintf("%2d. %s\n", i, perfect_dep$protein[i]))
  }
  if(nrow(perfect_dep) > 10) cat("... and", nrow(perfect_dep) - 10, "more\n")
}

cat("\n=== PERFECT NON-DEP-SPECIFIC PROTEINS ===\n")
cat("(Detected in 0 dep animals, all 5 non-dep animals)\n")
cat("Count:", nrow(perfect_non_dep), "\n")
if(nrow(perfect_non_dep) > 0) {
  for(i in 1:min(10, nrow(perfect_non_dep))) {
    cat(sprintf("%2d. %s\n", i, perfect_non_dep$protein[i]))
  }
  if(nrow(perfect_non_dep) > 10) cat("... and", nrow(perfect_non_dep) - 10, "more\n")
}

# # Save the basic results first
# write_csv(dep_specific, "dep_specific_proteins.csv")
# write_csv(non_dep_specific, "non_dep_specific_proteins.csv")
# write_csv(all_specific, "all_group_specific_proteins.csv")

cat("\n=== BASIC ANALYSIS COMPLETE ===\n")
cat("Files saved: dep_specific_proteins.csv, non_dep_specific_proteins.csv, all_group_specific_proteins.csv\n")

# Export key protein lists for pathway analysis
if(nrow(dep_specific) > 0) {
  top_dep_proteins <- head(dep_specific$protein, min(50, nrow(dep_specific)))
  writeLines(top_dep_proteins, "top_dep_specific_proteins.txt")
}

if(nrow(non_dep_specific) > 0) {
  top_non_dep_proteins <- head(non_dep_specific$protein, min(50, nrow(non_dep_specific)))
  writeLines(top_non_dep_proteins, "top_non_dep_specific_proteins.txt")
}

cat("Protein lists for pathway analysis saved.\n")

# ===== HEATMAP VISUALIZATION =====
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tibble)

cat("\n=== CREATING HEATMAP VISUALIZATION ===\n")

# Get top proteins from each group
top_dep <- head(dep_specific, 10)
top_non_dep <- head(non_dep_specific, 10)

# Combine proteins and maintain order
top_proteins <- c(top_dep$protein, top_non_dep$protein)

# Create detection matrix
detection_matrix <- protein_data %>%
  dplyr::filter(`...1` %in% top_proteins) %>%
  dplyr::select(1, all_of(c(Dep_cols, Non_dep_cols))) %>%
  column_to_rownames("...1") %>%
  mutate_all(~ifelse(is.na(.), 0, 1))

# Rename columns with correct animal numbers and sex
# Dep_cols c(2, 6, 7, 10) → M11, M19, M20, M8 
# From script: D1=M11(M), D2=M19(F), D3=M20(F), D4=M8(M)
# Non_dep_cols c(3, 4, 5, 8, 9) → M14, M16, M18, M3, M5
# From script: ND1=M14(F), ND2=M16(F), ND3=M18(F), ND4=M3(M), ND5=M5(M)
colnames(detection_matrix) <- c("D1-M", "D2-F", "D3-F", "D4-M", 
                                "ND1-F", "ND2-F", "ND3-F", "ND4-M", "ND5-M")

# Sort columns within each group by animal number
dep_cols <- grep("^D[0-9]", colnames(detection_matrix), value = TRUE)
nondep_cols <- grep("^ND", colnames(detection_matrix), value = TRUE)

dep_cols_sorted <- dep_cols[order(as.numeric(gsub("D([0-9]+)-[MF]", "\\1", dep_cols)))]
nondep_cols_sorted <- nondep_cols[order(as.numeric(gsub("ND([0-9]+)-[MF]", "\\1", nondep_cols)))]

detection_matrix <- detection_matrix[, c(dep_cols_sorted, nondep_cols_sorted)]

# Clean up protein names (remove _MOUSE suffix)
rownames(detection_matrix) <- gsub("_MOUSE", "", rownames(detection_matrix))

# Define protein name groups after cleaning
dep_protein_names <- gsub("_MOUSE", "", top_dep$protein)
nondep_protein_names <- gsub("_MOUSE", "", top_non_dep$protein)

# Sort rows by total detections (high to low) within each protein group
dep_rows <- detection_matrix[rownames(detection_matrix) %in% dep_protein_names, , drop = FALSE]
dep_rows <- dep_rows[order(rowSums(dep_rows), decreasing = TRUE), , drop = FALSE]

nondep_rows <- detection_matrix[rownames(detection_matrix) %in% nondep_protein_names, , drop = FALSE]
nondep_rows <- nondep_rows[order(rowSums(nondep_rows), decreasing = TRUE), , drop = FALSE]

# Recombine with Dep-specific on top
detection_matrix <- rbind(dep_rows, nondep_rows)

# Create column annotations with Group and Sex
# After sorting, column order is: D1-M, D2-F, D3-F, D4-M, ND1-F, ND2-F, ND3-F, ND4-M, ND5-M
col_annotations <- data.frame(
  Group = c(rep("Dep", 4), rep("Non-dep", 5)),
  Sex = c("M", "F", "F", "M", "F", "F", "F", "M", "M")
)
rownames(col_annotations) <- colnames(detection_matrix)

# Create row annotations to match sorted order
row_annotations <- data.frame(
  Specificity = c(rep("Dep-specific", nrow(dep_rows)), rep("Non-dep-specific", nrow(nondep_rows)))
)
rownames(row_annotations) <- rownames(detection_matrix)

# Define annotation colors
annotation_colors <- list(
  Group = c("Dep" = "#E31A1C", "Non-dep" = "#1F78B4"),
  Sex = c("M" = "#6A3D9A", "F" = "#FB9A99"),
  Specificity = c("Dep-specific" = "#FF7F00", "Non-dep-specific" = "#33A02C")
)

# Define heatmap colors
heatmap_colors <- colorRamp2(c(0, 1), c("white", "blue"))

# Create the heatmap
p2 <- Heatmap(as.matrix(detection_matrix),
              name = "Detection",
              col = heatmap_colors,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_row_names = TRUE,
              show_column_names = TRUE,
              row_names_gp = gpar(fontsize = 16),
              column_names_gp = gpar(fontsize = 16),
              heatmap_legend_param = list(title = NULL),
              show_heatmap_legend = FALSE,
              top_annotation = HeatmapAnnotation(
                df = col_annotations,
                col = annotation_colors,
                annotation_name_gp = gpar(fontsize = 14),
                annotation_legend_param = list(
                  title_gp = gpar(fontsize = 14),
                  labels_gp = gpar(fontsize = 12)
                )
              ),
              left_annotation = rowAnnotation(
                df = row_annotations,
                col = annotation_colors,
                annotation_name_gp = gpar(fontsize = 14),
                annotation_legend_param = list(
                  title_gp = gpar(fontsize = 14),
                  labels_gp = gpar(fontsize = 12)
                )
              ),
              row_split = factor(row_annotations$Specificity, 
                                 levels = c("Dep-specific", "Non-dep-specific")),
              row_gap = unit(5, "mm"),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.rect(x = x, y = y, width = width, height = height,
                          gp = gpar(col = "black", lwd = 0.5, fill = NA))
              })

# Draw the heatmap
draw(p2,
     annotation_legend_side = "right",
     padding = unit(c(2, 2, 2, 2), "cm"))

# Add main title
grid.text("Detection Patterns: Group-Specific Proteins",
          y = unit(1, "npc") - unit(1, "cm"),
          gp = gpar(fontsize = 20))

cat("\n=== HEATMAP VISUALIZATION COMPLETE ===\n")
