#DIANN to MSstats import----
#Part i: Load required Packages----
library(diann)
library(dplyr)
library(tidyselect)
library(DEqMS)
library(ggbiplot)
library(data.table)
library(readxl)
library(arrow)
library(data.table)
library(mixOmics)

#Part ii: Setting the working directory (wd)----
#Either navigate to top banner and click on 
#'session' and 'setwd' then select folder,
#or as below, manually type in wd location (fill in ...)

#PART 1: Loading Data----
df <- read_parquet("report.parquet",
                   col_select = NULL,
                   as_data_frame = TRUE,
                   mmap = TRUE)

# Change working dir to MSstats output folder

#Annotate file
annotations <- "annotations.xlsx"
df_2_ann = read_xlsx(annotations, sheet = "annotations")

#Remove outliers
#EV runs selected/other runs removed--
df_2_merged<-merge(df, df_2_ann, by='Run')
raw <- df_2_merged

#remove cRAP and outlier
raw2 <- subset(raw, Protein.Ids!= "cRAP" & Protein.Ids!= "Biognosys|iRT-Kit_WR_fusion" & Run!="CSF_M15_run2.raw")

# Convert raw2 to a data.table
setDT(raw2)
raw3 <- raw2

# Extract clean gene and protein names from diann output
gene.names <- as.data.frame(cbind(raw3$Genes, raw3$Protein.Names))
colnames(gene.names) <- c("Genes","Protein")
gene.names.final <- gene.names %>% distinct(Protein, Genes)
CSF_protein.names <- gsub('_MOUSE', '', raw3$Protein.Names)
CSF_protein.names.final <- as.data.frame(raw3$protein.names)

formatted_data <- diann_maxlfq(raw3 %>% filter(Q.Value < 0.01 & PG.Q.Value < 0.01 & Global.Q.Value < 0.01),
                               sample.header = "Run",
                               group.header = "Protein.Names", # can use Genes column instead to concatenate to one gene per protein
                               id.header = "Precursor.Id",
                               quantity.header = "Precursor.Normalised")


# Use formatted data from diann_maxlfq, before na.omit
log_data <- log2(formatted_data)

# rename columns with Condition and BioReplicates
colnames(log_data) <- c("Dep", "Non-Dep","Non-Dep","Non-Dep","Dep","Dep","Non-Dep","Non-Dep","Dep")
final_data <- as.data.frame(log_data)


# Save to file
final_data <- final_data[, order(colnames(final_data))]

# Select individual groups for Venn
Non_dep <- final_data[,c(5:9)]
Dep <- final_data[,c(1:4)]

# Filter to include proteins identified in at least 2 replicates within in group
Non_dep_filtered <- Non_dep[rowSums(!is.na(Non_dep)) >= 2, ]
Dep_filtered <- Dep[rowSums(!is.na(Dep)) >= 2, ]

# Save to file
write.csv(Non_dep_filtered, file = "Non-dep_filtered_PG.csv")
write.csv(Dep_filtered, file = "Dep_filtered_PG.csv")

# QC Plots
# Apply normalization - equalizeMedians
final_data2 <- na.omit(log_data)
norm_data <- normalizeMedianValues(final_data2)

# Create boxplot of data after normalization
post_norm <- as.data.frame(norm_data)

# Increase left margin (3rd value) to fit labels
par(las = 2, lwd = 1, mar = c(5, 10, 4, 2)) # c(bottom, left, top, right)
boxplot(post_norm, main = "Boxplot of Normalized data", horizontal = TRUE)

# Extract groupings from annotations for PCA plot
rows_to_keep <- (df_2_ann$Run[c(1:2,4:10)])
groups <- as.factor(df_2_ann$Condition[c(1:2,4:10)])
sex <- as.factor(df_2_ann$Other[c(1:2,4:10)])

# Create nested list for PCA
pca.data <- list(phenotype = groups,
                 sex = sex,
                 names = rows_to_keep,
                 data = t(norm_data), #transpose columns/rows so samples are in rows and proteins are in columns
                 protein.name = raw3$Protein.Names)

#classify X and Y
X <- pca.data$data
Y <- pca.data$phenotype

# Set row names before performing PCA
rownames(X) <- pca.data$names
pca.protein <- pca(X, ncomp = 5, center = TRUE, scale = TRUE)

# Then use ind.names = TRUE
plotIndiv(pca.protein, ncomp = 5, ind.names = TRUE,
          group = Y, style = "ggplot2", pch = pca.data$sex,
          point.lwd = 1.5,
          ellipse = TRUE, size.xlabel = rel(1.5), size.ylabel = rel(1.5), cex = 4,
          size.axis = rel(1.5), size.legend.title = rel(1.5),
          legend = TRUE, title = 'PCA on Mouse CSF')

# OPTIONAL - generate PCA plot with file names - used to identify outliers

# library(ggplot2)
# 
# # Extract PCA coordinates
# pca_coords <- as.data.frame(pca.protein$x[,1:2])
# pca_coords$names <- pca.data$names
# pca_coords$group <- Y
# pca_coords$sex <- pca.data$sex
# 
# # Get variance explained
# var_exp <- pca.protein$prop_expl_var$X[1:2] * 100
# 
# ggplot(pca_coords, aes(x = PC1, y = PC2, color = group, shape = sex)) +
#   geom_point(size = 4) +
#   geom_text(aes(label = names), vjust = -0.8, hjust = 0.5, size = 3) +
#   xlab(paste0("PC1 (", round(var_exp[1], 1), "%)")) +
#   ylab(paste0("PC2 (", round(var_exp[2], 1), "%)")) +
#   ggtitle("PCA on Mouse CSF") +
#   theme_bw() +
#   theme(legend.position = "bottom")

# Generate Venn function
df_list <- list(
 Non_dep = Non_dep,
  Dep = Dep
)

filter_dataframes <- function(df_list) {
  # Define a helper function to filter a single dataframe
  filter_single_df <- function(df) {
    # Remove rows with more than one NA value
    df <- df[rowSums(!is.na(df)) >= 2, ]
    return(df)
  }
  
  # Apply the filter function to each dataframe in the list and retain names
  filtered_list <- setNames(lapply(df_list, filter_single_df), names(df_list))
  
  return(filtered_list)
}

filtered_df_list <- filter_dataframes(df_list)

library(VennDiagram)

plot_venn <- function(filtered_list) {
  # Create a list of sets based on the rownames of each filtered dataframe
  sets <- lapply(filtered_list, function(df) rownames(df))
  
  # Remove empty sets
  sets <- sets[sapply(sets, length) > 0]
  
  # Check if the number of sets is between 2 and 5
  if (length(sets) < 2 || length(sets) > 5) {
    stop("The number of sets for a Venn diagram should be between 2 and 5.")
  }
  
  # Generate a Venn diagram
  venn.plot <- venn.diagram(
    x = sets,
    category.names = paste("Set", seq_along(sets)),
    filename = NULL,
    output = TRUE
  )
  
  # Plot the Venn diagram
  grid.draw(venn.plot)
}

# Example usage:
plot_venn(filtered_df_list)

library(VennDiagram)
library(grid)

plot_selected_venn <- function(filtered_list, selected_indices) {
  # Create a list of sets based on the rownames of each filtered dataframe
  sets <- lapply(filtered_list, function(df) rownames(df))
  
  # Select only the sets specified by the indices
  selected_sets <- sets[selected_indices]
  
  # Get the first column name from each selected dataframe
  first_column_names <- sapply(filtered_list[selected_indices], function(df) colnames(df)[1])
  
  # Check if the number of selected sets is between 2 and 5
  if (length(selected_sets) < 2 || length(selected_sets) > 5) {
    stop("The number of selected sets for a Venn diagram should be between 2 and 5.")
  }
  
  # Generate category names with the first column names
  category_names <- sapply(seq_along(selected_sets), function(j) {
    first_column_names[j]
  })
  
  # Define colors for the Venn diagram
  venn_colors <- c("#F8766D", "#00BFC4", "#7CAE00", "#C77CFF", "#00BA38")
  
  # Generate a Venn diagram with custom colors and font
  venn.plot <- venn.diagram(
    x = selected_sets,
    category.names = category_names,
    filename = NULL,
    output = TRUE,
    col = "transparent",  # Border color
    fill = venn_colors[1:length(selected_sets)],  # Fill colors
    alpha = 0.2,  # Increased transparency
    cat.col = venn_colors[1:length(selected_sets)],  # Category text colors
    cat.cex = 1.5,  # Category text size
    fontfamily = "sans",  # Font family for text
    cex = 1.5,  # Size of numbers
    fontface = "bold",  # Bold text
    cat.fontface = "bold",  # Bold text for category names
    cat.default.pos = "outer",  # Default position for category names
    margin = 0.1  # Margin around the plot
  )
  
  # Convert the Venn diagram to a grob object
  venn_grob <- gTree(children = venn.plot)
  
  # Create a ggplot object with the Venn diagram
  p <- ggplot() +
    annotation_custom(venn_grob) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = "white"),
      plot.margin = margin(1, 1, 1, 1, "cm")
    )
  
  # Plot the Venn diagram
  print(p)
}

# Select the indices of the dataframes you want to compare, e.g., 1, 3, and 5
selected_indices <- c(1,2)
plot_selected_venn(filtered_df_list, selected_indices)


# Create lookup table
protein_to_gene <- data.frame(
  protein = raw3$Protein.Names,
  gene = raw3$Genes
)

extract_venn_elements <- function(filtered_list, selected_indices, lookup_table) {
  # Create a list of sets based on the rownames of each filtered dataframe
  sets <- lapply(filtered_list, function(df) rownames(df))
  
  # Select only the sets specified by the indices
  selected_sets <- sets[selected_indices]
  
  # Get the names of the selected dataframes
  selected_names <- names(filtered_list)[selected_indices]
  
  # Create a named list to store unique and overlapping elements
  venn_elements <- list()
  
  # Find unique and overlapping elements
  for (i in seq_along(selected_sets)) {
    set_name <- selected_names[i]
    set_elements <- selected_sets[[i]]
    unique_elements <- set_elements[!set_elements %in% unlist(selected_sets[-i])]
    unique_gene_names <- lookup_table$gene[match(unique_elements, lookup_table$protein)]
    venn_elements[[paste0(set_name, "_unique")]] <- unique_gene_names
  }
  
  # Find overlapping elements
  overlap_elements <- Reduce(intersect, selected_sets)
  overlap_gene_names <- lookup_table$gene[match(overlap_elements, lookup_table$protein)]
  venn_elements$overlap <- overlap_gene_names
  
  return(venn_elements)
}

# Select the indices of the dataframes you want to compare, e.g., 1, 3, and 5
selected_indices <- c(1,2)
venn_elements <- extract_venn_elements(filtered_df_list, selected_indices, protein_to_gene)

write.csv(venn_elements$Non_dep_unique, file = "Venn_elements_Non_dep_unique.csv")
write.csv(venn_elements$Dep_unique, file = "Venn_elements_Dep_unique.csv")
write.csv(venn_elements$overlap, file = "Venn_elements_overlap.csv")

# now read in results to EnrichR

library(enrichR)

websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}

if (websiteLive) dbs <- listEnrichrDbs()

# Select GO dbs and read in Gene names from dataframes

dbs <- c("Mouse_Gene_Atlas","GO_Biological_Process_2023","GO_Cellular_Component_2023","GO_Molecular_Function_2023","Jensen_COMPARTMENTS","Jensen_DISEASES") #name of dbs
if (websiteLive) {
  enriched <- enrichr(venn_elements$Non_dep_unique, dbs) #Selects column in dataframe containing gene names
}
Non_dep_unique_enrichR <- if (websiteLive) enriched[[5]] #view tabular results

# Create plot of top 20 enriched terms
if (websiteLive) {
  plotEnrich(enriched[[5]], showTerms = 10, #edited to show 10 terms
             numChar = 80, y = "Count", orderBy = "P.value",
             title = "Proteins enriched in Mouse CSF (Non-dep)") #adjust numChar to show full terms in plot
}

# Repeat for Vapor proteins
if (websiteLive) {
  enriched <- enrichr(venn_elements$Dep_unique, dbs) #Selects column in dataframe containing gene names
}

Dep_unique_enrichR <- if (websiteLive) enriched[[5]] #view tabular results

# Create plot of top 20 enriched terms
if (websiteLive) {
  plotEnrich(enriched[[5]], showTerms = 10, #edited to show 10 terms
             numChar = 80, y = "Count", orderBy = "P.value",
             title = "Proteins enriched in Mouse CSF (Dep)") #adjust numChar to show full terms in plot
}

# Repeat for overlap proteins
if (websiteLive) {
  enriched <- enrichr(venn_elements$overlap, dbs) #Selects column in dataframe containing gene names
}

overlap_enrichR <- if (websiteLive) enriched[["Mouse_Gene_Atlas"]] #view tabular results

# Create plot of top 20 enriched terms
if (websiteLive) {
  plotEnrich(enriched[[5]], showTerms = 10, #edited to show 10 terms
             numChar = 80, y = "Count", orderBy = "P.value",
             title = "Proteins identified in Mouse CSF (Non-dep and Dep)") #adjust numChar to show full terms in plot
}
