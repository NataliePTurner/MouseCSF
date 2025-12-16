# =============================================================================
# CSF Proteomics Analysis - Alcohol Dependence Study
# IL-6 Antibody Treatment in CIE-2BC Model
# =============================================================================

# Load libraries
library(enrichR)

# =============================================================================
# PART 1: Pathway Enrichment Analysis
# =============================================================================

# --- Dependent-Unique Genes ---
genes_dep <- c("KRT5", "LMAN2", "LYNX1", "CMBL", "BCAM", "CHI3L1", "IGSF8", "GPD1", 
               "KRT10", "DCN", "ACTA1", "CBLN3", "CTSA", "GKN3", "MMP2", "LAMB2", 
               "GFAP", "ACYP2", "ACHE", "LAG3", "PRRT3", "FAM19A5", "NFASC", 
               "ACO1", "XXYLT1", "RELN", "CADM3", "PSMA3", "PTPRS", "CSPG5", 
               "NBL1", "CALR", "CLTC", "LAMP1", "DLST", "EFEMP2", "NPTX1", 
               "SPON1", "PCDHGC5", "SHISA6", "RNPEP", "KRT17", "CAV1", "KRT6A", 
               "C1QA", "ADGRL1", "C6", "PTPRD", "GGH", "FSTL1", "PTPRN2", "SGCE", 
               "EML2", "SPOCK2", "ATP1A1", "UCHL1", "AFP", "TGFB2", "PLXNB2", 
               "CADM1", "COL6A1", "ITIH5", "LAMC1", "CDH5", "FAM20C", "ACAN", 
               "LFNG", "NXPH1", "LINGO1", "KRT16", "HEXA", "FBLN5", "HSPA5", 
               "HSP90B1", "ARSB", "YIPF3", "GPC1", "FHL1", "MAPK3", "ANXA5", 
               "IGFBP4", "PCOLCE", "NUCB1", "IDH1", "PI16", "SERPINA7", "HEG1", 
               "TXNL1", "CNTN4", "VCAM1", "YWHAQ", "YWHAB", "UBA1", "PLEKHB1", 
               "AHCY", "WDR1", "PFKM", "RNH1", "RAB5C", "CX3CL1", "FRRS1L", 
               "USP5", "TGM2", "PTMS", "VCP", "CGREF1", "HTRA1", "PVALB", "TGFBI", 
               "LTBP4", "VCL", "LGMN", "GDA", "GSTM2", "P4HB", "PDIA3", "VSTM2A", 
               "ATP5F1B", "THBS1", "CADM2", "HAGH", "HSP90AA1", "SORD")

# Run enrichment - Dependent
results_dep <- enrichr(genes_dep, 
                       databases = c("GO_Biological_Process_2023",
                                     "KEGG_2021_Human",
                                     "Reactome_2022",
                                     "WikiPathway_2023_Human"))

# View top results - Dependent
head(results_dep$GO_Biological_Process_2023[, c("Term", "Adjusted.P.value", "Genes")], 25)
head(results_dep$KEGG_2021_Human[, c("Term", "Adjusted.P.value", "Genes")], 15)
head(results_dep$Reactome_2022[, c("Term", "Adjusted.P.value", "Genes")], 15)

# Filter for IL-6/inflammation terms - Dependent
bp_dep <- results_dep$GO_Biological_Process_2023
bp_dep[grep("IL-6|interleukin|cytokine|inflam|acute phase|JAK|STAT", bp_dep$Term, ignore.case = TRUE), 
       c("Term", "Adjusted.P.value", "Genes")]


# --- Non-Dependent Unique Genes ---
genes_nondep <- c("MCAM", "SUMO2", "MINPP1", "F9", "CALB1", "SMPDL3A", "CTSH", 
                  "A1BG", "DPYSL2", "CFHR1", "CPN1", "EPB42", "CALB2", "HSPE1", 
                  "ANP32A", "TAGLN3", "CSRP1", "SAA4", "CDH15", "LASP1", "SNCB", 
                  "CNTNAP2", "PRRT2", "PRL", "SPR", "DDAH2", "ABHD14B", "ACP1", 
                  "CFHR4", "F13A1", "IGFALS", "CRP", "PROC", "APOB", "CELA1", 
                  "PROZ", "ADIPOQ", "SPTB", "ANK1", "ESD", "APEH")

# Run enrichment - Non-Dependent
results_nondep <- enrichr(genes_nondep, 
                          databases = c("GO_Biological_Process_2023",
                                        "KEGG_2021_Human",
                                        "Reactome_2022"))

# View top results - Non-Dependent
head(results_nondep$GO_Biological_Process_2023[, c("Term", "Adjusted.P.value", "Genes")], 25)
head(results_nondep$KEGG_2021_Human[, c("Term", "Adjusted.P.value", "Genes")], 15)
head(results_nondep$Reactome_2022[, c("Term", "Adjusted.P.value", "Genes")], 15)

# Filter for IL-6/inflammation terms - Non-Dependent
bp_nondep <- results_nondep$GO_Biological_Process_2023
bp_nondep[grep("IL-6|interleukin|cytokine|inflam|acute phase|JAK|STAT", bp_nondep$Term, ignore.case = TRUE), 
          c("Term", "Adjusted.P.value", "Genes")]


# =============================================================================
# PART 2: Sample Data Setup
# =============================================================================

# Read protein abundance data
prot <- read.csv("final_data_outlier_rm.csv", row.names = 1)

# Transpose so samples are rows
prot_t <- as.data.frame(t(prot))

# Extract mouse IDs from column names
prot_t$mouse_id <- gsub("CSF_M|_run2.raw|.raw", "", rownames(prot_t))

# Create sample metadata
sample_data <- data.frame(
  sample = rownames(prot_t),
  mouse_id = prot_t$mouse_id
)

# Add group
sample_data$group <- ifelse(sample_data$mouse_id %in% c("8", "11", "19", "20"), "Dep", "NonDep")

# Add sex
sample_data$sex <- NA
sample_data$sex[sample_data$mouse_id %in% c("3", "5", "8", "11")] <- "M"
sample_data$sex[sample_data$mouse_id %in% c("14", "16", "18", "19", "20")] <- "F"

# Add study labels (D1-D4, ND1-ND5) in numerical order of mouse number
sample_data$study_label <- NA
sample_data$study_label[sample_data$mouse_id == "8"] <- "D1"
sample_data$study_label[sample_data$mouse_id == "11"] <- "D2"
sample_data$study_label[sample_data$mouse_id == "19"] <- "D3"
sample_data$study_label[sample_data$mouse_id == "20"] <- "D4"
sample_data$study_label[sample_data$mouse_id == "3"] <- "ND1"
sample_data$study_label[sample_data$mouse_id == "5"] <- "ND2"
sample_data$study_label[sample_data$mouse_id == "14"] <- "ND3"
sample_data$study_label[sample_data$mouse_id == "16"] <- "ND4"
sample_data$study_label[sample_data$mouse_id == "18"] <- "ND5"

# Add BAC for dependent animals (final timepoint - 7/28)
sample_data$bac_final <- NA
sample_data$bac_final[sample_data$mouse_id == "8"] <- 222.076
sample_data$bac_final[sample_data$mouse_id == "11"] <- 248.154
sample_data$bac_final[sample_data$mouse_id == "19"] <- 257.66
sample_data$bac_final[sample_data$mouse_id == "20"] <- 140.963

# Add drinking behavior data
sample_data$baseline_drinking <- NA
sample_data$cie_drinking <- NA

# Dependent animals
sample_data$baseline_drinking[sample_data$mouse_id == "8"] <- 2.66
sample_data$baseline_drinking[sample_data$mouse_id == "11"] <- 2.72
sample_data$baseline_drinking[sample_data$mouse_id == "19"] <- 14.34
sample_data$baseline_drinking[sample_data$mouse_id == "20"] <- 13.21

sample_data$cie_drinking[sample_data$mouse_id == "8"] <- 14.12
sample_data$cie_drinking[sample_data$mouse_id == "11"] <- 8.53
sample_data$cie_drinking[sample_data$mouse_id == "19"] <- 24.30
sample_data$cie_drinking[sample_data$mouse_id == "20"] <- 17.38

# Non-dependent controls
sample_data$baseline_drinking[sample_data$mouse_id == "3"] <- 4.69
sample_data$baseline_drinking[sample_data$mouse_id == "5"] <- 6.37
sample_data$baseline_drinking[sample_data$mouse_id == "14"] <- 16.18
sample_data$baseline_drinking[sample_data$mouse_id == "16"] <- 6.13
sample_data$baseline_drinking[sample_data$mouse_id == "18"] <- 17.18

sample_data$cie_drinking[sample_data$mouse_id == "3"] <- 11.74
sample_data$cie_drinking[sample_data$mouse_id == "5"] <- 10.90
sample_data$cie_drinking[sample_data$mouse_id == "14"] <- 16.29
sample_data$cie_drinking[sample_data$mouse_id == "16"] <- 14.46
sample_data$cie_drinking[sample_data$mouse_id == "18"] <- 14.70

# Calculate escalation ratio
sample_data$escalation_ratio <- sample_data$cie_drinking / sample_data$baseline_drinking


# =============================================================================
# PART 3: Add Proteins of Interest
# =============================================================================

# BBB/Endothelial markers
sample_data$Vcam1 <- as.numeric(prot_t$VCAM1_MOUSE)
sample_data$Mmp2 <- as.numeric(prot_t$MMP2_MOUSE)
sample_data$Cdh5 <- as.numeric(prot_t$CADH5_MOUSE)
sample_data$Lamc1 <- as.numeric(prot_t$LAMC1_MOUSE)
sample_data$Lamb2 <- as.numeric(prot_t$LAMB2_MOUSE)

# Complement
sample_data$C1qa <- as.numeric(prot_t$C1QA_MOUSE)

# Astrogliosis
sample_data$Gfap <- as.numeric(prot_t$GFAP_MOUSE)
sample_data$Chi3l1 <- as.numeric(prot_t$CH3L1_MOUSE)

# ER stress
sample_data$Hspa5 <- as.numeric(prot_t$BIP_MOUSE)
sample_data$Hsp90b1 <- as.numeric(prot_t$ENPL_MOUSE)
sample_data$Calr <- as.numeric(prot_t$CALR_MOUSE)
sample_data$P4hb <- as.numeric(prot_t$PDIA1_MOUSE)
sample_data$Pdia3 <- as.numeric(prot_t$PDIA3_MOUSE)


# =============================================================================
# PART 4: Correlation Analyses
# =============================================================================

# Subset dependent animals
dep_data <- sample_data[sample_data$group == "Dep", ]

# View dependent animal summary
cat("\n--- Dependent Animal Summary ---\n")
print(dep_data[, c("study_label", "mouse_id", "sex", "bac_final", "baseline_drinking", 
                   "cie_drinking", "escalation_ratio")])

# Check available protein data for dependent animals
cat("\n--- Protein Data Availability (Dependent) ---\n")
print(dep_data[, c("study_label", "Gfap", "Mmp2", "Lamb2", "Lamc1", "Chi3l1")])

# --- Correlation with BAC ---
markers <- c("Gfap", "Mmp2", "Vcam1", "Cdh5", "Chi3l1", "Hspa5", "Hsp90b1", "Calr", "Lamc1", "C1qa")

cor_bac <- data.frame(protein = markers, rho = NA, pval = NA)

for (i in seq_along(markers)) {
  y <- dep_data[[markers[i]]]
  if (sum(!is.na(y)) >= 3) {
    test <- cor.test(dep_data$bac_final, y, method = "spearman", exact = FALSE)
    cor_bac$rho[i] <- round(test$estimate, 3)
    cor_bac$pval[i] <- round(test$p.value, 3)
  }
}

cat("\n--- Correlations with Final BAC ---\n")
print(cor_bac)


# --- Correlation with Drinking Behavior ---
markers_drinking <- c("Gfap", "Mmp2", "Lamc1", "Lamb2", "Hspa5", "Calr", "Chi3l1")

cor_drinking <- data.frame(
  protein = markers_drinking, 
  rho_cie = NA, p_cie = NA, 
  rho_escalation = NA, p_escalation = NA
)

for (i in seq_along(markers_drinking)) {
  y <- dep_data[[markers_drinking[i]]]
  if (sum(!is.na(y)) >= 3) {
    # CIE drinking correlation
    test1 <- cor.test(dep_data$cie_drinking, y, method = "spearman", exact = FALSE)
    cor_drinking$rho_cie[i] <- round(test1$estimate, 3)
    cor_drinking$p_cie[i] <- round(test1$p.value, 3)
    
    # Escalation ratio correlation
    test2 <- cor.test(dep_data$escalation_ratio, y, method = "spearman", exact = FALSE)
    cor_drinking$rho_escalation[i] <- round(test2$estimate, 3)
    cor_drinking$p_escalation[i] <- round(test2$p.value, 3)
  }
}

cat("\n--- Correlations with Drinking Behavior ---\n")
print(cor_drinking)


# =============================================================================
# PART 5: Visualizations
# =============================================================================

# --- Boxplots: Dependent vs Non-Dependent ---
par(mfrow = c(2, 3))

for (marker in c("Gfap", "Mmp2", "Lamb2", "Lamc1", "Hspa5", "Chi3l1")) {
  dep_vals <- sample_data[[marker]][sample_data$group == "Dep"]
  nondep_vals <- sample_data[[marker]][sample_data$group == "NonDep"]
  
  boxplot(list(Dep = dep_vals, NonDep = nondep_vals),
          main = marker,
          ylab = "log2 intensity",
          xlab = "Group",
          col = c("lightblue", "orange"))
  
  stripchart(list(Dep = dep_vals, NonDep = nondep_vals),
             vertical = TRUE, method = "jitter", add = TRUE, pch = 19, cex = 1.2)
}


# --- Scatterplots: Correlations with all 4 animals where possible ---
par(mfrow = c(2, 2))

# Mmp2 vs CIE drinking (all 4 animals have Mmp2 data)
plot(dep_data$cie_drinking, dep_data$Mmp2,
     pch = 19, cex = 1.5,
     col = ifelse(dep_data$sex == "M", "blue", "red"),
     xlab = "CIE Drinking (g/kg/day)", ylab = "Mmp2 (log2)",
     main = "MMP2 vs CIE Drinking",
     xlim = c(5, 28), ylim = c(14.8, 16.2))
text(dep_data$cie_drinking, dep_data$Mmp2, labels = dep_data$study_label, pos = 3)
legend("topleft", legend = c("Male", "Female"), pch = 19, col = c("blue", "red"))

# Mmp2 vs Escalation (all 4 animals)
plot(dep_data$escalation_ratio, dep_data$Mmp2,
     pch = 19, cex = 1.5,
     col = ifelse(dep_data$sex == "M", "blue", "red"),
     xlab = "Escalation Ratio (CIE/Baseline)", ylab = "Mmp2 (log2)",
     main = "MMP2 vs Escalation",
     xlim = c(0.5, 6), ylim = c(14.8, 16.2))
text(dep_data$escalation_ratio, dep_data$Mmp2, labels = dep_data$study_label, pos = 3)
legend("topleft", legend = c("Male", "Female"), pch = 19, col = c("blue", "red"))

# Lamb2 vs CIE drinking (D1 is NA)
plot(dep_data$cie_drinking, dep_data$Lamb2,
     pch = 19, cex = 1.5,
     col = ifelse(dep_data$sex == "M", "blue", "red"),
     xlab = "CIE Drinking (g/kg/day)", ylab = "Lamb2 (log2)",
     main = "Lamb2 vs CIE Drinking (ρ = 1.0)\n(D1 not detected)",
     xlim = c(5, 28), ylim = c(15.1, 16.3))
text(dep_data$cie_drinking, dep_data$Lamb2, labels = dep_data$study_label, pos = 3)
legend("topleft", legend = c("Male", "Female"), pch = 19, col = c("blue", "red"))

# GFAP vs Escalation (D1 is NA)
plot(dep_data$escalation_ratio, dep_data$Gfap,
     pch = 19, cex = 1.5,
     col = ifelse(dep_data$sex == "M", "blue", "red"),
     xlab = "Escalation Ratio (CIE/Baseline)", ylab = "GFAP (log2)",
     main = "GFAP vs Escalation (ρ = 1.0)\n(D1 not detected)",
     xlim = c(0.5, 6), ylim = c(16.3, 17.5))
text(dep_data$escalation_ratio, dep_data$Gfap, labels = dep_data$study_label, pos = 3)
legend("topleft", legend = c("Male", "Female"), pch = 19, col = c("blue", "red"))

# =============================================================================
# PART 6: Summary Tables
# =============================================================================

# Sex differences in drinking
cat("\n--- Sex Differences in Drinking (Dependent Animals) ---\n")
sex_summary <- aggregate(cbind(baseline_drinking, cie_drinking, escalation_ratio, bac_final) ~ sex, 
                         data = dep_data, FUN = mean, na.rm = TRUE)
print(sex_summary)

# Full sample data summary
cat("\n--- Full Sample Data ---\n")
print(sample_data[order(sample_data$group, sample_data$mouse_id), 
                  c("study_label", "mouse_id", "group", "sex", "bac_final", 
                    "cie_drinking", "escalation_ratio")])