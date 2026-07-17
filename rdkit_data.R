library(rdkitpyr)
library(dplyr)
library(readr)

# Load drug dataset
drug_data <- readRDS("data/drug_data_complete.rds")

cat("Calculating RDKit descriptors for", nrow(drug_data), "drugs...\n")

# Two accumulators: one for the curated 5-descriptor comparison,
# one for the full descriptor set
rdkit_curated <- data.frame()
rdkit_full <- data.frame()

for (i in 1:nrow(drug_data)) {
  drug <- drug_data$Drug[i]
  smiles <- drug_data$SMILES[i]
  
  if (is.na(smiles)) {
    rdkit_curated <- bind_rows(rdkit_curated,
                               data.frame(Drug = drug, MW_rdkit = NA, PSA_rdkit = NA,
                                          HBD_rdkit = NA, HBA_rdkit = NA, LogP_rdkit = NA))
    rdkit_full <- bind_rows(rdkit_full, data.frame(Drug = drug))
    next
  }
  
  tryCatch({
    desc <- rdkitpyr::CalculateAllDescriptors(smiles)
    
    # Curated subset for main dataset comparison
    rdkit_curated <- bind_rows(rdkit_curated,
                               data.frame(
                                 Drug = drug,
                                 MW_rdkit = desc$MolWt,
                                 PSA_rdkit = desc$TPSA,
                                 HBD_rdkit = desc$NumHDonors,
                                 HBA_rdkit = desc$NumHAcceptors,
                                 LogP_rdkit = desc$MolLogP
                               ))
    
    # Full descriptor set for standalone dataset
    desc$Drug <- drug
    rdkit_full <- bind_rows(rdkit_full, desc)
    
  }, error = function(e) {
    message(paste("Error for", drug, ":", e$message))
    rdkit_curated <<- bind_rows(rdkit_curated,
                                data.frame(Drug = drug, MW_rdkit = NA, PSA_rdkit = NA,
                                           HBD_rdkit = NA, HBA_rdkit = NA, LogP_rdkit = NA))
    rdkit_full <<- bind_rows(rdkit_full, data.frame(Drug = drug))
  })
  
  if (i %% 20 == 0) cat("Processed", i, "/", nrow(drug_data), "\n")
}

cat("\n RDKIT EXTRACTION SUMMARY \n")
cat("Drugs with RDKit descriptors:", sum(!is.na(rdkit_curated$MW_rdkit)), "/", nrow(rdkit_curated), "\n")

# Merge curated columns into main dataset
drug_data <- drug_data |>
  left_join(rdkit_curated, by = "Drug") |>
  mutate(
    MW_chembl = as.numeric(MW),
    PSA_chembl = as.numeric(PSA),
    MW_diff_rdkit = abs(MW_chembl - MW_rdkit),
    PSA_diff_rdkit = abs(PSA_chembl - PSA_rdkit),
    HBD_diff_rdkit = abs(HBD - HBD_rdkit),
    HBA_diff_rdkit = abs(HBA - HBA_rdkit),
    LogP_diff_rdkit = abs(LogP_chembl - LogP_rdkit)
  )

saveRDS(drug_data, "data/drug_data_complete.rds")
write_csv(drug_data, "data/drug_data_complete.csv")
cat("\nSaved drug_data_complete with", ncol(drug_data), "columns\n")

# Save full descriptor set as standalone dataset 
rdkit_full <- rdkit_full |> select(Drug, everything())

saveRDS(rdkit_full, "data/rdkit_all_descriptors.rds")
write_csv(rdkit_full, "data/rdkit_all_descriptors.csv")
cat("Saved rdkit_all_descriptors with", ncol(rdkit_full), "columns (", nrow(rdkit_full), "drugs )\n")

# Comparison summaries
cat("\n MW COMPARISON (ChEMBL vs RDKit)\n")
cat("Mean abs diff:", round(mean(drug_data$MW_diff_rdkit, na.rm = TRUE), 4), "\n")
cat("Max diff:", round(max(drug_data$MW_diff_rdkit, na.rm = TRUE), 4), "\n")

cat("\nPSA COMPARISON (ChEMBL vs RDKit)\n")
cat("Mean abs diff:", round(mean(drug_data$PSA_diff_rdkit, na.rm = TRUE), 4), "\n")
cat("Max diff:", round(max(drug_data$PSA_diff_rdkit, na.rm = TRUE), 4), "\n")
cat("Drugs with diff > 1:", sum(drug_data$PSA_diff_rdkit > 1, na.rm = TRUE), "\n")

cat("\nHBD/HBA COMPARISON (ChEMBL vs RDKit)\n")
cat("Drugs where HBD differs:", sum(drug_data$HBD_diff_rdkit != 0, na.rm = TRUE), "\n")
cat("Drugs where HBA differs:", sum(drug_data$HBA_diff_rdkit != 0, na.rm = TRUE), "\n")

cat("\nLogP COMPARISON (ChEMBL vs RDKit)\n")
cat("Mean abs diff:", round(mean(drug_data$LogP_diff_rdkit, na.rm = TRUE), 3), "\n")
cat("Drugs with diff > 1:", sum(drug_data$LogP_diff_rdkit > 1, na.rm = TRUE), "\n")