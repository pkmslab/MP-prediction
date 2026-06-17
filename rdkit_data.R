install.packages("rdkitpyr")
library(rdkitpyr)

ls("package:rdkitpyr")

GetPythonInfo()

# Test with caffeine SMILES
caffeine_smiles <- "CN1C=NC2=C1C(=O)N(C)C(=O)N2C"
result <- rdkitpyr::CalculateAllDescriptors(caffeine_smiles)
print(result)

library(dplyr)
library(readr)

# Load drug dataset
drug_data <- readRDS("data/drug_data_complete.rds")

# Calculate RDKit descriptors for all SMILES
cat("Calculating RDKit descriptors for", nrow(drug_data), "drugs...\n")

rdkit_results <- data.frame()

for (i in 1:nrow(drug_data)) {
  drug <- drug_data$Drug[i]
  smiles <- drug_data$SMILES[i]
  
  if (is.na(smiles)) {
    rdkit_results <- bind_rows(rdkit_results,
                               data.frame(Drug = drug, MW_rdkit = NA, PSA_rdkit = NA,
                                          HBD_rdkit = NA, HBA_rdkit = NA, LogP_rdkit = NA))
    next
  }
  
  result <- tryCatch({
    desc <- rdkitpyr::CalculateAllDescriptors(smiles)
    data.frame(
      Drug = drug,
      MW_rdkit = desc$MolWt,
      PSA_rdkit = desc$TPSA,
      HBD_rdkit = desc$NumHDonors,
      HBA_rdkit = desc$NumHAcceptors,
      LogP_rdkit = desc$MolLogP
    )
  }, error = function(e) {
    message(paste("Error for", drug, ":", e$message))
    data.frame(Drug = drug, MW_rdkit = NA, PSA_rdkit = NA,
               HBD_rdkit = NA, HBA_rdkit = NA, LogP_rdkit = NA)
  })
  
  rdkit_results <- bind_rows(rdkit_results, result)
  
  if (i %% 20 == 0) cat("Processed", i, "/", nrow(drug_data), "\n")
}

cat("\nRDKIT EXTRACTION SUMMARY\n")
cat("Drugs with RDKit descriptors:", sum(!is.na(rdkit_results$MW_rdkit)), "/", nrow(rdkit_results), "\n")

# Merge with drug_data
drug_data <- drug_data |>
  left_join(rdkit_results, by = "Drug")

# Save
saveRDS(drug_data, "data/drug_data_complete.rds")
write_csv(drug_data, "data/drug_data_complete.csv")

cat("\nSaved! Columns:", ncol(drug_data), "\n")
print(names(drug_data))

# Compare ChEMBL vs RDKit
drug_data <- drug_data |>
  mutate(
    MW_chembl = as.numeric(MW),
    PSA_chembl = as.numeric(PSA),
    
    MW_diff_rdkit = abs(MW_chembl - MW_rdkit),
    PSA_diff_rdkit = abs(PSA_chembl - PSA_rdkit),
    HBD_diff_rdkit = abs(HBD - HBD_rdkit),
    HBA_diff_rdkit = abs(HBA - HBA_rdkit),
    LogP_diff_rdkit = abs(LogP_chembl - LogP_rdkit)
  )

cat("MW COMPARISON (ChEMBL vs RDKit)\n")
cat("Mean abs diff:", round(mean(drug_data$MW_diff_rdkit, na.rm = TRUE), 4), "\n")
cat("Max diff:", round(max(drug_data$MW_diff_rdkit, na.rm = TRUE), 4), "\n")
cat("Drugs with diff > 1:", sum(drug_data$MW_diff_rdkit > 1, na.rm = TRUE), "\n\n")

cat("PSA COMPARISON (ChEMBL vs RDKit)\n")
cat("Mean abs diff:", round(mean(drug_data$PSA_diff_rdkit, na.rm = TRUE), 4), "\n")
cat("Max diff:", round(max(drug_data$PSA_diff_rdkit, na.rm = TRUE), 4), "\n")
cat("Drugs with diff > 1:", sum(drug_data$PSA_diff_rdkit > 1, na.rm = TRUE), "\n\n")

cat("HBD COMPARISON (ChEMBL vs RDKit)\n")
cat("Drugs with different values:", sum(drug_data$HBD_diff_rdkit != 0, na.rm = TRUE), "\n\n")

cat("HBA COMPARISON (ChEMBL vs RDKit)\n")
cat("Drugs with different values:", sum(drug_data$HBA_diff_rdkit != 0, na.rm = TRUE), "\n\n")

cat("LogP COMPARISON (ChEMBL vs RDKit)\n")
cat("Mean abs diff:", round(mean(drug_data$LogP_diff_rdkit, na.rm = TRUE), 3), "\n")
cat("Max diff:", round(max(drug_data$LogP_diff_rdkit, na.rm = TRUE), 3), "\n")
cat("Drugs with diff > 1:", sum(drug_data$LogP_diff_rdkit > 1, na.rm = TRUE), "\n\n")

cat("TOP 10 DRUGS WITH LARGEST LogP DISCREPANCY (RDKit)\n")
drug_data |>
  filter(LogP_diff_rdkit > 1) |>
  select(Drug, LogP_chembl, LogP_opera, LogP_rdkit, LogP_diff_rdkit) |>
  arrange(desc(LogP_diff_rdkit)) |>
  head(10) |>
  print()

cat("\nDRUGS WHERE HBD DIFFERS\n")
drug_data |>
  filter(HBD_diff_rdkit != 0) |>
  select(Drug, HBD, HBD_rdkit) |>
  print()

cat("\nDRUGS WHERE HBA DIFFERS\n")
drug_data |>
  filter(HBA_diff_rdkit != 0) |>
  select(Drug, HBA, HBA_rdkit) |>
  print()

# Save
saveRDS(drug_data, "data/drug_data_complete.rds")
write_csv(drug_data, "data/drug_data_complete.csv")
cat("\nSaved with", ncol(drug_data), "columns\n")

drug_data |>
  filter(PSA_diff_rdkit > 1) |>
  select(Drug, PSA_chembl, PSA_rdkit, PSA_diff_rdkit) |>
  arrange(desc(PSA_diff_rdkit)) |>
  print()

drug_data <- readRDS("data/drug_data_complete.rds")

# Export SMILES with drug names for pkasolver
drug_data |>
  dplyr::filter(!is.na(SMILES)) |>
  dplyr::select(Drug, SMILES) |>
  readr::write_csv("data/pkasolver_input.csv")

cat("Exported", sum(!is.na(drug_data$SMILES)), "drugs\n")