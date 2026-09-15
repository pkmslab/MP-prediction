# MP-prediction : dataset extraction

Data extraction for the milk to plasma (M/P) ratio prediction work.

## Structure

```
MP-prediction/
├── MP-prediction.Rproj
├── README.md
├── 01_dataset_properties.Rmd     ChEMBL + PubChem + OPERA + RDKit  → main dataset
├── 02_pkasolver_logd.Rmd         pkasolver pKa, LogD7.4, fraction neutral → main dataset
├── 03_pka_comparison.Rmd         OPERA + MolGpKa + Uni-pKa → comparison dataset
├── pka_classification.py         SMARTS acid/base site classifier (used by 02)
├── unipka_predict.py             Uni-pKa macro-pKa predictions (used by 03)
└── data/
```

## Run order

```
01_dataset_properties.Rmd
   └─ pause at block 9a: upload data/new_opera_input.csv to EPA CompTox,
      save the result as data/new_opera_output.csv, then continue at 9b
*** restart R here ***
02_pkasolver_logd.Rmd            (conda env: pka_prediction)
unipka_predict.py                (conda env: the one with unipka installed)
03_pka_comparison.Rmd
```

**The restart between 01 and 02 is required.** Notebook 01 loads
`rdkitpyr`, which initialises reticulate against its own Python; reticulate
cannot switch environments afterwards, so `use_condaenv("pka_prediction")` in
02 fails with "failed to initialize requested version of Python". Each notebook
hands off through `data/new_drug_dataset.rds`, so a restart costs nothing. 02
checks for this and stops with a readable message.

Two different Python environments are involved and they cannot be merged:
pkasolver lives in `pka_prediction` (python 3.9), Uni-pKa needs python 3.13.
That is why Uni-pKa is a standalone script rather than a reticulate chunk.

## Outputs

**`data/new_drug_dataset.rds` / `.csv` — the main dataset (194 drugs)**

| Group | Columns |
|---|---|
| Identity | `Drug`, `ChEMBL_ID`, `SMILES`, `Source` |
| Properties (ChEMBL/PubChem) | `MW`, `PSA`, `HBD`, `HBA`, `LogP` + numeric `MW_chembl`, `PSA_chembl`, `LogP_chembl` |
| Fraction unbound | `Fu_median`, `N_fu`, `Fu_values`, `PPB_value`, `Fu_from_PPB`, `N_ppb`, `Fu_source`, `Fu_final` |
| Experimental (ChEMBL) | `pKa_exp`, `N_pka`, `pKa_values`, `LogD74_exp`, `N_logd`, `LogD_values` |
| OPERA | `pKa_acid_pred`, `pKa_base_pred`, `LogD74_pred`, `LogP_opera`, `Type` |
| RDKit | `MW_rdkit`, `PSA_rdkit`, `HBD_rdkit`, `HBA_rdkit`, `LogP_rdkit` + `*_diff_rdkit` |
| pkasolver | `pKa_pkasolver`, `N_pka_pkasolver`, `Type_pkasolver`, `pKa_acid_v1`, `pKa_base_v1`, `pKa_acid_v2`, `pKa_base_v2`, `FG_acid_v2`, `FG_base_v2` |
| Calculated from v2 | `f_neutral_v2`, `fraction_ionized_v2`, `LogD74_calc_pkasolver_v2` |

**`data/rdkit_all_descriptors.rds` / `.csv`** - the full RDKit descriptor table,
one row per drug, kept separate because of its width.

**`data/pka_comparison.rds` / `.csv`** - acid and base pKa per drug from
pkasolver v2, OPERA, MolGpKa and Uni-pKa side by side, plus `pKa_exp`. Nothing
from this file feeds back into the main dataset.

