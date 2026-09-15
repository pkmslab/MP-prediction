"""
UnipKa macro-pKa predictions for the MP-prediction dataset.

Run standalone from the project root, in the environment where `unipka` is
installed (python 3.13 base env  NOT the `pka_prediction` env that holds
pkasolver, which is python 3.9):

    python unipka_predict.py

Input:  data/new_drug_dataset.csv   (needs the Drug and SMILES columns)
Output: data/unipka_macro_pka.csv

Raw predictions are kept alongside the filtered ones. Filtering matches the
pkasolver and MolGpKa convention: values outside pKa 3-11 are set to NA,
because outside that window the site is effectively fully ionised or fully
neutral at physiological pH and the prediction is not informative.

The output is written after every molecule, so an interrupted run keeps its
progress.
"""

import time

import numpy as np
import pandas as pd
import unipka

INPUT_FILE = "data/new_drug_dataset.csv"
OUTPUT_FILE = "data/unipka_macro_pka.csv"
BATCH_SIZE = 32
PKA_MIN, PKA_MAX = 3.0, 11.0


def main():
    df = pd.read_csv(INPUT_FILE)[["Drug", "SMILES"]].copy()
    print(f"Molecules: {len(df)}")

    model = unipka.UnipKa(batch_size=BATCH_SIZE)
    print("Uni-pKa model loaded.")

    df["unipka_acidic_raw"] = np.nan
    df["unipka_basic_raw"] = np.nan
    df["unipka_error"] = ""

    for i in range(len(df)):
        drug = df.loc[i, "Drug"]
        smiles = df.loc[i, "SMILES"]
        print(f"[{i + 1}/{len(df)}] {drug}")

        if pd.isna(smiles):
            df.loc[i, "unipka_error"] = "No SMILES"
            continue

        for mode, column, getter in (
            ("Acidic", "unipka_acidic_raw", model.get_acidic_macro_pka),
            ("Basic", "unipka_basic_raw", model.get_basic_macro_pka),
        ):
            try:
                start = time.time()
                value = getter(smiles)
                df.loc[i, column] = value
                print(f"  {mode} macro-pKa: {value}  ({time.time() - start:.1f}s)")
            except Exception as e:
                error = f"{mode}: {type(e).__name__}: {e}"
                existing = df.loc[i, "unipka_error"]
                df.loc[i, "unipka_error"] = f"{existing} | {error}" if existing else error
                print(f"  FAILED {error}")

        df.to_csv(OUTPUT_FILE, index=False)

    # Filtered columns . these are what 03_pka_comparison.Rmd reads
    for raw, filtered in (
        ("unipka_acidic_raw", "unipka_acidic_macro_pKa"),
        ("unipka_basic_raw", "unipka_basic_macro_pKa"),
    ):
        df[filtered] = df[raw].where(df[raw].between(PKA_MIN, PKA_MAX))

    df = df[
        [
            "Drug",
            "SMILES",
            "unipka_acidic_raw",
            "unipka_basic_raw",
            "unipka_acidic_macro_pKa",
            "unipka_basic_macro_pKa",
            "unipka_error",
        ]
    ]
    df.to_csv(OUTPUT_FILE, index=False)

    print("\nDONE ->", OUTPUT_FILE)
    print("Acidic (raw / in 3-11):",
          int(df["unipka_acidic_raw"].notna().sum()), "/",
          int(df["unipka_acidic_macro_pKa"].notna().sum()))
    print("Basic  (raw / in 3-11):",
          int(df["unipka_basic_raw"].notna().sum()), "/",
          int(df["unipka_basic_macro_pKa"].notna().sum()))
    print("Molecules with errors:", int((df["unipka_error"] != "").sum()))


if __name__ == "__main__":
    main()
