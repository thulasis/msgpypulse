#!/usr/bin/env python3
"""
peptide_crosstab_generator.py

Generates a peptide × sample intensity cross-tab from filtered PSM TSV files.

Features:
- Input: per-sample *_filtered.tsv files with columns:
    Peptide, Protein (or ProteinID), ParentIonIntensity
    Optional: PeptideFlanked
- Peptide intensities:
    • Raw intensity preserved per peptide per sample
    • Multiple PSMs for the same peptide in a sample are summed
    • Shared peptides are kept as-is
- Output: peptide × sample matrix with missing values filled as 0
- Suitable for metaproteomics, simple bacterial, and probe-based proteomics

Usage:
python peptide_crosstab_generator.py -i /path/to/filtered_psms/ -o peptide_crosstab.tsv
"""

import argparse
from pathlib import Path
import pandas as pd

def detect_protein_col(df):
    """Detect protein column (Protein or ProteinID)."""
    for col in ("Protein", "ProteinID"):
        if col in df.columns:
            return col
    raise ValueError(f"Protein column not found. Available columns: {list(df.columns)}")

def load_peptides(tsv_path):
    """Load a single filtered PSM file and clean columns."""
    df = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    protein_col = detect_protein_col(df)

    # Standardize and clean columns
    df["Peptide"] = df["Peptide"].astype(str).str.strip()
    if "PeptideFlanked" in df.columns:
        df["PeptideFlanked"] = df["PeptideFlanked"].astype(str).str.strip()
    else:
        df["PeptideFlanked"] = df["Peptide"]  # fallback
    df[protein_col] = df[protein_col].astype(str).str.strip()
    df["ParentIonIntensity"] = pd.to_numeric(df["ParentIonIntensity"], errors="coerce")

    df = df.dropna(subset=["Peptide", protein_col, "ParentIonIntensity"])
    df = df.rename(columns={protein_col: "Protein", "ParentIonIntensity": "Intensity"})

    
    df = df.groupby(["Peptide", "PeptideFlanked", "Protein"], as_index=False)["Intensity"].max() # changed from sum to max for biological relevance

    return df
    
def merge_peptide_ctab(tsvs):
    """Merge multiple sample PSM files into a peptide × sample matrix."""
    wide = None
    for tsv in tsvs:
        sample_name = tsv.stem.replace("_filtered", "")
        df = load_peptides(tsv)
        df = df.rename(columns={"Intensity": sample_name})

        if wide is None:
            wide = df
        else:
            wide = wide.merge(df, on=["Peptide", "PeptideFlanked", "Protein"], how="outer")

    # Fill missing values with 0
    sample_cols = [c for c in wide.columns if c not in ("Peptide", "PeptideFlanked", "Protein")]
    wide[sample_cols] = wide[sample_cols].fillna(0)

    # Sort for readability
    wide = wide.sort_values(["Protein", "PeptideFlanked", "Peptide"], kind="mergesort").reset_index(drop=True)

    return wide

def main():
    parser = argparse.ArgumentParser(
        description="Generate peptide × sample cross-tab from filtered PSM TSVs."
    )
    parser.add_argument("-i", "--input_dir", required=True,
                        help="Directory containing per-sample *_filtered.tsv files.")
    parser.add_argument("-o", "--output_tsv", required=True,
                        help="Path to write the peptide cross-tab TSV.")
    args = parser.parse_args()

    in_dir = Path(args.input_dir)
    tsvs = sorted([p for p in in_dir.iterdir() if p.is_file() and p.name.endswith("_filtered.tsv")])
    if not tsvs:
        raise FileNotFoundError(f"No *_filtered.tsv files found in {in_dir}")

    peptide_crosstab = merge_peptide_ctab(tsvs)

    # Save output
    out_path = Path(args.output_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    peptide_crosstab.to_csv(out_path, sep="\t", index=False)

    print(f"Peptide cross-tab written: {out_path} "
          f"(peptides: {len(peptide_crosstab)}, samples: {len(tsvs)})")

if __name__ == "__main__":
    main()