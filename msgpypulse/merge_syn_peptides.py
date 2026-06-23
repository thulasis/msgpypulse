#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import sys
import pandas as pd

"""
Author: Harrison Hall, (PhD)
"""

def is_decoy(protein):
    return isinstance(protein, str) and protein.startswith("XXX_")

def load_safe_tsv(path: Path, sep="\t") :
    try:
        return pd.read_csv(path, sep=sep)
    except Exception as e:
        raise RuntimeError(f"Failed reading {path}: {e}")

def pick_template_rows(fdr) :
    """Pick one template row per peptide (lowest QValue if available)."""
    tmp = fdr.copy()
    if "QValue" in tmp.columns:
        tmp["__Q__"] = pd.to_numeric(tmp["QValue"], errors="coerce")
        tmp = tmp.sort_values(["Peptide", "__Q__"], kind="mergesort")
    else:
        tmp["__Q__"] = pd.NA
        tmp = tmp.sort_values(["Peptide"], kind="mergesort")
    best = tmp.groupby("Peptide", as_index=False).head(1).drop(columns="__Q__", errors="ignore")
    return best

def append_syn_to_fdr(fdr_path: Path, syn_path: Path, out_suffix="_withsyn.tsv") -> tuple[int, int]:
    fdr = load_safe_tsv(fdr_path, sep="\t")
    syn = load_safe_tsv(syn_path, sep="\t")

    # Required columns
    for col in ("Peptide", "Protein"):
        if col not in fdr.columns:
            raise ValueError(f"{fdr_path.name} missing required column: {col}")
        if col not in syn.columns:
            raise ValueError(f"{syn_path.name} missing required column: {col}")

    # Build sets of peptide–protein pairs
    fdr_pairs = fdr[["Peptide", "Protein"]].drop_duplicates()
    syn_pairs = syn[["Scan","Peptide", "Protein"]].drop_duplicates() #added Scan

    # Only consider peptides already present in FDR (so we can template-copy)
    syn_pairs = syn_pairs[syn_pairs["Peptide"].isin(fdr["Peptide"])]

    # New pairs = in SYN but not in FDR
    new_pairs = syn_pairs.merge(
        fdr_pairs.assign(_in_fdr=True),
        on=["Peptide", "Protein"],
        how="left"
    )
    new_pairs = new_pairs[new_pairs["_in_fdr"].isna()][["Peptide", "Protein"]]

    if new_pairs.empty:
        out_path = fdr_path.with_name(fdr_path.stem.replace("_fdrstats", "") + out_suffix)
        fdr.to_csv(out_path, sep="\t", index=False)
        return (0, fdr_pairs.shape[0])

    # Template row per peptide from FDR
    template = pick_template_rows(fdr)

    # Merge new pairs to their template rows
    # Left side has 'Protein' (the SYN protein we want to insert)
    # Right side is the template row (its 'Protein' will be suffixed to 'Protein_tmpl')
    new_rows = new_pairs.merge(template, on="Peptide", how="left", suffixes=("", "_tmpl"))

    # Replace template protein with the SYN protein (already in 'Protein')
    # Drop the right-side 'Protein_tmpl' if present
    if "Protein_tmpl" in new_rows.columns:
        new_rows = new_rows.drop(columns=["Protein_tmpl"])

    # Ensure/compute IsDecoy for the appended rows only (existing FDR rows are untouched)
    if "IsDecoy" not in new_rows.columns:
        # If the template lacked IsDecoy, create it
        new_rows["IsDecoy"] = new_rows["Protein"].astype(str).map(is_decoy)
    else:
        # Overwrite for appended rows per your spec
        new_rows["IsDecoy"] = new_rows["Protein"].astype(str).map(is_decoy)

    # Align appended-row columns to FDR column order
    new_rows = new_rows[[c for c in fdr.columns if c in new_rows.columns]]

    # Append; DO NOT drop duplicates globally (we must not remove any existing FDR rows)
    combined = pd.concat([fdr, new_rows], ignore_index=True)

    # Write output alongside the original FDR file
    out_name = fdr_path.stem.replace("_fdrstats", "") + out_suffix
    out_path = fdr_path.with_name(out_name)
    combined.to_csv(out_path, sep="\t", index=False)

    # Report: how many appended and now how many unique peptide–protein pairs
    return (len(new_rows), combined[["Peptide", "Protein"]].drop_duplicates().shape[0])

def find_pairs(sic_dir: Path, fdr_dir="fdr_estd", syn_dir="results/PHRPOut"):
    """
    Find matching FDR/SYN file pairs.

    Parameters:
        sic_dir: Path to SIC directory (user-provided)
        fdr_dir: Directory containing FDRStats (inside sic_dir)
        syn_dir: Name of SYN directory relative to SICdir's parent
    """
    sic_dir = sic_dir.resolve()
    fdr_root = sic_dir / fdr_dir
    syn_root = sic_dir.parent / syn_dir  # SYN next to SICdir

    if not fdr_root.is_dir():
        raise FileNotFoundError(f"Missing FDR directory: {fdr_root}")
    if not syn_root.is_dir():
        raise FileNotFoundError(f"Missing SYN directory: {syn_root}")

    pairs = []
    for fdr_file in fdr_root.glob("*_fdrstats.tsv"):
        base = fdr_file.name[:-len("_fdrstats.tsv")]
        syn_file = syn_root / f"{base}_syn.txt"
        if syn_file.exists():
            pairs.append((fdr_file, syn_file))
        else:
            print(f"[WARN] No matching SYN for {fdr_file.name} (looked for {syn_file.name})", file=sys.stderr)
    return pairs

def main():
    ap = argparse.ArgumentParser(description="Append peptide→protein combos from PHRP SYN into FDRStats.")
    ap.add_argument("-s", "--sic-dir", type=Path, required=True,
                    help="SIC directory (FDRStats inside; SYN directory expected next to it)")
    ap.add_argument("--fdr-dir", default="fdr_estd", help="FDRStats folder inside SIC directory")
    ap.add_argument("--syn-dir", default="results/PHRPOut",
                    help="SYN folder relative to parent of SIC directory")
    ap.add_argument("--suffix", default="_withsyn.tsv", help="Suffix for output files")
    args = ap.parse_args()

    pairs = find_pairs(args.sic_dir, fdr_dir=args.fdr_dir, syn_dir=args.syn_dir)
    if not pairs:
        print("No file pairs found. Nothing to do.")
        return

    total_added = 0
    for fdr_path, syn_path in pairs:
        added, unique_pairs = append_syn_to_fdr(fdr_path, syn_path, out_suffix=args.suffix)
        print(f"[OK] {fdr_path.name}  +{added} rows  -> unique peptide–protein pairs now: {unique_pairs}")
        total_added += added

    print(f"\nDone. Total new rows appended across all files: {total_added}")


if __name__ == "__main__":
    main()