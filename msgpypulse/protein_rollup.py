#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
from scipy import stats
from sklearn.preprocessing import quantile_transform
import warnings
import re
warnings.filterwarnings("ignore")

# ---------- Helpers ----------
def _coerce_bool_series(s):
    """Convert series to boolean, handling various input formats."""
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False)
    
    # Handle string and numeric representations
    return s.fillna(False).apply(lambda x: str(x).strip().lower() in ("true", "1", "t", "yes", "y"))

def _is_boolish_col(col):
    """Check if column contains boolean-like data."""
    if pd.api.types.is_bool_dtype(col):
        return True
    
    # Check numeric 0/1 columns
    if pd.api.types.is_numeric_dtype(col):
        unique_vals = set(pd.unique(col.dropna()))
        return unique_vals <= {0, 1}
    
    # Check string boolean representations
    uniq_vals = set(str(v).strip().lower() for v in pd.unique(col.dropna()))
    return uniq_vals <= {"true", "false", "t", "f", "yes", "no", "y", "n", "0", "1"}

def _infer_sample_cols(df, exclude=()):
    """Infer which columns contain sample intensity data."""
    exclude = set(exclude)
    
    # First try: numeric columns that aren't boolean-like
    numeric = [c for c in df.columns 
               if c not in exclude 
               and pd.api.types.is_numeric_dtype(df[c]) 
               and not _is_boolish_col(df[c])]
    
    if numeric:
        return numeric
    
    # Fallback: any non-boolean columns
    return [c for c in df.columns if c not in exclude and not _is_boolish_col(df[c])]

def _pick_peptide_key(df, prefer_flanked=True):
    """Select the best peptide identifier column."""
    if prefer_flanked and "PeptideFlanked" in df.columns:
        return "PeptideFlanked"
    if "Peptide" in df.columns:
        return "Peptide"
    
    # Look for columns with "peptide" in the name (case insensitive)
    peptide_cols = [col for col in df.columns if 'peptide' in col.lower()]
    if peptide_cols:
        return peptide_cols[0]
    
    raise ValueError("Input must contain a peptide identifier column (Peptide, PeptideFlanked, or similar)")

def _filter_to_coverage_ids(rolled_df, coverage_tsv):
    """Filter results to proteins present in grouped coverage file."""
    if not coverage_tsv.is_file():
        raise FileNotFoundError(f"Grouped coverage file not found: {coverage_tsv}")
    
    try:
        cov = pd.read_csv(coverage_tsv, sep="\t")
    except Exception as e:
        raise SystemExit(f"[ERROR] Failed reading coverage file {coverage_tsv}: {e}")
    
    cov_ids = set(cov["Protein"].astype(str))
    filtered_df = rolled_df[rolled_df["Protein"].astype(str).isin(cov_ids)].copy()
    return filtered_df

def _natural_sort_key(s):
    """Natural sorting key for sample column names."""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', str(s))]

def _validate_input_data(df, sample_cols):
    """Validate input data quality."""
    if df.empty:
        raise ValueError("Input dataframe is empty")
    
    if not sample_cols:
        raise ValueError("No sample columns found in the data")
    
    # Check for required columns
    if "Protein" not in df.columns:
        raise ValueError("Required column 'Protein' is missing")

# ---------- RRollup internals ----------
def _pick_reference(peptab, sample_cols):
    """Select reference peptide with least missing values and highest total intensity."""
    if peptab.empty:
        raise ValueError("Empty peptide table for reference selection")
    
    miss = peptab[sample_cols].isna().sum(axis=1)
    totals = peptab[sample_cols].sum(axis=1, skipna=True)
    
    order = pd.DataFrame({
        'idx': peptab.index, 
        'miss': miss, 
        'tot': totals
    })
    
    best = order.sort_values(['miss', 'tot'], ascending=[True, False]).iloc[0]
    return int(best['idx'])

def _median_ratio(ref_vec, pep_vec):
    """Calculate median ratio between reference and peptide vectors."""
    mask = (~np.isnan(ref_vec)) & (~np.isnan(pep_vec)) & (pep_vec != 0) & (ref_vec != 0)
    
    if not mask.any():
        return np.nan
    
    ratios = ref_vec[mask] / pep_vec[mask]
    ratios = ratios[np.isfinite(ratios)]
    
    return float(np.median(ratios)) if len(ratios) > 0 else np.nan

def _grubbs_filter(values, alpha=0.05):
    """Apply Grubbs test to filter outliers."""
    x = values.astype(float)
    mask = ~np.isnan(x)
    n = mask.sum()
    
    if n < 3:
        return mask
    
    vals = x[mask]
    mean_val = vals.mean()
    sd = vals.std(ddof=1)
    
    if sd == 0:
        return mask
    
    z = np.abs(vals - mean_val) / sd
    i = z.argmax()
    G = z[i]
    
    t_crit = stats.t.ppf(1 - alpha / (2 * n), n - 2)
    G_crit = ((n - 1) / np.sqrt(n)) * np.sqrt(t_crit**2 / (n - 2 + t_crit**2))
    
    if G > G_crit:
        outlier_idx = np.flatnonzero(mask)[i]
        mask[outlier_idx] = False
    
    return mask

def _first_non_null(x):
    """Get first non-null value from series."""
    for v in x:
        if pd.notna(v) and str(v).strip() != "":
            return v
    return ""

# ---------- Core Rollup ----------
def rollup_from_annotated(
    pep_annot_tsv,
    rollup="sum",
    mode="all_matches",
    outlier_alpha=None
):
    """
    Perform protein rollup from annotated peptide data.
    
    Parameters:
    -----------
    pep_annot_tsv
        Path to annotated peptide TSV file
    rollup
        Rollup method ('sum', 'rrollup', 'zrollup', 'qrollup', 'weighted')
    mode
        Peptide selection mode ('all_matches', 'unique_only', 'requires_unique')
    outlier_alpha
        Alpha level for Grubbs outlier test (RRollup only)
    """
    
    # Read and validate input
    try:
        df = pd.read_csv(pep_annot_tsv, sep="\t")
    except Exception as e:
        raise SystemExit(f"[ERROR] Failed reading {pep_annot_tsv}: {e}")
    
    pep_key = _pick_peptide_key(df, prefer_flanked=True)
    
    # Handle boolean columns properly
    if "Unique" not in df.columns:
        df["Unique"] = False
    df["Unique"] = _coerce_bool_series(df["Unique"])
    
    # Identify annotation columns
    gene_cols = [c for c in ("Gene", "gene") if c in df.columns]
    func_cols = [c for c in ("Function", "function") if c in df.columns]
    
    # Identify sample columns
    meta_cols = [pep_key, "Protein", "Unique", "Cluster"] + gene_cols + func_cols
    sample_cols = _infer_sample_cols(df, exclude=meta_cols)
    
    # Validate input
    _validate_input_data(df, sample_cols)
    
    # Prepare working dataframe
    required_cols = ["Protein", "Unique", pep_key]
    available_cols = [col for col in required_cols if col in df.columns]
    
    # Handle peptide columns
    have_plain = ("Peptide" in df.columns) and (pep_key != "Peptide")
    if not have_plain and pep_key != "Peptide":
        df["Peptide"] = df[pep_key]
    
    working_cols = available_cols + ["Peptide"] + sample_cols
    if "Cluster" in df.columns:
        working_cols.append("Cluster")
    working_cols.extend(gene_cols + func_cols)
    
    slim = df[working_cols].copy()
    slim = slim.rename(columns={pep_key: "Peptide_pref"})
    if have_plain:
        slim["Peptide_plain"] = df["Peptide"].astype(str)
    else:
        slim["Peptide_plain"] = slim["Peptide"].astype(str)
    
    # Apply mode filtering
    if mode == "unique_only":
        working = slim[slim["Unique"]]
    elif mode == "requires_unique":
        proteins_with_unique = slim.groupby("Protein")["Unique"].any()
        keep_proteins = proteins_with_unique[proteins_with_unique].index
        working = slim[slim["Protein"].isin(keep_proteins)]
    else:  # "all_matches"
        working = slim
    
    if working.empty:
        print("Warning: No peptides remaining after filtering")
        return pd.DataFrame(columns=["Protein"] + sample_cols + [
            "Peptide Number", "Peptides", "Flanked Peptides", 
            "Unique Peptide(s)", "Shared Peptide(s)"
        ])
    
    # Build annotation maps
    cluster_map = working.groupby("Protein")["Cluster"].first() if "Cluster" in working.columns else pd.Series(dtype=object)
    gene_map = {c: working.groupby("Protein")[c].apply(_first_non_null) for c in gene_cols}
    func_map = {c: working.groupby("Protein")[c].apply(_first_non_null) for c in func_cols}
    
    uniq_map = working.groupby("Protein")["Unique"].any()
    shared_map = working.groupby("Protein")["Unique"].apply(lambda s: (~s).any())
    
    # Build peptide maps for ALL methods (FIX for Bug 1)
    flanked_peptide_map = working.groupby("Protein")["Peptide_pref"].apply(
        lambda x: ";".join(sorted(set(x.dropna().astype(str))))
    )
    plain_peptide_map = working.groupby("Protein")["Peptide_plain"].apply(
        lambda x: ";".join(sorted({s for s in x.dropna().astype(str) if s}))
    )
    
    # Filter by grouped coverage and process each protein
    proteins_filtered_out = 0
    
    if rollup.lower() == "sum":
        # Use pandas groupby for sum aggregation
        agg_fn = {c: "sum" for c in sample_cols}
        
        grouped = working.groupby("Protein", as_index=False).agg({
            **agg_fn,
            "Peptide_pref": lambda x: ";".join(sorted(set(map(str, x)))),
            "Peptide_plain": lambda x: ";".join(sorted({s for s in map(str, x) if s})),
        })
    
    else:
        # Handle advanced rollup methods
        prot_rows = []
        
        for prot, peptab in working.groupby("Protein"):
            vals = peptab[sample_cols].values.astype(float)
            
            if rollup.lower() == "rrollup":
                if peptab.shape[0] == 1:
                    prot_vals = vals[0]
                else:
                    try:
                        ref_idx = _pick_reference(peptab, sample_cols)
                        ref = peptab.loc[ref_idx, sample_cols].values.astype(float)
                        
                        scaled = []
                        for _, row in peptab.iterrows():
                            v = row[sample_cols].values.astype(float)
                            if row.name == ref_idx:
                                v_scaled = v.copy()
                            else:
                                sf = _median_ratio(ref, v)
                                if np.isfinite(sf) and sf > 0:
                                    v_scaled = v * sf
                                else:
                                    v_scaled = np.full_like(v, np.nan)
                            scaled.append(v_scaled)
                        
                        scaled = np.vstack(scaled)
                        
                        # Apply outlier filtering if requested
                        if outlier_alpha is not None:
                            for j in range(scaled.shape[1]):
                                col_mask = _grubbs_filter(scaled[:, j], alpha=outlier_alpha)
                                scaled[~col_mask, j] = np.nan
                        
                        prot_vals = np.nanmedian(scaled, axis=0)
                        
                    except Exception as e:
                        warnings.warn(f"RRollup failed for protein {prot}: {e}. Using median.")
                        prot_vals = np.nanmedian(vals, axis=0)
            
            elif rollup.lower() == "zrollup":
                    # Strict validation with warnings instead of silent fixes
                if vals.shape[0] == 1:
                    warnings.warn(f"ZRollup not applicable for protein {prot}: only one peptide present")
                    prot_vals = np.full(vals.shape[1], np.nan)
                elif np.all(np.isnan(vals)):
                    warnings.warn(f"ZRollup not applicable for protein {prot}: all peptide values are NaN")
                    prot_vals = np.full(vals.shape[1], np.nan)
                else:
                    means = np.nanmean(vals, axis=0)
                    stds = np.nanstd(vals, axis=0, ddof=1)

                    if np.any(stds < 1e-10):
                        warnings.warn(f"ZRollup not applicable for protein {prot}: near-zero variance in peptides")
                        prot_vals = np.full(vals.shape[1], np.nan)
                    else:
                        vals_z = (vals - means) / stds
                        if not np.isfinite(vals_z).all():
                            warnings.warn(f"ZRollup not applicable for protein {prot}: non-finite z-scores encountered")
                            prot_vals = np.full(vals.shape[1], np.nan)
                        else:
                            prot_vals = np.nanmedian(vals_z, axis=0)
            
            elif rollup.lower() == "qrollup":
                if np.all(np.isnan(vals)):
                    prot_vals = np.full(vals.shape[1], np.nan)
                else:
                    vals_q = quantile_transform(vals, axis=0, copy=True, 
                                              output_distribution='uniform')
                    prot_vals = np.nanmedian(vals_q, axis=0)
            
            elif rollup.lower() == "weighted":
                peptide_counts = peptab.groupby("Peptide_pref").size()
                weights = peptab["Peptide_pref"].map(lambda x: 1.0 / peptide_counts[x])
                weighted_vals = vals * weights.values[:, np.newaxis]
                prot_vals = np.nansum(weighted_vals, axis=0)
            
            else:
                raise ValueError(f"Unsupported rollup method: {rollup}")
            
            prot_rows.append((prot, prot_vals, len(peptab)))
        
        if proteins_filtered_out > 0:
            print(f"Filtered out {proteins_filtered_out} proteins")
        
        if not prot_rows:
            print("Warning: No proteins remaining after filtering")
            return pd.DataFrame(columns=["Protein"] + sample_cols + [
                "Peptide Number", "Peptides", "Flanked Peptides",
                "Unique Peptide(s)", "Shared Peptide(s)"
            ])
        
        # Build result dataframe
        grouped = pd.DataFrame({"Protein": [p for p, _, _ in prot_rows]})
        
        # Add sample columns
        for j, col in enumerate(sample_cols):
            grouped[col] = [vals[j] for _, vals, _ in prot_rows]
        
        # Add peptide count
        grouped["Peptide Number"] = [count for _, _, count in prot_rows]
    
    # Add common annotations for ALL methods (FIX for Bug 1)
    if not cluster_map.empty:
        grouped["Cluster"] = grouped["Protein"].map(cluster_map)
    
    for c in gene_cols:
        grouped[c] = grouped["Protein"].map(gene_map[c])
    
    for c in func_cols:
        grouped[c] = grouped["Protein"].map(func_map[c])
    
    grouped["Unique Peptide(s)"] = grouped["Protein"].map(uniq_map).fillna(False).astype(bool)
    grouped["Shared Peptide(s)"] = grouped["Protein"].map(shared_map).fillna(False).astype(bool)
    
    # Add peptide information for ALL methods (FIX for Bug 1)
    grouped["Flanked Peptides"] = grouped["Protein"].map(flanked_peptide_map).fillna("")
    grouped["Peptides"] = grouped["Protein"].map(plain_peptide_map).fillna("")
    
    # Calculate peptide counts if not already present
    if "Peptide Number" not in grouped.columns:
        base_list_col = "Flanked Peptides" if grouped["Flanked Peptides"].astype(str).str.len().gt(0).any() else "Peptides"
        grouped["Peptide Number"] = grouped[base_list_col].apply(
            lambda s: 0 if pd.isna(s) or str(s) == "" else len(set(str(s).split(";")))
        )
    
    # Clean up internal columns
    for col in ["Peptide_pref", "Peptide_plain"]:
        if col in grouped.columns:
            grouped = grouped.drop(columns=[col])
    
    return grouped

def _finalize_column_order(df):
    """Organize columns in a logical order with naturally sorted sample names."""
    meta_cols = {"Protein", "Cluster", "Gene", "Function", "Peptide Number", 
                "Unique Peptide(s)", "Shared Peptide(s)", "Peptides", "Flanked Peptides"}
    
    # Identify and sort sample columns naturally
    sample_cols = [c for c in df.columns 
                   if c not in meta_cols 
                   and pd.api.types.is_numeric_dtype(df[c]) 
                   and not _is_boolish_col(df[c])]
    sample_cols = sorted(sample_cols, key=_natural_sort_key)
    
    # Define desired order
    desired_order = (["Protein"] + sample_cols + 
                    ["Cluster", "Gene", "Function", "Peptide Number",
                     "Unique Peptide(s)", "Shared Peptide(s)", "Peptides", "Flanked Peptides"])
    
    # Build final column list
    final_cols = [c for c in desired_order if c in df.columns]
    remaining_cols = [c for c in df.columns if c not in final_cols]
    
    return df[final_cols + remaining_cols]

def run_rollup(input_tsv, output_protein_tsv=None,
               rollup="sum", mode="all_matches",
               outlier_alpha=None, coverage_tsv=None):
    """
    Programmatic API that runs the same pipeline as main() but returns the rolled DataFrame.
    Raises exceptions to caller on failure.
    """
    input_path = Path(input_tsv).resolve()
    input_dir = input_path.parent

    rolled = rollup_from_annotated(
        pep_annot_tsv=input_path,
        rollup=rollup,
        mode=mode,
        outlier_alpha=outlier_alpha
    )

    # 2) Filter to coverage IDs
    cov_path = Path(coverage_tsv).resolve() if coverage_tsv else (input_dir / "map_files" / "grouped_coverage.tsv")
    rolled = _filter_to_coverage_ids(rolled, cov_path)

    # 3) Finalize column order
    rolled = _finalize_column_order(rolled)

    # 4) Optionally write output
    if output_protein_tsv is not None:
        out_path = Path(output_protein_tsv)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        rolled.to_csv(out_path, sep="\t", index=False, na_rep="NA")

    return rolled

# ---------- CLI ----------
def main():
    parser = argparse.ArgumentParser(
        description="Protein rollup from annotated peptide crosstab",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input_tsv", required=True, 
                       help="Annotated peptide crosstab TSV")
    parser.add_argument("-o", "--output_protein_tsv", required=True, 
                       help="Output protein-level TSV")
    parser.add_argument("--rollup", 
                       choices=["sum", "rrollup", "zrollup", "qrollup", "weighted"], 
                       default="sum",
                       help="Rollup method")
    parser.add_argument("--mode", 
                       choices=["all_matches", "unique_only", "requires_unique"], 
                       default="all_matches",
                       help="Peptide selection mode")
    parser.add_argument("--outlier_alpha", type=float, default=None, 
                       help="Grubbs alpha for outlier detection (RRollup only)")
    parser.add_argument("--coverage_tsv", default=None,
                       help="Path to grouped_coverage.tsv. If not provided, uses <input_dir>/map_files/grouped_coverage.tsv")
    
    args = parser.parse_args()
    
    input_path = Path(args.input_tsv).resolve()

    try:
        rolled = run_rollup(
            input_tsv=input_path,
            output_protein_tsv=args.output_protein_tsv,
            rollup=args.rollup,
            mode=args.mode,
            outlier_alpha=args.outlier_alpha,
            coverage_tsv=args.coverage_tsv
        )
        print(f"[INFO] Protein rollup completed successfully.")
        print(f"[INFO] Output: {args.output_protein_tsv}")
        print(f"[INFO] Proteins in output: {len(rolled)}")
        return 0
    
    except Exception as e:
        print(f"Error during rollup: {e}")
        return 1
    
if __name__ == "__main__":
    exit(main())