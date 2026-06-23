#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
import re

# ---- Helpers ----
def get_gene_from_info(info: str):
    if not isinstance(info, str) or not info:
        return "Unknown"
    m = re.search(r'(?:^|\s|;|\|)(?:GN|Gene)=([^\s;|]+)', info, flags=re.IGNORECASE)
    return m.group(1) if m else "Unknown"

def _remove_cluster_tokens(text: str) -> str:
    if not isinstance(text, str):
        return ""
    cleaned = re.sub(r'(?:^|\s|;|\|)Cluster=[^\s;|]+', ' ', text, flags=re.IGNORECASE)
    return " ".join(cleaned.split())

def get_protein_function(protein_id: str, proteins: dict) -> dict:
    protein_info = proteins.get(protein_id)
    if not protein_info:
        return {'function': 'Unknown', 'gene': 'Unknown'}
    function_string = protein_info.get('Function', '') or ''
    if "OS" in function_string:
        parts = function_string.split('OS', 1)
        function = _remove_cluster_tokens(parts[0].strip())
        info = "OS" + parts[1].strip() if len(parts) > 1 else ''
    else:
        function = _remove_cluster_tokens(function_string.strip())
        info = ''
    return {'function': function or 'Unknown', 'gene': get_gene_from_info(info) or 'Unknown'}

def build_protein_info_map(fasta_paths, cluster_tag: str = "Cluster="):
    """
    Parse one or more FASTA files and return a combined map:
        { protein_id -> {"Cluster": str or None, "Function": header tail} }
    """
    prot2info = {}
    cluster_re = re.compile(r'(?:^|\s|;|\|)' + re.escape(cluster_tag) + r'([^\s;|]+)', flags=re.IGNORECASE)

    # Ensure fasta_paths is a list
    if isinstance(fasta_paths, (str, Path)):
        fasta_paths = [fasta_paths]

    for fasta_path in fasta_paths:
        fasta_path = Path(fasta_path)
        if not fasta_path.exists():
            print(f"Warning: FASTA not found, skipping: {fasta_path}")
            continue

        # If directory, glob for fasta/faa files
        if fasta_path.is_dir():
            files = sorted(fasta_path.glob("*.fasta")) + sorted(fasta_path.glob("*.faa"))
            if not files:
                print(f"Warning: No FASTA files found in directory: {fasta_path}")
                continue
        else:
            files = [fasta_path]

        for f in files:
            with open(f, "r", encoding="utf-8", errors="ignore") as fh:
                for line in fh:
                    if not line.startswith(">"):
                        continue
                    header = line[1:].rstrip("\n")
                    protein_id = header.split()[0]
                    rest = header[len(protein_id):].lstrip()
                    m = cluster_re.search(header)
                    clust = m.group(1) if m else None
                    prot2info[protein_id] = {"Cluster": clust, "Function": rest}

    return prot2info

# ---- Main ----
def main():
    parser = argparse.ArgumentParser(
        description="Annotate peptide–protein crosstab with Cluster/Gene/Function and uniqueness metrics."
    )
    parser.add_argument("-i", "--crosstab-tsv", required=True, help="Input peptide–protein crosstab TSV.")
    parser.add_argument("-f", "--fasta", required=True,
                        help="Path to a FASTA file or a directory containing FASTA (*.fasta, *.faa) files.")
    parser.add_argument("-o", "--output-tsv", required=True, help="Output annotated TSV path.")
    args = parser.parse_args()

    in_path = Path(args.crosstab_tsv)
    fasta_path = Path(args.fasta)
    out_path = Path(args.output_tsv)

    # Load peptide cross-tab
    df = pd.read_csv(in_path, sep="\t", dtype=str)
    required_cols = ["Peptide", "PeptideFlanked", "Protein"]
    for c in required_cols:
        if c not in df.columns:
            raise ValueError(f"Missing required column: {c}")
        df[c] = df[c].astype(str).str.strip()

    sample_cols = [c for c in df.columns if c not in set(required_cols + ["Unique", "Shared", "Unique to Cluster", "Gene", "Function", "Cluster"])]

    # Build protein info map (handles single FASTA or directory)
    prot2info = build_protein_info_map(fasta_path)

    # Map Cluster, Function, Gene
    df["Cluster"] = df["Protein"].map(lambda x: prot2info.get(str(x), {}).get("Cluster"))
    info_series = df["Protein"].apply(lambda pid: get_protein_function(str(pid), prot2info))
    df["Function"] = info_series.apply(lambda d: d["function"])
    df["Gene"] = info_series.apply(lambda d: d["gene"])

    # Compute uniqueness
    df["_nprot_"] = df.groupby("PeptideFlanked")["Protein"].transform("nunique")
    df["Unique"] = df["_nprot_"] == 1
    df.drop(columns=["_nprot_"], inplace=True)

    df["_nclus_"] = df.groupby("PeptideFlanked")["Cluster"].transform("nunique")
    df["Unique to Cluster"] = df["_nclus_"] == 1
    df.drop(columns=["_nclus_"], inplace=True)

    df["Shared"] = ~df["Unique"]

    # Coerce sample columns
    for c in sample_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)

    # Reorder columns
    front = ["Peptide", "PeptideFlanked", "Protein"]
    tail = ["Cluster", "Unique to Cluster", "Unique", "Shared", "Gene", "Function"]
    middle = [c for c in df.columns if c not in set(front + tail)]
    df = df[front + middle + tail]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)
    print(f"Wrote annotations → {out_path} (rows: {len(df)}, samples: {len(sample_cols)})")


if __name__ == "__main__":
    main()