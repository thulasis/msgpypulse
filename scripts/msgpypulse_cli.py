#!/usr/bin/env python3
import os
import argparse
import sys
from typing import Dict, Optional
import time
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import msgpypulse
from pathlib import Path
from glob import glob

class MSGFPDOWNPipeline:
    def __init__(self):
        self.log("Initializing MSGFPDOWN Pipeline")

    def log(self, message):
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{timestamp}] MSGFPDOWN: {message}")

    def prepare_output_dir(self, path):
        """Remove existing directory if it exists, then recreate it"""
        if os.path.exists(path) and os.path.isdir(path):
            self.log(f"Removing existing directory: {path}")
            shutil.rmtree(path)
        os.makedirs(path, exist_ok=True)

    def prepare_sics(self, input_dir, output_dir="SICs"):
        self.log("Preparing SICs")
        self.prepare_output_dir(output_dir)
        try:
            msgpypulse.prepare_SICs(input_dir=input_dir, output_dir=output_dir) #fixed
            return True
        except Exception as e:
            self.log(f"prepare_SICs failed: {e}")
            return False

    def fdr_estimator(self, input_dir="SICs", output_dir="fdr_estd"):
        self.log("FDR estimation")
        self.prepare_output_dir(output_dir)
        try:
            msgpypulse.process_directory(input_dir, output_dir) #fixed
            return True
        except Exception as e:
            self.log(f"fdr_estimator failed: {e}")
            return False

    def merge_syn_peptides(self, sic_dir, fdr_dir="fdr_estd", syn_dir="results/PHRPOut", suffix="_withsyn.tsv"):
        """
        Merge shared peptides from PHRPOut located next to SICdir.
        """
        self.log("Merging shared peptides from PHRPOut")
        try:
            sic_dir = Path(sic_dir).resolve()
            parent_dir = sic_dir.parent  # Parent directory of SICdir
            fdr_dir = parent_dir / fdr_dir
            syn_dir = parent_dir / syn_dir  # PHRPOut location

            pairs = msgpypulse.find_pairs(sic_dir, fdr_dir=fdr_dir, syn_dir=syn_dir)

            if not pairs:
                self.log(f"No matching FDR/SYN file pairs found in {syn_dir}")
                return False

            total_added = 0
            for fdr_path, syn_path in pairs:
                added, unique_pairs = msgpypulse.append_syn_to_fdr(
                    fdr_path, syn_path, out_suffix=suffix
                )
                self.log(f"{fdr_path.name}: +{added} rows (unique pairs={unique_pairs})")
                total_added += added

            self.log(f"Done. Total new rows appended: {total_added}")
            return True
        except Exception as e:
            self.log(f"merge_syn_peptides failed: {e}")
            return False

    def true_pep_filter(self, input_dir="fdr_estd", output_dir="fdr_filt", config: Dict = None):
        """Filter peptides from FDR results using configurable metrics and thresholds."""
        self.log("True peptide filtering")
        self.prepare_output_dir(output_dir)

        # Default values
        metric1, threshold1, metric2, threshold2 = None, None, None, None

        if config:
            if config.get("score_field"):
                metric1 = config["score_field"]
            if config.get("threshold"):
                threshold1 = float(config["threshold"])
            if config.get("score_field2"):
                metric2 = config["score_field2"]
            if config.get("threshold2"):
                threshold2 = float(config["threshold2"])

        try:                                                #fixed
            msgpypulse.process_fdrdir(
                input_directory=input_dir,
                output_directory=output_dir,
                metric1=metric1,
                threshold1=threshold1,
                metric2=metric2,
                threshold2=threshold2
            )
            return True
        except Exception as e:
            self.log(f"true_pep_filter failed: {e}")
            return False

    def peptide_ctab_generator(self, input_dir="fdr_filt", output_file="peptide_crosstab.tsv"):
        self.log("Generating peptide crosstab")
        try:
            in_dir = Path(input_dir)
            tsvs = sorted([p for p in in_dir.iterdir() if p.is_file() and p.name.endswith("_filtered.tsv")])
            if not tsvs:
                raise FileNotFoundError(f"No *_filtered.tsv files found in {in_dir}")

            # call the existing function from your module
            peptide_crosstab = msgpypulse.merge_peptide_ctab(tsvs) #fixed

            out_path = Path(output_file)
            out_path.parent.mkdir(parents=True, exist_ok=True)
            peptide_crosstab.to_csv(out_path, sep="\t", index=False)

            self.log(f"Peptide cross-tab written: {out_path} "
                    f"(peptides: {len(peptide_crosstab)}, samples: {len(tsvs)})")
            return True
        except Exception as e:
            self.log(f"peptide_ctab_generator failed: {e}")
            return False

    def fasta_cluster_generator(self, input_path=None, cdhit_path=None):
        self.log("Generating FASTA clusters")
        try:
            # Resolve cd-hit path using module function
            cdhit_exec = msgpypulse.resolve_cdhit_path(cdhit_path) #fixed

            # Determine input FASTA files
            if input_path:
                input_path = os.path.abspath(input_path)
                if os.path.isfile(input_path):
                    input_fastas = [input_path]
                elif os.path.isdir(input_path):
                    input_fastas = glob(os.path.join(input_path, "*.fasta")) + glob(os.path.join(input_path, "*.faa"))
                    if not input_fastas:
                        raise FileNotFoundError(f"No FASTA/FAA files found in {input_path}")
                else:
                    raise FileNotFoundError(f"Input path not found: {input_path}")
            else:
                input_fastas = msgpypulse.find_input_fastas() #fixed

            # Process each FASTA
            for fasta in input_fastas:
                msgpypulse.process_single_fasta(fasta, cdhit_exec) #fixed

            return True
        except Exception as e:
            self.log(f"fasta_cluster_generator failed: {e}")
            return False

    def peptide_ctab_annotator(self, input_file, output_file, database_path):
        self.log("Annotating peptide crosstab")
        try:
            in_path = Path(input_file)
            fasta_path = Path(database_path)
            out_path = Path(output_file)

            # Load peptide cross-tab
            df = pd.read_csv(in_path, sep="\t", dtype=str)
            required_cols = ["Peptide", "PeptideFlanked", "Protein"]
            for c in required_cols:
                if c not in df.columns:
                    raise ValueError(f"Missing required column: {c}")
                df[c] = df[c].astype(str).str.strip()

            sample_cols = [c for c in df.columns if c not in set(required_cols + 
                            ["Unique", "Shared", "Unique to Cluster", "Gene", "Function", "Cluster"])]

            # Build protein info map
            prot2info = msgpypulse.build_protein_info_map(fasta_path) #fixed

            # Map Cluster, Function, Gene
            df["Cluster"] = df["Protein"].map(lambda x: prot2info.get(str(x), {}).get("Cluster"))
            info_series = df["Protein"].apply(lambda pid: msgpypulse.get_protein_function(str(pid), prot2info))
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

            self.log(f"Wrote annotations → {out_path} (rows: {len(df)}, samples: {len(sample_cols)})")
            return True

        except Exception as e:
            self.log(f"peptide_ctab_annotator failed: {e}")
            return False
    
    def peptide_protein_mapbuilder(
            self,
            input_file,
            database_dir,
            num_pep: int = 2,
            mode="all_matches",
            output_dir=None,
            sample_cols=None,
            sic_dir = None
        ):
        self.log("Building peptide-protein maps")
        try:
            in_path = Path(input_file)
            db_dir = Path(database_dir)
            outdir = Path(output_dir) if output_dir else Path(sic_dir).parent / "map_files"
            self.prepare_output_dir(str(outdir))

            if not in_path.is_file():
                raise FileNotFoundError(f"Input table not found: {in_path}")

            # Load input table and detect samples
            df = pd.read_csv(in_path, sep="\t", low_memory=False, dtype=str)
            msgpypulse.ensure_required_columns(df)
            if sample_cols is None:
                sample_cols = msgpypulse.detect_sample_columns(df)
        
            if not sample_cols:
                raise ValueError(f"No sample columns detected. Available columns: {list(df.columns)}")
            self.log(f"Detected sample columns: {sample_cols}")

            if sic_dir is None:
                raise ValueError("Failed to locate inputfile.tsv")
            inputfile_path = Path(sic_dir).parent / "inputfile.tsv"

            # Load FASTA(s) for samples
            fasta_dict, valid_sample_cols = msgpypulse.load_fasta_for_samples(db_dir, sample_cols, inputfile_path=inputfile_path) #fixed
            self.log(f"Successfully loaded FASTA data for {len(valid_sample_cols)} samples")

            # Per-sample coverage
            per_sample_tables = []
            for sample_col in valid_sample_cols:
                sub = msgpypulse.prepare_subtable_for_sample(df, sample_col, mode) #fixed
                cov_tbl = msgpypulse.compute_coverage_for_sample(sub, fasta_dict[sample_col], num_pep, sample_col) #fixed
                if not cov_tbl.empty:
                    out_map = outdir / f"{sample_col}_mapfile.tsv"
                    cov_tbl.rename(columns={sample_col: "Coverage"}).to_csv(out_map, sep="\t", index=False)
                    per_sample_tables.append(cov_tbl[["Protein", sample_col]].copy())
                    self.log(f"Wrote mapfile for {sample_col}: {len(cov_tbl)} proteins")
                else:
                    self.log(f"Warning: No proteins found for sample {sample_col}")

            if not per_sample_tables:
               self.log("No proteins passed filters for any samples; nothing to write.")
               return True

            # Merge per-sample coverage
            merged = None
            for tbl in per_sample_tables:
                merged = tbl.copy() if merged is None else merged.merge(tbl, on="Protein", how="outer")

            if merged is None or merged.empty:
                self.log("No proteins passed filters; nothing to write.")
                return True

            # Tidy & write merged matrix
            sample_cols_sorted = sorted([c for c in merged.columns if c != "Protein"])
            merged = merged[["Protein"] + sample_cols_sorted]
            for c in sample_cols_sorted:
                merged[c] = pd.to_numeric(merged[c], errors="coerce").fillna(0)

            out_matrix = outdir / "grouped_coverage.tsv"
            merged.to_csv(out_matrix, sep="\t", index=False)

            # Generate heatmap
            try:
                pdf_path = msgpypulse.generate_coverage_heatmap(merged, outdir)
                if pdf_path:
                    self.log(f"Wrote heatmap: {pdf_path}")
                else:
                    self.log("Heatmap skipped (not enough data)")
            except Exception as e:
                self.log(f"Heatmap generation skipped due to error: {e}")

            self.log(f"Wrote per-sample mapfiles: {outdir}")
            self.log(f"Wrote merged coverage matrix: {out_matrix}  (proteins: {len(merged)})")
            return True

        except Exception as e:
            self.log(f"peptide_protein_mapbuilder failed: {e}")
            return False
    
    def protein_rollup(
        self,
        input_tsv,
        output_protein_tsv,
        rollup="sum",
        mode="all_matches",
        outlier_alpha=None,
        coverage_tsv=None
    ):
        """
        Wrapper for protein rollup using the core functions.
        """
        try:
            msgpypulse.run_rollup(
                input_tsv=input_tsv,
                output_protein_tsv=output_protein_tsv,
                rollup=rollup,
                mode=mode,
                outlier_alpha=outlier_alpha,
                coverage_tsv=coverage_tsv,
            )
            return True
        except Exception as e:
            self.log(f"protein_rollup failed: {e}")
            return False

    def cluster_ctab_generator(self, input_file, output_file, rollup="sum", mode="all_matches", outlier_alpha=None, coverage_tsv=None):
        self.log("Generating cluster crosstab")

        try:
            # Call the core function directly with keyword arguments
            msgpypulse.cluster_ctab_generator(
                input_tsv=input_file,
                output_cluster_tsv=output_file,
                rollup=rollup,
                mode=mode,
                outlier_alpha=outlier_alpha,
                coverage_tsv=coverage_tsv
            )
            return True
        except Exception as e:
            self.log(f"cluster_ctab_generator failed: {e}")
            return False

    def run_pipeline(self, config):
        start_time = time.time()
        try:
            # Resolve SICdir and parent
            sic_dir = Path(config["sic_input_dir"]).resolve()
            parent_dir = sic_dir.parent

            # Step 1: Prepare SICs
            sic_output_dir = parent_dir / "SICs"
            if not self.prepare_sics(sic_dir, output_dir=sic_output_dir):
                raise RuntimeError("SICs preparation failed")

            # Step 2: FDR estimation
            fdr_dir = parent_dir / "fdr_estd"
            if not self.fdr_estimator(input_dir=sic_output_dir, output_dir=fdr_dir):
                raise RuntimeError("FDR estimation failed")

            # Step 3: Merge shared peptides
            syn_dir = parent_dir / "results/PHRPOut"
            if not self.merge_syn_peptides(sic_dir=sic_output_dir, fdr_dir=fdr_dir, syn_dir=syn_dir):
                raise RuntimeError("Merging shared peptides failed")

            # Step 4: True peptide filtering
            fdr_filt_dir = parent_dir / "fdr_filt"
            if not self.true_pep_filter(input_dir=fdr_dir, output_dir=fdr_filt_dir, config=config):
                raise RuntimeError("Filtering failed")

            # Step 5: Generate peptide crosstab
            peptide_crosstab_file = parent_dir / "peptide_crosstab.tsv"
            if not self.peptide_ctab_generator(input_dir=fdr_filt_dir, output_file=peptide_crosstab_file):
                raise RuntimeError("Peptide crosstab generation failed")

            # Step 6: Generate FASTA clusters
            if not self.fasta_cluster_generator(input_path=config["database_path"]):
                raise RuntimeError("Clustering failed")

            # Step 7: Annotate peptide crosstab
            annotated_file = parent_dir / "peptide_crosstab_annotated.tsv"
            if not self.peptide_ctab_annotator(peptide_crosstab_file, annotated_file, config["database_path"]):
                raise RuntimeError("Peptide annotation failed")

            # Detect sample columns for mapping
            df_annot = pd.read_csv(annotated_file, sep="\t", dtype=str)
            required_cols = ["Peptide", "PeptideFlanked", "Protein", "Cluster",
                             "Unique to Cluster", "Unique", "Shared", "Gene", "Function"]
            sample_cols = [c for c in df_annot.columns if c not in required_cols]
            self.log(f"Detected sample columns for mapping: {sample_cols}")

            # Step 8: Build peptide-protein maps
            if not self.peptide_protein_mapbuilder(
                input_file=annotated_file,
                database_dir=config["database_path"],
                num_pep=config.get("num_pep", 2),
                mode=config.get("mode", "all_matches"),
                sample_cols=sample_cols,
                sic_dir=sic_dir
            ):
                raise RuntimeError("Peptide Protein Mapping failed")

            # Step 9: Protein rollup
            if not self.protein_rollup(
                annotated_file,
                config["final_output"],
                rollup=config.get("rollup", "sum"),
                mode=config.get("mode", "all_matches"),
                outlier_alpha=config.get("outlier_alpha"),
                coverage_tsv=config.get("coverage_tsv"),
            ):
                raise RuntimeError("Rollup failed")

            # Step 10: Optional cluster crosstab generation
            cluster_out = None
            if config.get("run_cluster_step", False):
                cluster_out = config.get("cluster_output")
                if not cluster_out:
                    raise RuntimeError("Cluster rollup requested but no cluster_output specified")
                self.log("Running optional cluster crosstab generation...")

                try:
                    msgpypulse.run_rollup_cluster(
                        input_tsv=annotated_file,
                        output_cluster_tsv=cluster_out,
                        rollup=config.get("rollup", "sum"),
                        mode=config.get("mode", "all_matches"),
                        outlier_alpha=config.get("outlier_alpha"),
                        coverage_tsv=config.get("coverage_tsv")
                    )
                    self.log(f"Cluster crosstab generated successfully at: {cluster_out}")

                except Exception as e:
                    raise RuntimeError(f"Cluster rollup failed: {e}")
            
            else:
                self.log("Cluster rollup step skipped (run_cluster_step=False).")

            total_time = time.time() - start_time
            self.log(f"Pipeline completed successfully! Total time: {total_time:.2f} seconds")
            return True

        except Exception as e:
            total_time = time.time() - start_time
            self.log(f"Pipeline failed after {total_time:.2f} seconds: {str(e)}")
            return False

def main():
    print("Running MS-GF+ Donwstream CLI pipeline")
    parser = argparse.ArgumentParser(description="MS-GF+ Downstream Pipeline Controller")
    parser.add_argument("-i", "--input", required=True, help="Input SIC directory")
    parser.add_argument("-o", "--output", required=True, help="Final output file path")
    parser.add_argument("--database", required=True, help="Database directory path")
    parser.add_argument("--score_field", default="MSMSScore", help="Score field for Step 4 filtering (default: MSMSScore)")
    parser.add_argument("--threshold", type=int, default=10, help="Threshold for Step 4 filtering (default: 10)")
    parser.add_argument("--score_field2", help="Optional second score field for Step 4 filtering (-f2)")
    parser.add_argument("--threshold2", type=float, help="Optional second threshold for Step 4 filtering (-t2)")
    parser.add_argument("--num_pep", type=int, default=2, help="Minimum number of peptides for Step 8")
    parser.add_argument("--mode", default="all_matches", choices=["unique_only", "requires_unique", "all_matches"], help="Mode for Steps 8 and 9")
    parser.add_argument("--rollup", default="sum", choices=["sum", "rrollup", "zrollup", "qrollup", "weighted"], help="Rollup method for Step 9")
    parser.add_argument("--outlier_alpha", type=float, help="Grubbs alpha for outlier detection (RRrollup only)")
    parser.add_argument("--coverage_tsv", help="Path to grouped_coverage.tsv (optional)")
    parser.add_argument("--cluster_rollup", action="store_true", help="Run optional Step 10 (cluster crosstab)")
    parser.add_argument("-co", "--cluster_output", help="Output file for Step 10 cluster crosstab")

    args = parser.parse_args()

    if args.cluster_rollup and not args.cluster_output:
        parser.error("--cluster_output (-co) is required when --cluster_rollup is specified")

    config = {
        "sic_input_dir": args.input,
        "database_path": args.database,
        "score_field": args.score_field,
        "threshold": args.threshold,
        "score_field2": args.score_field2,
        "threshold2": args.threshold2,
        "num_pep": args.num_pep,
        "mode": args.mode,
        "rollup": args.rollup,
        "outlier_alpha": args.outlier_alpha,
        "coverage_tsv": args.coverage_tsv,
        "run_cluster_step": args.cluster_rollup,
        "cluster_output": args.cluster_output,
        "final_output": args.output,
    }

    pipeline = MSGFPDOWNPipeline()
    success = pipeline.run_pipeline(config)
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()