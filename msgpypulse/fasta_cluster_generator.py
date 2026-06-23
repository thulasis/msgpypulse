#!/usr/bin/env python3
import argparse
import subprocess
import os
import shutil
from Bio import SeqIO
from glob import glob
from pathlib import Path
import sys
import re
import tempfile

CDHIT_PATH = "bin/cd-hit"

def find_input_fastas(path=None):
    """Return list of FASTA/FAA files to process."""
    if path is None:
        path = "./database"
    path = os.path.abspath(path)
    if os.path.isfile(path):
        return [path]
    elif os.path.isdir(path):
        files = glob(os.path.join(path, "*.fasta")) + glob(os.path.join(path, "*.faa"))
        if not files:
            raise FileNotFoundError(f"No FASTA/FAA files found in {path}")
        return files
    else:
        raise FileNotFoundError(f"Input path not found: {path}")

def is_already_processed(fasta_path, max_checks=1000):
    count = 0
    try:
        with open(fasta_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    count += 1
                    if 'Cluster=' in line:
                        return True
                    if count >= max_checks:
                        break
    except Exception as e:
        print(f"[ERROR] Failed to check for existing cluster annotations: {e}")
    return False

def resolve_cdhit_path(provided_path=None):
    if provided_path:
        if os.path.isfile(provided_path) and os.access(provided_path, os.X_OK):
            return provided_path
        else:
            sys.exit(f"[ERROR] Provided cd-hit path '{provided_path}' is not executable.")
    found = shutil.which("cd-hit")
    if found:
        return found
    if CDHIT_PATH and os.path.isfile(CDHIT_PATH) and os.access(CDHIT_PATH, os.X_OK):
        return CDHIT_PATH
    sys.exit("[ERROR] cd-hit binary not found in PATH and no --cdhit-path provided.")

def run_cdhit_clustering(input_fasta, identity=0.95, cdhit_path=None):
    tmp_out = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta")
    tmp_out.close()
    cmd = [
        cdhit_path,
        "-i", input_fasta,
        "-o", tmp_out.name,
        "-c", str(identity),
        "-d", "0"
    ]
    print(f"[INFO] Running cd-hit: {' '.join(cmd)}")
    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, cmd)
    clstr_path = tmp_out.name + ".clstr"
    if not os.path.exists(clstr_path):
        raise FileNotFoundError(f"Expected cd-hit output not found: {clstr_path}")
    return tmp_out.name, clstr_path

def parse_clstr_file(clstr_path):
    cluster_map = {}
    cluster_counter = 0
    current_cluster = None
    with open(clstr_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">Cluster"):
                cluster_counter += 1
                current_cluster = f"Cluster_{cluster_counter:04d}"
                continue
            if current_cluster and '>' in line and '...' in line:
                seq_id = line.split('>',1)[1].split('...')[0].strip()
                cluster_map[seq_id] = current_cluster
    return cluster_map

def _norm_uniprot(description: str) -> str | None:
    if not description:
        return None
    first_token = description.strip().split(" ")[0]
    match = re.match(r"^[a-z]{2}\|([A-Z0-9_]+)\|", first_token, re.IGNORECASE)
    return match.group(1) if match else None

def _norm_pipe_accession(s):
    if '|' in s:
        parts = s.split('|')
        if len(parts) >= 2 and parts[1]:
            return parts[1]
    return None

def _norm_raw_id(s):
    return s

def _norm_description(s):
    return s

def _norm_strip_version(s):
    for sep in ('/', '.'):
        if sep in s:
            base, tail = s.rsplit(sep,1)
            if tail.isdigit():
                return base
    return s

def _norm_no_whitespace(s):
    return s.replace(" ", "")

NORMALIZERS = [
    ("uniprot", lambda rec: _norm_uniprot(rec.description)),
    ("pipe_accession", lambda rec: _norm_pipe_accession(rec.id)),
    ("raw_id", lambda rec: _norm_raw_id(rec.id)),
    ("description", lambda rec: _norm_description(rec.description)),
    ("strip_version", lambda rec: _norm_strip_version(rec.id)),
    ("no_whitespace", lambda rec: _norm_no_whitespace(rec.id)),
]

def annotate_fasta_with_clusters(input_fasta, cluster_map, output_fasta):
    def choose_normalizer_for_record(rec, cluster_keys):
        for name, func in NORMALIZERS:
            try:
                key = func(rec)
            except Exception:
                key = None
            if key and key in cluster_keys:
                return name, func
        return None, None

    cluster_keys = set(cluster_map.keys())
    with open(input_fasta, "r") as fh:
        first_record = next(SeqIO.parse(fh, "fasta"), None)
    if first_record is None:
        raise ValueError("Input FASTA appears empty.")
    name, func = choose_normalizer_for_record(first_record, cluster_keys)
    if name is None:
        name, func = "raw_id", lambda rec: rec.id
    print(f"[INFO] Initial normalizer chosen: {name}")

    out_handle = open(output_fasta, "w")
    temp_in_path = input_fasta
    current_normalizer = func
    iteration = 0

    while True:
        iteration += 1
        tmp_fd, temp_out_path = tempfile.mkstemp(suffix=".fasta", prefix="unmatched_")
        os.close(tmp_fd)
        unmatched_count = 0

        with open(temp_in_path, "r") as in_fh, open(temp_out_path, "w") as unmatched_fh:
            for rec in SeqIO.parse(in_fh, "fasta"):
                try:
                    key = current_normalizer(rec)
                except Exception:
                    key = None
                cluster_label = None
                if key and key in cluster_map:
                    cluster_label = cluster_map[key]
                else:
                    if rec.id in cluster_map:
                        cluster_label = cluster_map[rec.id]
                    elif rec.description in cluster_map:
                        cluster_label = cluster_map[rec.description]
                if cluster_label:
                    rec.description = f"{rec.description} Cluster={cluster_label.replace('Cluster_','')}"
                    SeqIO.write(rec, out_handle, "fasta")
                else:
                    SeqIO.write(rec, unmatched_fh, "fasta")
                    unmatched_count += 1

        print(f"[INFO] Iteration {iteration}: unmatched_count = {unmatched_count}")
        if unmatched_count == 0:
            break

        with open(temp_out_path, "r") as ufh:
            first_unmatched = next(SeqIO.parse(ufh, "fasta"), None)
        if first_unmatched is None:
            break
        new_name, new_func = choose_normalizer_for_record(first_unmatched, cluster_keys)
        if new_name is None or new_func == current_normalizer:
            with open(temp_out_path, "r") as ufh:
                for rec in SeqIO.parse(ufh, "fasta"):
                    rec.description = f"{rec.description} Cluster=Unassigned"
                    SeqIO.write(rec, out_handle, "fasta")
            break
        current_normalizer = new_func
        if temp_in_path != input_fasta and os.path.exists(temp_in_path):
            os.remove(temp_in_path)
        temp_in_path = temp_out_path

    out_handle.close()
    if os.path.exists(temp_out_path):
        os.remove(temp_out_path)

def move_original_fasta(original_fasta):
    # Parent of the database folder containing the fasta
    input_dir = os.path.dirname(os.path.abspath(original_fasta))
    parent_dir = os.path.abspath(os.path.join(input_dir, ".."))  # one level up
    destination = os.path.join(parent_dir, os.path.basename(original_fasta))
    try:
        shutil.move(original_fasta, destination)
        print(f"[INFO] Moved original FASTA to: {destination}")
    except Exception as e:
        print(f"[WARNING] Could not move the FASTA file. Error: {e}")

def process_single_fasta(input_fasta, cdhit_exec):
    input_dir = os.path.dirname(os.path.abspath(input_fasta)) or "."
    prefix = os.path.splitext(os.path.basename(input_fasta))[0]
    output_fasta = os.path.join(input_dir, f"{prefix}_95percent.fasta")

    if is_already_processed(input_fasta):
        print(f"[INFO] File '{input_fasta}' already contains cluster annotations. Skipping.")
        return

    sorted_fasta_path = os.path.join(input_dir, f"{prefix}_sorted.fasta")
    records = list(SeqIO.parse(input_fasta, "fasta"))
    records.sort(key=lambda r: r.id)
    with open(sorted_fasta_path, "w") as sorted_fh:
        SeqIO.write(records, sorted_fh, "fasta")

    tmp_fa, clstr_path = run_cdhit_clustering(sorted_fasta_path, identity=0.95, cdhit_path=cdhit_exec)
    cluster_map = parse_clstr_file(clstr_path)

    annotate_fasta_with_clusters(input_fasta, cluster_map, output_fasta)

    # Cleanup temporary cd-hit files
    for f in [tmp_fa, clstr_path, sorted_fasta_path]:
        if os.path.exists(f):
            os.remove(f)

    print(f"[DONE] Clustered FASTA written to: {output_fasta}")

    move_original_fasta(input_fasta)

def main():
    parser = argparse.ArgumentParser(description="Cluster FASTA with cd-hit and annotate headers with cluster IDs.")
    parser.add_argument("-i", "--input", help="Input FASTA/FAA file or directory (default: all files in ./database/)", default=None)
    parser.add_argument("--cdhit-path", type=str, default=None, help="Path to cd-hit binary")
    args = parser.parse_args()

    cdhit_exec = resolve_cdhit_path(args.cdhit_path)

    # Properly handle input argument
    if args.input:
        input_path = os.path.abspath(args.input)
        if os.path.isfile(input_path):
            input_fastas = [input_path]
        elif os.path.isdir(input_path):
            # expand to all FASTA/FAA files in directory
            input_fastas = glob(os.path.join(input_path, "*.fasta")) + glob(os.path.join(input_path, "*.faa"))
            if not input_fastas:
                raise FileNotFoundError(f"No FASTA/FAA files found in {input_path}")
        else:
            raise FileNotFoundError(f"Input path not found: {input_path}")
    else:
        # default to ./database/
        input_fastas = find_input_fastas()

    for fasta in input_fastas:
        process_single_fasta(fasta, cdhit_exec)

if __name__ == "__main__":
    main()