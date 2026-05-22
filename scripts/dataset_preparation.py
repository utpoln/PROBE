"""
Dataset Preparation Pipeline for:
"Cross-Organism Evaluation of General-Purpose LLMs for GO-Based
Protein Function Prediction: A Zero-Shot Benchmark"

Steps:
  1. Parse uniprot_sprot.dat
  2. Filter to experimental GO evidence only
  3. Stratify by organism (5 organisms)
  4. Sample 1000 proteins per organism
  5. Split by GO namespace (BP, MF, CC)
  6. Export to CSV files ready for LLM prompting
  7. Print dataset statistics

Requirements:
  pip install biopython pandas
"""

import random
import json
from pathlib import Path
from collections import defaultdict

import pandas as pd
from Bio import SwissProt

# ─────────────────────────────────────────────
# CONFIG
# ─────────────────────────────────────────────

DAT_FILE = "uniprot_sprot.dat"       # <-- put your .dat file path here
OUTPUT_DIR = Path("benchmark_dataset")
RANDOM_SEED = 42
SAMPLE_PER_ORGANISM = 200
MIN_GO_TERMS = 3                      # minimum GO terms per entry

# Experimental evidence codes only (no computational)
EXPERIMENTAL_CODES = {"IDA", "IMP", "IPI", "IGI", "IEP", "EXP", "HDA", "HMP", "HGI", "HEP"}

# Target organisms: (Swiss-Prot organism string fragment, label)
TARGET_ORGANISMS = {
    "Homo sapiens":           "human",
    "Mus musculus":           "mouse",
    "Danio rerio":            "zebrafish",
    "Saccharomyces cerevisiae": "yeast",
    "Escherichia coli":       "ecoli",
}

# GO namespace codes
NAMESPACE_MAP = {
    "F": "molecular_function",
    "P": "biological_process",
    "C": "cellular_component",
}

# ─────────────────────────────────────────────
# STEP 1: PARSE + FILTER
# ─────────────────────────────────────────────

def parse_swissport(dat_file):
    """Parse Swiss-Prot .dat and return filtered records."""
    print(f"[1/5] Parsing {dat_file} ...")
    
    all_entries = []
    skipped = 0
    total = 0

    with open(dat_file, encoding="utf-8") as f:
        for record in SwissProt.parse(f):
            total += 1

            # Match organism
            organism_label = None
            for org_key, label in TARGET_ORGANISMS.items():
                if org_key in record.organism:
                    organism_label = label
                    break
            if organism_label is None:
                skipped += 1
                continue

            # Extract experimental GO terms only
            go_terms = []
            for ref in record.cross_references:
                if ref[0] != "GO":
                    continue
                go_id    = ref[1]                          # e.g. GO:0003700
                go_label = ref[2]                          # e.g. F:DNA binding
                evidence = ref[3].split(":")[0]            # e.g. IDA
                if evidence not in EXPERIMENTAL_CODES:
                    continue
                namespace_code = go_label[0]               # F, P, or C
                namespace = NAMESPACE_MAP.get(namespace_code, "unknown")
                term_name = go_label[2:]                   # strip "F:" prefix
                go_terms.append({
                    "go_id":      go_id,
                    "term_name":  term_name,
                    "namespace":  namespace,
                    "evidence":   evidence,
                })

            if len(go_terms) < MIN_GO_TERMS:
                skipped += 1
                continue

            # Extract function comment
            func_comment = ""
            for comment in record.comments:
                if comment.startswith("FUNCTION"):
                    func_comment = comment.replace("FUNCTION: ", "").strip()
                    break

            # Gene name — record.gene_name is a list of dicts in BioPython
            gene_name = ""
            if record.gene_name:
                try:
                    if isinstance(record.gene_name, list):
                        # e.g. [{'Name': ['TP53'], 'Synonyms': ['P53']}]
                        first = record.gene_name[0]
                        if isinstance(first, dict) and "Name" in first:
                            gene_name = first["Name"][0]
                        elif isinstance(first, str):
                            gene_name = first
                    elif isinstance(record.gene_name, str):
                        parts = record.gene_name.split(";")
                        for p in parts:
                            if p.strip().startswith("Name="):
                                gene_name = p.strip().replace("Name=", "").strip()
                                break
                except Exception:
                    gene_name = ""

            # Fallback to entry_name prefix if gene_name missing or too short
            if len(gene_name) <= 1:
                gene_name = record.entry_name.split("_")[0]
            all_entries.append({
                "accession":      record.accessions[0],
                "entry_name":     record.entry_name,
                "gene_name":      gene_name,
                "organism":       organism_label,
                "organism_full":  record.organism.strip(),
                "protein_name":   record.description,
                "function_text":  func_comment,
                "go_terms":       go_terms,
                "n_go_terms":     len(go_terms),
            })

    print(f"    Total records parsed : {total:,}")
    print(f"    After filtering      : {len(all_entries):,}")
    print(f"    Skipped              : {skipped:,}")
    return all_entries


# ─────────────────────────────────────────────
# STEP 2: STRATIFIED SAMPLING
# ─────────────────────────────────────────────

def stratified_sample(entries, n=SAMPLE_PER_ORGANISM, seed=RANDOM_SEED):
    """Sample N proteins per organism."""
    print(f"\n[2/5] Stratified sampling ({n} per organism) ...")
    
    by_organism = defaultdict(list)
    for e in entries:
        by_organism[e["organism"]].append(e)

    sampled = []
    random.seed(seed)
    for org, records in by_organism.items():
        if len(records) < n:
            print(f"    WARNING: {org} has only {len(records)} entries (< {n}), using all.")
            sampled.extend(records)
        else:
            sampled.extend(random.sample(records, n))
        print(f"    {org:12s}: {min(len(records), n):,} sampled  (pool: {len(records):,})")

    print(f"    Total sampled: {len(sampled):,}")
    return sampled


# ─────────────────────────────────────────────
# STEP 3: BUILD FLAT DATAFRAME
# ─────────────────────────────────────────────

def build_dataframe(sampled):
    """Flatten go_terms list into one row per protein with joined labels."""
    print(f"\n[3/5] Building dataframe ...")

    rows = []
    for e in sampled:
        # Group GO terms by namespace
        by_ns = defaultdict(list)
        for gt in e["go_terms"]:
            by_ns[gt["namespace"]].append(f"{gt['go_id']}|{gt['term_name']}")

        rows.append({
            "accession":          e["accession"],
            "entry_name":         e["entry_name"],
            "gene_name":          e["gene_name"],
            "organism":           e["organism"],
            "organism_full":      e["organism_full"],
            "protein_name":       e["protein_name"],
            "function_text":      e["function_text"],
            "n_go_terms_total":   e["n_go_terms"],
            # Ground truth per namespace (semicolon-separated)
            "go_molecular_function":  "; ".join(by_ns["molecular_function"]),
            "go_biological_process":  "; ".join(by_ns["biological_process"]),
            "go_cellular_component":  "; ".join(by_ns["cellular_component"]),
            # All GO terms as JSON string (for evaluation scripts)
            "go_terms_json":      json.dumps(e["go_terms"]),
        })

    df = pd.DataFrame(rows)
    print(f"    Dataframe shape: {df.shape}")
    return df


# ─────────────────────────────────────────────
# STEP 4: EXPORT
# ─────────────────────────────────────────────

def export(df, output_dir):
    """Save master CSV + per-organism + per-namespace splits."""
    output_dir.mkdir(exist_ok=True)
    print(f"\n[4/5] Exporting to {output_dir}/ ...")

    # Master file
    master_path = output_dir / "benchmark_master.csv"
    df.to_csv(master_path, index=False)
    print(f"    Saved: {master_path}  ({len(df):,} rows)")

    # Per-organism files
    org_dir = output_dir / "by_organism"
    org_dir.mkdir(exist_ok=True)
    for org in df["organism"].unique():
        subset = df[df["organism"] == org]
        path = org_dir / f"{org}.csv"
        subset.to_csv(path, index=False)
        print(f"    Saved: {path}  ({len(subset):,} rows)")

    # Per-namespace files (proteins that have at least 1 GO term in that namespace)
    ns_dir = output_dir / "by_namespace"
    ns_dir.mkdir(exist_ok=True)
    for ns_col, ns_name in [
        ("go_molecular_function", "MF"),
        ("go_biological_process", "BP"),
        ("go_cellular_component", "CC"),
    ]:
        subset = df[df[ns_col].str.len() > 0].copy()
        path = ns_dir / f"{ns_name}.csv"
        subset.to_csv(path, index=False)
        print(f"    Saved: {path}  ({len(subset):,} rows)")

    # Also save as JSON (useful for prompting scripts)
    json_path = output_dir / "benchmark_master.json"
    df.to_json(json_path, orient="records", indent=2)
    print(f"    Saved: {json_path}")


# ─────────────────────────────────────────────
# STEP 5: STATISTICS REPORT
# ─────────────────────────────────────────────

def print_stats(df):
    print(f"\n[5/5] Dataset Statistics")
    print("=" * 50)
    print(f"Total proteins          : {len(df):,}")
    print(f"\nBy organism:")
    print(df["organism"].value_counts().to_string())
    print(f"\nGO terms per protein (mean): {df['n_go_terms_total'].mean():.1f}")
    print(f"GO terms per protein (max) : {df['n_go_terms_total'].max()}")
    print(f"GO terms per protein (min) : {df['n_go_terms_total'].min()}")
    
    for ns_col, ns_name in [
        ("go_molecular_function", "Molecular Function (MF)"),
        ("go_biological_process", "Biological Process (BP)"),
        ("go_cellular_component", "Cellular Component (CC)"),
    ]:
        has_ns = df[df[ns_col].str.len() > 0]
        print(f"\nProteins with {ns_name}: {len(has_ns):,}")
    print("=" * 50)
    print("\n✅ Dataset preparation complete!")
    print(f"   Output folder: {OUTPUT_DIR}/")


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

if __name__ == "__main__":
    entries  = parse_swissport(DAT_FILE)
    sampled  = stratified_sample(entries)
    df       = build_dataframe(sampled)
    export(df, OUTPUT_DIR)
    print_stats(df)
