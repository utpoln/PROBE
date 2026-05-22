"""
PROBE: Evaluation Script
Cross-Organism Evaluation of General-Purpose LLMs for GO-Based
Protein Function Prediction

Computes:
  - Precision, Recall, F1 per model/organism/namespace/prompt
  - Hallucination Rate: predicted GO IDs absent from the FULL
    reference set (union of all ground truth GO IDs in the benchmark)
  - Empty Response Rate
  - Exports to CSV + Excel

Usage:
  python evaluation.py
  python evaluation.py --results_dir probe_results --output_dir probe_evaluation
"""

import re
import json
import argparse
from pathlib import Path

import pandas as pd
from tqdm import tqdm

# ── Config ────────────────────────────────────────────────────
RESULTS_DIR = Path("probe_results")
OUTPUT_DIR  = Path("probe_evaluation")
GO_PATTERN  = re.compile(r"GO:\d{7}")

NS_SHORT = {
    "molecular_function": "MF",
    "biological_process": "BP",
    "cellular_component": "CC",
}

# Map any known filename variants to a clean display name
MODEL_DISPLAY = {
    "llama3.1:8b":          "Llama 3.1 8B",
    "llama3.3:70b":         "Llama 3.3 70B",
    "mistral:7b":           "Mistral 7B",
    "mistral-large:123b":   "Mistral Large 123B",
    "qwen2.5:7b":           "Qwen2.5 7B",
    "qwen2.5:72b":          "Qwen2.5 72B",
    "Qwen/Qwen3-32B":       "Qwen3 32B",
    "gemma3:12b":           "Gemma3 12B",
    "phi4:14b":             "Phi-4 14B",
    "microsoft/phi-4":      "Phi-4 14B",
    "mixtral:8x7b":         "Mixtral 8x7B",
}

# ── Step 1: Load ──────────────────────────────────────────────
def load_results(results_dir: Path) -> list:
    records = []
    files = sorted(results_dir.glob("*.jsonl"))
    if not files:
        print(f"No .jsonl files in {results_dir}/"); exit(1)

    print(f"[1/5] Loading results from {len(files)} model file(s) ...")
    for f in files:
        count = 0
        with open(f) as fh:
            for line in fh:
                try:
                    rec = json.loads(line.strip())
                    if rec.get("success"):
                        # Normalise model name to display name
                        raw_model = rec.get("model", f.stem)
                        rec["model"] = MODEL_DISPLAY.get(raw_model, raw_model)
                        records.append(rec)
                        count += 1
                except Exception:
                    continue
        print(f"    {f.name}: {count:,} records")

    print(f"    Total records loaded: {len(records):,}")
    return records

# ── Step 2: Build reference set ───────────────────────────────
def build_reference_set(records: list) -> set:
    """
    Union of ALL ground truth GO IDs across the entire benchmark.
    A predicted GO ID absent from this set is a hallucination.
    """
    ref = set()
    for rec in records:
        ref.update(GO_PATTERN.findall(rec.get("ground_truth", "")))
    return ref

# ── Step 3: Parse ─────────────────────────────────────────────
def parse_go(text: str) -> set:
    if not text: return set()
    return set(GO_PATTERN.findall(text))

def compute_metrics(predicted: set, ground_truth: set) -> dict:
    if not ground_truth:
        return {"precision": 0.0, "recall": 0.0, "f1": 0.0,
                "tp": 0, "fp": 0, "fn": 0}
    tp = len(predicted & ground_truth)
    fp = len(predicted - ground_truth)
    fn = len(ground_truth - predicted)
    prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    rec  = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1   = 2*prec*rec / (prec+rec) if (prec+rec) > 0 else 0.0
    return {"precision": round(prec,4), "recall": round(rec,4),
            "f1": round(f1,4), "tp": tp, "fp": fp, "fn": fn}

# ── Step 4: Evaluate ──────────────────────────────────────────
def build_eval_df(records: list, ref_set: set) -> pd.DataFrame:
    print(f"\n[2/5] Evaluating {len(records):,} records ...")
    rows = []
    for rec in tqdm(records):
        predicted    = parse_go(rec["response"])
        ground_truth = parse_go(rec["ground_truth"])
        metrics      = compute_metrics(predicted, ground_truth)

        n_pred = len(predicted)

        # Hallucination: predicted IDs absent from the FULL reference set
        hallucinated = predicted - ref_set
        hall_rate = len(hallucinated) / n_pred if n_pred > 0 else 0.0

        rows.append({
            "model":             rec["model"],
            "accession":         rec["accession"],
            "gene_name":         rec.get("gene_name", ""),
            "organism":          rec["organism"],
            "namespace":         rec["namespace"],
            "namespace_short":   NS_SHORT.get(rec["namespace"], rec["namespace"]),
            "prompt_id":         rec["prompt_id"],
            "n_predicted":       n_pred,
            "n_ground_truth":    len(ground_truth),
            "tp":                metrics["tp"],
            "fp":                metrics["fp"],
            "fn":                metrics["fn"],
            "precision":         metrics["precision"],
            "recall":            metrics["recall"],
            "f1":                metrics["f1"],
            # CORRECT hallucination: absent from full reference set
            "hallucination_rate": round(hall_rate, 4),
            "empty_response":    int(n_pred == 0),
        })

    df = pd.DataFrame(rows)
    print(f"    Evaluation dataframe: {df.shape}")
    return df

# ── Step 5: Aggregate ─────────────────────────────────────────
def aggregate(df: pd.DataFrame) -> dict:
    print(f"\n[3/5] Aggregating results ...")
    tables = {}
    agg = {
        "precision": "mean", "recall": "mean", "f1": "mean",
        "hallucination_rate": "mean", "empty_response": "mean",
        "n_predicted": "mean", "n_ground_truth": "mean",
    }
    def fmt(d): return d.round(4)

    tables["leaderboard"] = fmt(
        df.groupby("model").agg(agg).reset_index()
          .sort_values("f1", ascending=False).reset_index(drop=True))

    tables["model_x_namespace"] = fmt(
        df.groupby(["model","namespace_short"]).agg(agg).reset_index())

    tables["model_x_organism"] = fmt(
        df.groupby(["model","organism"]).agg(agg).reset_index())

    tables["model_x_prompt"] = fmt(
        df.groupby(["model","prompt_id"]).agg(agg).reset_index())

    tables["organism_x_namespace"] = fmt(
        df.groupby(["organism","namespace_short"]).agg(agg).reset_index())

    tables["full_detail"] = fmt(
        df.groupby(["model","organism","namespace_short","prompt_id"])
          .agg(agg).reset_index())

    for name, t in tables.items():
        print(f"    Table '{name}': {t.shape}")
    return tables

# ── Step 6: Export ────────────────────────────────────────────
def export(df_raw: pd.DataFrame, tables: dict, output_dir: Path):
    output_dir.mkdir(exist_ok=True)
    print(f"\n[4/5] Exporting to {output_dir}/ ...")
    df_raw.to_csv(output_dir/"evaluation_raw.csv", index=False)
    print(f"    Saved: evaluation_raw.csv")
    for name, table in tables.items():
        table.to_csv(output_dir/f"{name}.csv", index=False)
        print(f"    Saved: {name}.csv")
    excel_path = output_dir/"PROBE_results.xlsx"
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        df_raw.to_excel(writer, sheet_name="raw", index=False)
        for name, table in tables.items():
            table.to_excel(writer, sheet_name=name[:31], index=False)
    print(f"    Saved: PROBE_results.xlsx")

# ── Step 7: Print leaderboard ─────────────────────────────────
def print_leaderboard(tables: dict):
    print(f"\n[5/5] PROBE Leaderboard")
    print("=" * 95)
    lb = tables["leaderboard"]
    print(f"{'Rank':<5} {'Model':<25} {'F1':>7} {'Prec':>7} "
          f"{'Recall':>7} {'Halluc%':>9} {'Empty%':>8}")
    print("-" * 95)
    for rank, (_, row) in enumerate(lb.iterrows(), 1):
        print(f"{rank:<5} {row['model']:<25} "
              f"{row['f1']:>7.4f} "
              f"{row['precision']:>7.4f} "
              f"{row['recall']:>7.4f} "
              f"{row['hallucination_rate']*100:>8.1f}% "
              f"{row['empty_response']*100:>7.1f}%")
    print("=" * 95)

    print(f"\n  Best model : {lb.iloc[0]['model']} (F1={lb.iloc[0]['f1']:.4f})")
    print(f"  Worst model: {lb.iloc[-1]['model']} (F1={lb.iloc[-1]['f1']:.4f})")

    print(f"\n--- Avg F1 by Organism ---")
    org = tables["model_x_organism"].groupby("organism")["f1"].mean().sort_values(ascending=False)
    for o, f in org.items():
        print(f"  {o:<15}: {f:.4f}")

    print(f"\n--- Avg F1 by Prompt ---")
    prom = tables["model_x_prompt"].groupby("prompt_id")["f1"].mean().sort_values(ascending=False)
    for p, f in prom.items():
        print(f"  {p:<20}: {f:.4f}")

    print(f"\n--- Avg F1 by Namespace ---")
    ns = tables["model_x_namespace"].groupby("namespace_short")["f1"].mean().sort_values(ascending=False)
    for n, f in ns.items():
        print(f"  {n}: {f:.4f}")

    print("\n✅ Evaluation complete!")
    print(f"   Results: {OUTPUT_DIR}/PROBE_results.xlsx")

# ── Main ──────────────────────────────────────────────────────
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--results_dir", default="probe_results")
    parser.add_argument("--output_dir",  default="probe_evaluation")
    args = parser.parse_args()

    RESULTS_DIR = Path(args.results_dir)
    OUTPUT_DIR  = Path(args.output_dir)

    records  = load_results(RESULTS_DIR)
    ref_set  = build_reference_set(records)
    print(f"    Reference GO set: {len(ref_set):,} unique experimental terms")
    df_raw   = build_eval_df(records, ref_set)
    tables   = aggregate(df_raw)
    export(df_raw, tables, OUTPUT_DIR)
    print_leaderboard(tables)