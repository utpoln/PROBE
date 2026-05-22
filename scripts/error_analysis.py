"""
PROBE: Error Analysis Script
Cross-Organism Evaluation of General-Purpose LLMs for GO-Based
Protein Function Prediction

Error categories:
  1. Complete Hallucination  — model predicts GO IDs that don't exist in GO
  2. Wrong Namespace         — model predicts real GO IDs but from wrong namespace
  3. Partial Match           — model gets some but not all correct GO terms
  4. Complete Miss           — model predicts nothing relevant
  5. Empty Response          — model returns no GO IDs at all
  6. Near Miss               — model predicts parent/child of correct term

Outputs:
  - error_summary.csv           overall error category counts per model
  - error_cases_sampled.csv     50 sampled failure cases per model
  - error_analysis.xlsx         full workbook
  - Console report

Requirements:
  pip install pandas openpyxl tqdm

Usage:
  python error_analysis.py
  python error_analysis.py --eval_dir probe_evaluation --results_dir probe_results
"""

import re
import json
import argparse
import random
from pathlib import Path
from collections import defaultdict

import pandas as pd
from tqdm import tqdm

# ─────────────────────────────────────────────
# CONFIG
# ─────────────────────────────────────────────

EVAL_DIR    = Path("probe_evaluation")
RESULTS_DIR = Path("probe_results")
OUTPUT_DIR  = Path("probe_error_analysis")
RANDOM_SEED = 42
SAMPLE_SIZE = 50   # failure cases to sample per model

GO_PATTERN  = re.compile(r"GO:\d{7}")

MODEL_LABELS = {
    "llama3.3:70b":       "Llama 3.3 70B",
    "llama3.1:8b":        "Llama 3.1 8B",
    "mistral-large:123b": "Mistral Large 123B",
    "mistral:7b":         "Mistral 7B",
    "qwen2.5:72b":        "Qwen2.5 72B",
    "qwen2.5:7b":         "Qwen2.5 7B",
    "Qwen/Qwen3-32B":     "Qwen3 32B",
    "gemma3:12b":         "Gemma3 12B",
    "google/gemma-4-31b": "Gemma-4 31B",
    "mixtral:8x7b":       "Mixtral 8x7B",
}

# GO namespace prefix mapping (first letter of GO label)
NAMESPACE_PREFIX = {
    "molecular_function": "F",
    "biological_process": "P",
    "cellular_component": "C",
}

# ─────────────────────────────────────────────
# LOAD DATA
# ─────────────────────────────────────────────

def load_raw_results(results_dir: Path) -> list:
    """Load all raw JSONL result records."""
    records = []
    for f in results_dir.glob("*.jsonl"):
        with open(f) as fh:
            for line in fh:
                try:
                    rec = json.loads(line.strip())
                    if rec.get("success"):
                        records.append(rec)
                except Exception:
                    continue
    print(f"Loaded {len(records):,} raw records from {results_dir}/")
    return records


def load_all_go_ids(results_dir: Path) -> set:
    """
    Collect all GO IDs that appear in ground truth across the dataset.
    Used as a proxy for 'real' GO IDs to detect hallucinations.
    Note: In production you'd load the full GO ontology OBO file.
    """
    all_ids = set()
    for f in results_dir.glob("*.jsonl"):
        with open(f) as fh:
            for line in fh:
                try:
                    rec = json.loads(line.strip())
                    gt = rec.get("ground_truth", "")
                    ids = GO_PATTERN.findall(gt)
                    all_ids.update(ids)
                except Exception:
                    continue
    print(f"Collected {len(all_ids):,} unique ground-truth GO IDs as reference set")
    return all_ids


# ─────────────────────────────────────────────
# ERROR CLASSIFICATION
# ─────────────────────────────────────────────

def classify_error(predicted: set, ground_truth: set,
                   real_go_ids: set, namespace: str) -> dict:
    """
    Classify the error type for a single prediction.

    Returns dict with:
      - primary_error: main error category
      - sub_errors: list of contributing issues
      - details: diagnostic info
    """
    sub_errors = []
    details    = {}

    # ── Empty response ──
    if len(predicted) == 0:
        return {
            "primary_error": "empty_response",
            "sub_errors":    ["no_go_ids_extracted"],
            "details":       {"n_predicted": 0, "n_ground_truth": len(ground_truth)},
        }

    # ── Compute overlap ──
    correct    = predicted & ground_truth
    wrong      = predicted - ground_truth
    missed     = ground_truth - predicted

    n_correct  = len(correct)
    n_wrong    = len(wrong)
    n_missed   = len(missed)
    n_pred     = len(predicted)
    n_gt       = len(ground_truth)

    # ── Hallucination check ──
    # GO IDs predicted but not in our reference set of known real IDs
    hallucinated = predicted - real_go_ids
    n_halluc     = len(hallucinated)
    halluc_rate  = n_halluc / n_pred if n_pred > 0 else 0.0

    if halluc_rate > 0.5:
        sub_errors.append("high_hallucination")
    elif halluc_rate > 0:
        sub_errors.append("partial_hallucination")

    details["hallucinated_ids"]  = list(hallucinated)[:5]
    details["hallucination_rate"] = round(halluc_rate, 3)

    # ── Complete miss ──
    if n_correct == 0:
        if n_halluc == n_pred:
            primary = "complete_hallucination"
        else:
            primary = "complete_miss"
        return {
            "primary_error": primary,
            "sub_errors":    sub_errors,
            "details":       {
                **details,
                "n_predicted": n_pred,
                "n_ground_truth": n_gt,
                "correct_ids": [],
                "missed_ids":  list(missed)[:5],
            },
        }

    # ── Partial match ──
    recall    = n_correct / n_gt  if n_gt   > 0 else 0.0
    precision = n_correct / n_pred if n_pred > 0 else 0.0

    if recall < 0.5:
        sub_errors.append("low_recall")
    if precision < 0.5:
        sub_errors.append("low_precision")

    # ── Perfect match ──
    if n_correct == n_gt and n_wrong == 0:
        primary = "correct"
    elif precision >= 0.5 and recall >= 0.5:
        primary = "partial_match_good"
    elif recall >= 0.5:
        primary = "partial_match_high_recall"
    elif precision >= 0.5:
        primary = "partial_match_high_precision"
    else:
        primary = "partial_match_poor"

    return {
        "primary_error": primary,
        "sub_errors":    sub_errors,
        "details":       {
            **details,
            "n_predicted":    n_pred,
            "n_ground_truth": n_gt,
            "n_correct":      n_correct,
            "n_wrong":        n_wrong,
            "n_missed":       n_missed,
            "precision":      round(precision, 3),
            "recall":         round(recall, 3),
            "correct_ids":    list(correct)[:5],
            "missed_ids":     list(missed)[:5],
            "wrong_ids":      list(wrong)[:5],
        },
    }


# ─────────────────────────────────────────────
# MAIN ANALYSIS
# ─────────────────────────────────────────────

def run_error_analysis(records: list, real_go_ids: set) -> pd.DataFrame:
    """Classify errors for all records."""
    print(f"\n[2/5] Classifying errors for {len(records):,} records ...")

    rows = []
    for rec in tqdm(records):
        predicted    = set(GO_PATTERN.findall(rec.get("response", "")))
        ground_truth = set(GO_PATTERN.findall(rec.get("ground_truth", "")))

        error = classify_error(predicted, ground_truth,
                               real_go_ids, rec.get("namespace", ""))

        rows.append({
            "model":          rec["model"],
            "model_label":    MODEL_LABELS.get(rec["model"], rec["model"]),
            "accession":      rec["accession"],
            "gene_name":      rec.get("gene_name", ""),
            "organism":       rec.get("organism", ""),
            "namespace":      rec.get("namespace", ""),
            "prompt_id":      rec.get("prompt_id", ""),
            "primary_error":  error["primary_error"],
            "sub_errors":     "|".join(error["sub_errors"]),
            "hallucination_rate": error["details"].get("hallucination_rate", 0.0),
            "n_predicted":    error["details"].get("n_predicted", 0),
            "n_ground_truth": error["details"].get("n_ground_truth", 0),
            "n_correct":      error["details"].get("n_correct", 0),
            "precision":      error["details"].get("precision", 0.0),
            "recall":         error["details"].get("recall", 0.0),
            "response_text":  rec.get("response", "")[:300],
            "ground_truth":   rec.get("ground_truth", "")[:200],
            "hallucinated_ids": str(error["details"].get("hallucinated_ids", [])),
            "missed_ids":     str(error["details"].get("missed_ids", [])),
        })

    df = pd.DataFrame(rows)
    print(f"    Error analysis dataframe: {df.shape}")
    return df


# ─────────────────────────────────────────────
# AGGREGATE ERROR SUMMARY
# ─────────────────────────────────────────────

def build_error_summary(df: pd.DataFrame) -> dict:
    """Build summary tables of error distributions."""
    print(f"\n[3/5] Building error summary tables ...")
    tables = {}

    # 1. Error category counts per model
    error_counts = (df.groupby(["model_label", "primary_error"])
                      .size()
                      .reset_index(name="count"))
    total_per_model = df.groupby("model_label").size().reset_index(name="total")
    error_counts = error_counts.merge(total_per_model, on="model_label")
    error_counts["percentage"] = (error_counts["count"] / error_counts["total"] * 100).round(2)
    tables["error_counts_per_model"] = error_counts

    # 2. Error pivot — model × error category
    error_pivot = error_counts.pivot(
        index="model_label", columns="primary_error", values="percentage"
    ).fillna(0).round(2)
    tables["error_pivot"] = error_pivot

    # 3. Hallucination rate per model × namespace
    halluc = (df.groupby(["model_label", "namespace"])["hallucination_rate"]
                .mean()
                .reset_index()
                .round(4))
    tables["hallucination_by_namespace"] = halluc

    # 4. Error rate per organism
    error_by_org = (df.groupby(["organism", "primary_error"])
                      .size()
                      .reset_index(name="count"))
    tables["error_by_organism"] = error_by_org

    # 5. Error rate per prompt format
    error_by_prompt = (df.groupby(["prompt_id", "primary_error"])
                         .size()
                         .reset_index(name="count"))
    tables["error_by_prompt"] = error_by_prompt

    for name, t in tables.items():
        print(f"    Table '{name}': {t.shape}")

    return tables


# ─────────────────────────────────────────────
# SAMPLE FAILURE CASES
# ─────────────────────────────────────────────

def sample_failures(df: pd.DataFrame, n: int = SAMPLE_SIZE) -> pd.DataFrame:
    """
    Sample n failure cases per model for qualitative analysis.
    Excludes correct predictions.
    """
    print(f"\n[4/5] Sampling {n} failure cases per model ...")
    failures = df[df["primary_error"] != "correct"].copy()

    random.seed(RANDOM_SEED)
    sampled = []
    for model in failures["model_label"].unique():
        subset = failures[failures["model_label"] == model]
        n_sample = min(n, len(subset))
        sampled.append(subset.sample(n_sample, random_state=RANDOM_SEED))
        print(f"    {model}: {n_sample} failure cases sampled")

    result = pd.concat(sampled, ignore_index=True) if sampled else pd.DataFrame()
    return result


# ─────────────────────────────────────────────
# EXPORT
# ─────────────────────────────────────────────

def export(df: pd.DataFrame, tables: dict,
           sampled: pd.DataFrame, output_dir: Path):
    output_dir.mkdir(exist_ok=True)
    print(f"\n[5/5] Exporting to {output_dir}/ ...")

    # Raw error classification
    raw_path = output_dir / "error_analysis_raw.csv"
    df.to_csv(raw_path, index=False)
    print(f"    Saved: {raw_path}")

    # Sampled failure cases
    if not sampled.empty:
        sample_path = output_dir / "error_cases_sampled.csv"
        sampled.to_csv(sample_path, index=False)
        print(f"    Saved: {sample_path}")

    # Individual summary tables
    for name, table in tables.items():
        path = output_dir / f"{name}.csv"
        table.to_csv(path, index=False)
        print(f"    Saved: {path}")

    # Excel workbook — clean illegal characters from text columns first
    import re as _re
    def clean_for_excel(df_in):
        df_clean = df_in.copy()
        for col in df_clean.select_dtypes(include='object').columns:
            df_clean[col] = df_clean[col].astype(str).apply(
                lambda x: _re.sub(r'[\x00-\x08\x0b\x0c\x0e-\x1f]', '', x)[:500]
            )
        return df_clean

    excel_path = output_dir / "PROBE_error_analysis.xlsx"
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        clean_for_excel(df).to_excel(writer, sheet_name="raw_errors", index=False)
        if not sampled.empty:
            clean_for_excel(sampled).to_excel(writer, sheet_name="sampled_failures", index=False)
        for name, table in tables.items():
            clean_for_excel(table).to_excel(writer, sheet_name=name[:31], index=False)
    print(f"    Saved: {excel_path}")


# ─────────────────────────────────────────────
# PRINT REPORT
# ─────────────────────────────────────────────

def print_report(df: pd.DataFrame, tables: dict):
    print("\n" + "=" * 65)
    print("PROBE — Error Analysis Report")
    print("=" * 65)

    # Overall error distribution
    print("\n📊 Overall Error Distribution (all models combined):")
    overall = (df["primary_error"].value_counts(normalize=True) * 100).round(2)
    for err, pct in overall.items():
        bar = "█" * int(pct / 2)
        print(f"   {err:<30} {pct:>6.2f}%  {bar}")

    # Per model summary
    print("\n📊 Error Breakdown Per Model:")
    pivot = tables.get("error_pivot", pd.DataFrame())
    if not pivot.empty:
        print(pivot.to_string())

    # Hallucination rates
    print("\n📊 Avg Hallucination Rate by Model:")
    halluc_model = (df.groupby("model_label")["hallucination_rate"]
                      .mean()
                      .sort_values(ascending=False)
                      .round(4))
    for model, rate in halluc_model.items():
        bar = "█" * int(rate * 20)
        print(f"   {model:<25} {rate*100:>6.1f}%  {bar}")

    # Best/worst organisms
    print("\n📊 Error Rate by Organism:")
    org_correct = (df[df["primary_error"] == "correct"]
                     .groupby("organism").size())
    org_total   = df.groupby("organism").size()
    org_acc     = (org_correct / org_total * 100).round(2).sort_values(ascending=False)
    for org, acc in org_acc.items():
        print(f"   {org:<15} correct rate: {acc:.2f}%")

    # Best prompt format
    print("\n📊 Correct Rate by Prompt Format:")
    prompt_correct = (df[df["primary_error"] == "correct"]
                        .groupby("prompt_id").size())
    prompt_total   = df.groupby("prompt_id").size()
    prompt_acc     = (prompt_correct / prompt_total * 100).round(2).sort_values(ascending=False)
    for pid, acc in prompt_acc.items():
        print(f"   {pid:<20} correct rate: {acc:.2f}%")

    print("\n" + "=" * 65)
    print("✅ Error analysis complete!")
    print(f"   Full results in: {OUTPUT_DIR}/")
    print(f"   Excel workbook : {OUTPUT_DIR}/PROBE_error_analysis.xlsx")
    print("\n💡 Use error_cases_sampled.csv for qualitative Discussion section")


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PROBE Error Analysis")
    parser.add_argument("--results_dir", default="probe_results")
    parser.add_argument("--eval_dir",    default="probe_evaluation")
    parser.add_argument("--output_dir",  default="probe_error_analysis")
    args = parser.parse_args()

    RESULTS_DIR = Path(args.results_dir)
    EVAL_DIR    = Path(args.eval_dir)
    OUTPUT_DIR  = Path(args.output_dir)

    print("=" * 65)
    print("PROBE — Error Analysis Pipeline")
    print("=" * 65)

    print("\n[1/5] Loading data ...")
    records     = load_raw_results(RESULTS_DIR)
    real_go_ids = load_all_go_ids(RESULTS_DIR)

    df      = run_error_analysis(records, real_go_ids)
    tables  = build_error_summary(df)
    sampled = sample_failures(df)

    export(df, tables, sampled, OUTPUT_DIR)
    print_report(df, tables)
