"""
PROBE: Open-Source Reference Row for the Closed-Source Comparison Table

closed_source_eval.py scores closed-source models by averaging per-instance
F1 over every protein x namespace pair in its 50-protein sample, scoring 0
for pairs with no experimental ground truth in that namespace (the same
convention evaluation.py uses for the main benchmark: see compute_metrics()
in both scripts). This script applies that identical convention to compute
the open-source reference row -- Mistral Large 123B under P5 candidate
selection across the full 1,000-protein benchmark -- so the comparison
table (paper Table 14) is computed under one consistent methodology
end-to-end.

Usage:
    python closed_source_reference.py --model mistral-large-123b

Output:
    probe_evaluation/closed_source_reference.csv
"""

import argparse
import csv
import json
import re
from pathlib import Path

RESULTS_DIR = Path("probe_results")
CLOSED_RESULTS_DIR = Path("probe_results_closed")
OUTPUT_DIR = Path("probe_evaluation")
GO_PATTERN = re.compile(r"GO:\d{7}")

MODEL_LABELS = {
    "claude-sonnet-4-6": "Claude Sonnet 4.6",
    "gemini-2-5-pro": "Gemini 2.5 Pro",
    "gpt-5-4-2026-03-05": "GPT-5.4",
    "gemini-2-5-flash": "Gemini 2.5 Flash",
}


def compute_metrics(predicted: set, ground_truth: set):
    """Identical convention to evaluation.py: scores 0/0/0 whenever the
    protein/namespace pair has no experimental ground truth, and includes
    every such instance in the macro-average (rather than skipping it)."""
    if not ground_truth:
        return 0.0, 0.0, 0.0
    tp = len(predicted & ground_truth)
    fp = len(predicted - ground_truth)
    fn = len(ground_truth - predicted)
    prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    rec = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0.0
    return prec, rec, f1


def macro_average(filepath: Path, prompt_filter: str = None):
    precs, recs, f1s = [], [], []
    with open(filepath) as f:
        for line in f:
            try:
                rec = json.loads(line)
            except Exception:
                continue
            if not rec.get("success", False):
                continue
            if prompt_filter and rec.get("prompt_id") != prompt_filter:
                continue
            predicted = set(GO_PATTERN.findall(rec.get("response", "") or ""))
            ground_truth = set(GO_PATTERN.findall(rec.get("ground_truth", "") or ""))
            p, r, f1 = compute_metrics(predicted, ground_truth)
            precs.append(p)
            recs.append(r)
            f1s.append(f1)
    n = len(f1s)
    if n == 0:
        return None
    return {
        "n": n,
        "precision": round(sum(precs) / n, 4),
        "recall": round(sum(recs) / n, 4),
        "f1": round(sum(f1s) / n, 4),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", default="mistral-large-123b",
                         help="Open-source results filename stem in probe_results/ "
                              "(without .jsonl); default matches the paper's reference model.")
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(exist_ok=True)

    open_path = RESULTS_DIR / f"{args.model}.jsonl"
    if not open_path.exists():
        # fall back to the exact filename used in this repo
        open_path = RESULTS_DIR / "mistral-large-123b.jsonl"
    open_metrics = macro_average(open_path, prompt_filter="P5_selection")
    if open_metrics is None:
        raise SystemExit(f"No P5 instances found in {open_path}")

    rows = [{
        "model": "Mistral Large 123B",
        "type": f"Open, {1000}-protein reference",
        "n": open_metrics["n"],
        "precision": open_metrics["precision"],
        "recall": open_metrics["recall"],
        "f1": open_metrics["f1"],
    }]

    for stem, label in MODEL_LABELS.items():
        fp = CLOSED_RESULTS_DIR / f"{stem}.jsonl"
        if not fp.exists():
            continue
        m = macro_average(fp)
        if m is None:
            continue
        rows.append({
            "model": label,
            "type": "Closed (50-protein sample)",
            "n": m["n"],
            "precision": m["precision"],
            "recall": m["recall"],
            "f1": m["f1"],
        })

    rows.sort(key=lambda r: (-r["f1"] if r["type"].startswith("Closed") else -999))

    with open(OUTPUT_DIR / "closed_source_reference.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["model", "type", "n", "precision", "recall", "f1"])
        w.writeheader()
        w.writerows(rows)

    print("=== Closed-source vs. open-source reference (Table 14, one consistent convention) ===")
    for r in rows:
        print(f"  {r['model']:<20s} {r['type']:<32s} N={r['n']:5d}  "
              f"Prec={r['precision']:.3f}  Rec={r['recall']:.3f}  F1={r['f1']:.3f}")
    closed_f1s = [r["f1"] for r in rows if r["type"].startswith("Closed")]
    open_f1 = rows[0]["f1"]
    if closed_f1s:
        best_closed = max(closed_f1s)
        print(f"\nBest closed-source F1: {best_closed:.3f}; open-source reference F1: {open_f1:.3f}; "
              f"gap = {best_closed - open_f1:.3f}")
    print(f"\nSaved: {OUTPUT_DIR / 'closed_source_reference.csv'}")


if __name__ == "__main__":
    main()
