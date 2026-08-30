"""
PROBE: Temporal-Holdout Validation of Unmatched Predictions

Motivation: GO and UniProt/Swiss-Prot annotations are incomplete and
continuously updated. A predicted GO term absent from PROBE's benchmark
reference set R (built from the Swiss-Prot 2024_03 release) is not
automatically biologically incorrect -- it could be a true annotation
that simply had not yet been curated as of that release.

This script tests that possibility directly. It samples "out-of-benchmark"
and "misattribution" predictions (the two GO/Swiss-Prot-incompleteness-
sensitive categories from hallucination_decompose.py) and checks whether
the predicted term now appears among the target protein's annotations in
the CURRENT UniProt release (i.e., a release strictly newer than the one
used to build the benchmark). If the large majority of sampled predictions
remain unsupported under a substantially more recent and complete
reference, that is evidence the observed unmatched-prediction rate mostly
reflects genuine errors rather than annotation-database incompleteness.

Usage:
    python temporal_holdout_validation.py [--n-per-category 120] [--seed 42]

Requires network access to https://rest.uniprot.org.

Output:
    probe_statistics/temporal_holdout_sample.csv       (per-prediction detail)
    probe_statistics/temporal_holdout_summary.csv       (aggregate statistics)
"""

import argparse
import csv
import json
import random
import re
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

from scipy import stats

RESULTS_DIR = Path("probe_results")
OBO_PATH = Path("go-basic.obo")
OUTPUT_DIR = Path("probe_statistics")
GO_PATTERN = re.compile(r"GO:\d{7}")

EXPERIMENTAL_CODES = {
    "IDA", "IMP", "IPI", "IGI", "IEP", "EXP", "HDA", "HMP", "HGI", "HEP",
}


def load_valid_go_ids(obo_path: Path) -> set:
    valid = set()
    with open(obo_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("id: GO:"):
                valid.add(line.replace("id: ", "").strip())
    return valid


def build_reference_and_ground_truth(results_dir: Path):
    """Union of all ground-truth GO IDs (R), per-(accession, namespace)
    ground truth sets, and an accession -> organism lookup."""
    R = set()
    protein_gt = defaultdict(set)
    acc_organism = {}
    for filepath in results_dir.glob("*.jsonl"):
        with open(filepath) as f:
            for line in f:
                try:
                    rec = json.loads(line)
                except Exception:
                    continue
                gt = rec.get("ground_truth", "")
                if not gt:
                    continue
                terms = set(GO_PATTERN.findall(gt))
                key = (rec["accession"], rec["namespace"])
                protein_gt[key].update(terms)
                R.update(terms)
                acc_organism[rec["accession"]] = rec.get("organism", "")
    return R, protein_gt, acc_organism


def collect_unmatched_tuples(results_dir: Path, valid_go_ids: set, R: set, protein_gt: dict):
    """Unique (accession, namespace, go_id) tuples for the two
    incompleteness-sensitive failure categories: out-of-benchmark
    (valid GO ID, absent from R) and misattribution (valid GO ID,
    present in R but not the correct annotation for this protein)."""
    out_of_benchmark, misattribution = {}, {}
    for filepath in results_dir.glob("*.jsonl"):
        with open(filepath) as f:
            for line in f:
                try:
                    rec = json.loads(line)
                except Exception:
                    continue
                if not rec.get("success", False):
                    continue
                key = (rec["accession"], rec["namespace"])
                gt_terms = protein_gt.get(key, set())
                if not gt_terms:
                    continue
                predicted = set(GO_PATTERN.findall(rec.get("response", "") or ""))
                for go_id in predicted:
                    if go_id in gt_terms or go_id not in valid_go_ids:
                        continue  # correct, or syntactic (not incompleteness-sensitive)
                    tup = (rec["accession"], rec["namespace"], go_id)
                    if go_id not in R:
                        out_of_benchmark[tup] = True
                    else:
                        misattribution[tup] = True
    return out_of_benchmark, misattribution


def fetch_current_uniprot_go(accessions, sleep_sec: float = 0.15):
    """Fetch current GO cross-references (any evidence code) for each
    accession from the live UniProt REST API."""
    cache = {}
    for i, acc in enumerate(sorted(accessions)):
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
        for attempt in range(3):
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "PROBE-holdout-check/1.0"})
                with urllib.request.urlopen(req, timeout=20) as resp:
                    data = json.loads(resp.read().decode())
                entries = []
                for ref in data.get("uniProtKBCrossReferences", []):
                    if ref.get("database") != "GO":
                        continue
                    evtype = ""
                    for prop in ref.get("properties", []):
                        if prop.get("key") == "GoEvidenceType":
                            evtype = prop.get("value", "")
                    entries.append({"go_id": ref.get("id"), "evidence": evtype})
                cache[acc] = entries
                break
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    cache[acc] = None
                    break
                time.sleep(1.5)
            except Exception:
                time.sleep(1.5)
        else:
            cache[acc] = None
        if (i + 1) % 25 == 0:
            print(f"  fetched {i + 1}/{len(accessions)} accessions...")
        time.sleep(sleep_sec)
    return cache


def classify(sample_rows, uniprot_cache):
    out = []
    for acc, ns, go_id, category in sample_rows:
        entries = uniprot_cache.get(acc)
        status = "protein_not_found" if entries is None else "absent_current"
        evidence = ""
        if entries is not None:
            for e in entries:
                if e["go_id"] == go_id:
                    evidence = e.get("evidence", "")
                    code = evidence.split(":")[0] if evidence else ""
                    status = "present_experimental" if code in EXPERIMENTAL_CODES else "present_non_experimental"
                    break
        out.append((acc, ns, go_id, category, status, evidence))
    return out


def clopper_pearson(k: int, n: int, alpha: float = 0.05):
    lo = 0.0 if k == 0 else stats.beta.ppf(alpha / 2, k, n - k + 1)
    hi = 1.0 if k == n else stats.beta.ppf(1 - alpha / 2, k + 1, n - k)
    return lo, hi


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-per-category", type=int, default=120)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(exist_ok=True)

    print("Loading GO OBO file...")
    valid_go_ids = load_valid_go_ids(OBO_PATH)
    print(f"  {len(valid_go_ids):,} valid GO IDs")

    print("Building reference set R and per-protein ground truth...")
    R, protein_gt, acc_organism = build_reference_and_ground_truth(RESULTS_DIR)
    print(f"  |R| = {len(R):,} unique experimentally validated GO terms")

    print("Collecting unique out-of-benchmark / misattribution tuples...")
    oob, misattr = collect_unmatched_tuples(RESULTS_DIR, valid_go_ids, R, protein_gt)
    print(f"  out-of-benchmark: {len(oob):,} unique tuples")
    print(f"  misattribution:   {len(misattr):,} unique tuples")

    random.seed(args.seed)
    oob_sample = random.sample(sorted(oob.keys()), min(args.n_per_category, len(oob)))
    mis_sample = random.sample(sorted(misattr.keys()), min(args.n_per_category, len(misattr)))
    sample_rows = (
        [(a, n, g, "out_of_benchmark") for a, n, g in oob_sample]
        + [(a, n, g, "misattribution") for a, n, g in mis_sample]
    )
    accessions = {a for a, _, _, _ in sample_rows}
    print(f"Sampled {len(sample_rows)} predictions across {len(accessions)} unique accessions.")

    print("Querying current UniProt REST API (this may take a few minutes)...")
    uniprot_cache = fetch_current_uniprot_go(accessions)

    classified = classify(sample_rows, uniprot_cache)

    with open(OUTPUT_DIR / "temporal_holdout_sample.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["accession", "organism", "namespace", "go_id", "category", "current_status", "current_evidence"])
        for acc, ns, go_id, category, status, evidence in classified:
            w.writerow([acc, acc_organism.get(acc, ""), ns, go_id, category, status, evidence])

    summary_rows = []
    for category in ["out_of_benchmark", "misattribution"]:
        cat_rows = [r for r in classified if r[3] == category]
        total = len(cat_rows)
        counts = Counter(r[4] for r in cat_rows)
        validated = counts.get("present_experimental", 0) + counts.get("present_non_experimental", 0)
        lo, hi = clopper_pearson(validated, total)
        summary_rows.append({
            "category": category,
            "n_sampled": total,
            "n_present_experimental": counts.get("present_experimental", 0),
            "n_present_non_experimental": counts.get("present_non_experimental", 0),
            "n_absent_current": counts.get("absent_current", 0),
            "n_protein_not_found": counts.get("protein_not_found", 0),
            "pct_validated": round(100 * validated / total, 2) if total else 0.0,
            "ci95_lower_pct": round(100 * lo, 2),
            "ci95_upper_pct": round(100 * hi, 2),
        })

    all_rows = classified
    total_all = len(all_rows)
    counts_all = Counter(r[4] for r in all_rows)
    validated_all = counts_all.get("present_experimental", 0) + counts_all.get("present_non_experimental", 0)
    lo_all, hi_all = clopper_pearson(validated_all, total_all)
    summary_rows.append({
        "category": "combined",
        "n_sampled": total_all,
        "n_present_experimental": counts_all.get("present_experimental", 0),
        "n_present_non_experimental": counts_all.get("present_non_experimental", 0),
        "n_absent_current": counts_all.get("absent_current", 0),
        "n_protein_not_found": counts_all.get("protein_not_found", 0),
        "pct_validated": round(100 * validated_all / total_all, 2) if total_all else 0.0,
        "ci95_lower_pct": round(100 * lo_all, 2),
        "ci95_upper_pct": round(100 * hi_all, 2),
    })

    with open(OUTPUT_DIR / "temporal_holdout_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        w.writeheader()
        w.writerows(summary_rows)

    print("\n=== Temporal-Holdout Validation Summary ===")
    for row in summary_rows:
        print(
            f"  {row['category']:>18s}: {row['n_present_experimental'] + row['n_present_non_experimental']}"
            f"/{row['n_sampled']} validated ({row['pct_validated']}%, "
            f"95% CI [{row['ci95_lower_pct']}%, {row['ci95_upper_pct']}%])"
        )
    print(f"\nSaved: {OUTPUT_DIR / 'temporal_holdout_sample.csv'}")
    print(f"Saved: {OUTPUT_DIR / 'temporal_holdout_summary.csv'}")


if __name__ == "__main__":
    main()
