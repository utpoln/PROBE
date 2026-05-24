import json
import re
from pathlib import Path
from collections import defaultdict
import csv

# Step 1: Parse GO OBO file
print("Loading GO OBO file...")
valid_go_ids = set()
with open("go-basic.obo") as f:
    for line in f:
        line = line.strip()
        if line.startswith("id: GO:"):
            go_id = line.replace("id: ", "").strip()
            valid_go_ids.add(go_id)
print(f"Total valid GO IDs in ontology: {len(valid_go_ids):,}")

# Step 2: Build reference set R from all ground truth
print("Building reference set R from ground truth...")
results_dir = Path("probe_results")
all_files = list(results_dir.glob("*.jsonl"))

R = set()
protein_gt = defaultdict(set)

for filepath in all_files:
    with open(filepath) as f:
        for line in f:
            try:
                rec = json.loads(line)
                gt = rec.get("ground_truth", "")
                if gt:
                    terms = re.findall(r"GO:\d{7}", gt)
                    accession = rec["accession"]
                    namespace = rec["namespace"]
                    key = (accession, namespace)
                    protein_gt[key].update(terms)
                    R.update(terms)
            except:
                continue

print(f"Reference set R size: {len(R):,} unique GO terms")
print(f"Proteins with ground truth: {len(protein_gt):,}")

# Step 3: Decompose hallucination for each model
print("Decomposing hallucination rates...")
results = []

for filepath in sorted(all_files):
    model_name = filepath.stem
    total_predicted = 0
    syntactic = 0
    out_of_benchmark = 0
    misattribution = 0
    correct = 0

    with open(filepath) as f:
        for line in f:
            try:
                rec = json.loads(line)
                if not rec.get("success", False):
                    continue
                response = rec.get("response", "")
                accession = rec["accession"]
                namespace = rec["namespace"]
                key = (accession, namespace)
                gt_terms = protein_gt.get(key, set())
                if not gt_terms:
                    continue
                predicted = set(re.findall(r"GO:\d{7}", response))
                if not predicted:
                    continue
                for go_id in predicted:
                    total_predicted += 1
                    if go_id in gt_terms:
                        correct += 1
                    elif go_id not in valid_go_ids:
                        syntactic += 1
                    elif go_id not in R:
                        out_of_benchmark += 1
                    else:
                        misattribution += 1
            except:
                continue

    if total_predicted > 0:
        results.append({
            "model": model_name,
            "total_predicted": total_predicted,
            "correct_pct": round(100 * correct / total_predicted, 1),
            "syntactic_pct": round(100 * syntactic / total_predicted, 1),
            "out_of_benchmark_pct": round(100 * out_of_benchmark / total_predicted, 1),
            "misattribution_pct": round(100 * misattribution / total_predicted, 1),
        })

# Step 4: Print results
print("\n" + "="*85)
print(f"{'Model':<30} {'Correct%':>9} {'(a)Syntax%':>11} {'(b)OutBench%':>13} {'(c)Misattr%':>12}")
print("="*85)
for r in sorted(results, key=lambda x: x["model"]):
    print(f"{r['model']:<30} {r['correct_pct']:>8}% {r['syntactic_pct']:>10}% {r['out_of_benchmark_pct']:>12}% {r['misattribution_pct']:>11}%")
print("="*85)

# Save CSV
with open("hallucination_decomposed.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=results[0].keys())
    writer.writeheader()
    writer.writerows(results)
print("Results saved to hallucination_decomposed.csv")
