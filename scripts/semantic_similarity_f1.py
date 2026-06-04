"""
PROBE: Semantic Similarity F1 — no pronto needed
Parses go-basic.obo directly, no external dependencies beyond
what is already in requirements.txt.

Run from PROBE root:
    python scripts/semantic_similarity_f1.py

Output:
    probe_statistics/semantic_f1_leaderboard.csv
"""

import json
import re
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict

# ── Find paths ─────────────────────────────────────────────────
_script = Path(__file__).resolve()
_root   = _script.parent.parent if (_script.parent.parent / "probe_results").exists() \
          else Path(".")

OBO_PATH     = _root / "go-basic.obo"
RESULTS_DIR  = _root / "probe_results"
OUTPUT_DIR   = _root / "probe_statistics"
OUTPUT_DIR.mkdir(exist_ok=True)

print(f"Root:        {_root}")
print(f"OBO file:    {OBO_PATH}  exists={OBO_PATH.exists()}")
print(f"Results dir: {RESULTS_DIR}  exists={RESULTS_DIR.exists()}")

# ── Parse go-basic.obo ─────────────────────────────────────────
def parse_obo(obo_path):
    """
    Parse go-basic.obo and return:
      parents: dict {go_id: set of parent go_ids}
    """
    parents = defaultdict(set)
    current = None

    print(f"\nParsing {obo_path} ...")
    with open(obo_path) as f:
        for line in f:
            line = line.strip()
            if line == "[Term]":
                current = None
            elif line.startswith("id: GO:"):
                current = line.split("id: ")[1].strip()
            elif line.startswith("is_a:") and current:
                # is_a: GO:0000001 ! ...
                parent = line.split("is_a:")[1].strip().split()[0]
                if parent.startswith("GO:"):
                    parents[current].add(parent)
            elif line.startswith("relationship: part_of") and current:
                parent = line.split("part_of")[1].strip().split()[0]
                if parent.startswith("GO:"):
                    parents[current].add(parent)

    print(f"Parsed {len(parents):,} GO terms with parent relationships")
    return parents

# ── Build ancestor cache ───────────────────────────────────────
def build_ancestors(parents):
    """Compute full ancestor set for each GO term via BFS."""
    print("Building ancestor cache...")
    ancestors = {}

    def get_anc(term, visited=None):
        if term in ancestors:
            return ancestors[term]
        if visited is None:
            visited = set()
        if term in visited:
            return set()
        visited.add(term)
        anc = set()
        for p in parents.get(term, set()):
            anc.add(p)
            anc.update(get_anc(p, visited))
        ancestors[term] = anc
        return anc

    all_terms = set(parents.keys())
    for t in all_terms:
        get_anc(t)

    print(f"Cached ancestors for {len(ancestors):,} terms")
    return ancestors

# ── Semantic similarity (Jaccard of ancestor sets) ─────────────
def jaccard_sim(go_a, go_b, ancestors):
    if go_a == go_b:
        return 1.0
    anc_a = ancestors.get(go_a, set()) | {go_a}
    anc_b = ancestors.get(go_b, set()) | {go_b}
    inter = len(anc_a & anc_b)
    union = len(anc_a | anc_b)
    return inter / union if union else 0.0

def semantic_prec(pred, gt, ancestors):
    if not pred:
        return 0.0
    return np.mean([max((jaccard_sim(p, g, ancestors) for g in gt), default=0.0)
                    for p in pred])

def semantic_rec(pred, gt, ancestors):
    if not gt:
        return 0.0
    return np.mean([max((jaccard_sim(g, p, ancestors) for p in pred), default=0.0)
                    for g in gt])

def sem_f1(pred, gt, ancestors):
    if not pred and not gt:
        return 1.0
    if not pred or not gt:
        return 0.0
    p = semantic_prec(pred, gt, ancestors)
    r = semantic_rec(pred, gt, ancestors)
    return 2*p*r/(p+r) if p+r else 0.0

def exact_f1(pred, gt):
    if not pred and not gt:
        return 1.0
    if not pred or not gt:
        return 0.0
    tp = len(pred & gt)
    p = tp/len(pred) if pred else 0
    r = tp/len(gt)   if gt   else 0
    return 2*p*r/(p+r) if p+r else 0.0

def parse_go_ids(text):
    if not text:
        return set()
    return set(re.findall(r'GO:\d{7}', str(text)))

def parse_gt(gt_str):
    if not gt_str:
        return set()
    ids = set()
    for item in str(gt_str).split('|'):
        item = item.strip()
        if item.startswith('GO:'):
            ids.add(item.split()[0])
    return ids

# ── Model labels ───────────────────────────────────────────────
MODEL_LABELS = {
    "llama3.1-8b":        "Llama 3.1 8B",
    "llama3.3-70b":       "Llama 3.3 70B",
    "mistral-7b":         "Mistral 7B",
    "mistral-large-123b": "Mistral Large 123B",
    "mixtral-8x7b":       "Mixtral 8x7B",
    "qwen2.5-7b":         "Qwen2.5 7B",
    "qwen2.5-72b":        "Qwen2.5 72B",
    "Qwen_Qwen3-32B":     "Qwen3 32B",
    "gemma3-12b":         "Gemma3 12B",
    "google_gemma-4-31b": "Gemma-4 31B",
}

# ── Main ───────────────────────────────────────────────────────
if not OBO_PATH.exists():
    print(f"\nERROR: {OBO_PATH} not found")
    print("Download with: curl -L -o go-basic.obo https://purl.obolibrary.org/obo/go/go-basic.obo")
    exit(1)

parents   = parse_obo(OBO_PATH)
ancestors = build_ancestors(parents)

jsonl_files = list(RESULTS_DIR.glob("*.jsonl"))
print(f"\nFound {len(jsonl_files)} result files")

rows = []
for jfile in sorted(jsonl_files):
    label = MODEL_LABELS.get(jfile.stem, jfile.stem)
    print(f"Processing: {jfile.name} → {label}")

    exact_scores, sem_scores = [], []
    n_processed = 0

    with open(jfile) as f:
        for line in f:
            try:
                rec = json.loads(line)
                if not rec.get("success", False):
                    continue
                pred = parse_go_ids(rec.get("response", ""))
                gt   = parse_gt(rec.get("ground_truth", ""))
                if not gt:
                    continue

                ef = exact_f1(pred, gt)
                sf = sem_f1(pred, gt, ancestors)

                exact_scores.append(ef)
                sem_scores.append(sf)
                n_processed += 1

            except Exception:
                continue

    if exact_scores:
        em = np.mean(exact_scores)
        sm = np.mean(sem_scores)
        rows.append({
            "model":        label,
            "n":            n_processed,
            "exact_f1":     round(em, 4),
            "semantic_f1":  round(sm, 4),
            "gain":         round(sm - em, 4),
            "gain_pct":     round((sm - em) / em * 100, 1) if em > 0 else 0,
        })
        print(f"  n={n_processed:,}  exact={em:.4f}  semantic={sm:.4f}  gain={sm-em:+.4f}")

# ── Results ────────────────────────────────────────────────────
df = pd.DataFrame(rows).sort_values("exact_f1", ascending=False)

print("\n" + "="*70)
print("SEMANTIC vs EXACT F1 COMPARISON")
print("="*70)
print(f"\n{'Model':25s} {'Exact F1':>10} {'Semantic F1':>12} {'Gain':>8} {'Gain%':>8}")
print("-"*65)
for _, r in df.iterrows():
    print(f"  {r['model']:23s} {r['exact_f1']:>10.4f} "
          f"{r['semantic_f1']:>12.4f} {r['gain']:>+8.4f} {r['gain_pct']:>+7.1f}%")

df.to_csv(OUTPUT_DIR / "semantic_f1_leaderboard.csv", index=False)
print(f"\nSaved: {OUTPUT_DIR}/semantic_f1_leaderboard.csv")

print("\n── FOR PAPER ──")
mean_gain = df['gain'].mean()
mean_gain_pct = df['gain_pct'].mean()
print(f"Mean semantic gain: {mean_gain:+.4f} ({mean_gain_pct:+.1f}%)")
print(f"Best model exact F1:    {df.iloc[0]['exact_f1']:.4f}")
print(f"Best model semantic F1: {df.iloc[0]['semantic_f1']:.4f}")