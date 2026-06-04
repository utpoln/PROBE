"""
PROBE: Confidence Intervals — fixed model matching
Run from PROBE root: python scripts/confidence_intervals.py
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats

OUTPUT_DIR = Path("probe_statistics")
OUTPUT_DIR.mkdir(exist_ok=True)

N_BOOTSTRAP = 10000
SEED        = 42
rng         = np.random.default_rng(SEED)

print("Loading evaluation_raw.csv...")
raw = pd.read_csv("probe_evaluation/evaluation_raw.csv")
print(f"Loaded {len(raw):,} records")

# Show exact model names in the file
print("\nExact model names in data:")
for m in sorted(raw['model'].unique()):
    print(f"  '{m}'")

def bootstrap_ci(values, n_boot=N_BOOTSTRAP):
    values = np.array(values)
    boot_means = [
        np.mean(rng.choice(values, size=len(values), replace=True))
        for _ in range(n_boot)
    ]
    mean  = np.mean(values)
    lower = np.percentile(boot_means, 2.5)
    upper = np.percentile(boot_means, 97.5)
    return mean, lower, upper

# ── Use exact model names from data ───────────────────────────
rows = []
print(f"\n{'Model':25s} {'N':>7} {'F1':>8} {'95% CI':>22} {'HR':>8} {'95% CI HR':>22}")
print("-"*100)

for model_name in sorted(raw['model'].unique()):
    sub = raw[raw['model'] == model_name]
    n   = len(sub)

    f1_vals = sub['f1'].dropna().values
    hr_vals = sub['hallucination_rate'].dropna().values

    f1_mean, f1_lo, f1_hi = bootstrap_ci(f1_vals)
    hr_mean, hr_lo, hr_hi = bootstrap_ci(hr_vals)

    rows.append({
        "model":        model_name,
        "n":            n,
        "f1_mean":      round(f1_mean, 4),
        "f1_ci_lower":  round(f1_lo,   4),
        "f1_ci_upper":  round(f1_hi,   4),
        "hr_mean":      round(hr_mean,  4),
        "hr_ci_lower":  round(hr_lo,   4),
        "hr_ci_upper":  round(hr_hi,   4),
    })

    print(f"  {model_name:23s} {n:>7,} {f1_mean:>8.4f} "
          f"[{f1_lo:.4f}–{f1_hi:.4f}] "
          f"{hr_mean:>8.4f} [{hr_lo:.4f}–{hr_hi:.4f}]")

df = pd.DataFrame(rows).sort_values("f1_mean", ascending=False)
df.to_csv(OUTPUT_DIR / "confidence_intervals.csv", index=False)
print(f"\nSaved: {OUTPUT_DIR}/confidence_intervals.csv")

# ── Top-2 significance test ────────────────────────────────────
print("\n── TOP-2 MODEL COMPARISON ──")
models_sorted = df['model'].tolist()
if len(models_sorted) >= 2:
    m1_name = models_sorted[0]
    m2_name = models_sorted[1]
    m1_f1   = raw[raw['model']==m1_name]['f1'].dropna().values
    m2_f1   = raw[raw['model']==m2_name]['f1'].dropna().values
    n_min   = min(len(m1_f1), len(m2_f1))
    stat, p = stats.wilcoxon(m1_f1[:n_min], m2_f1[:n_min])
    gap     = df.iloc[0]['f1_mean'] - df.iloc[1]['f1_mean']
    print(f"  {m1_name} vs {m2_name}")
    print(f"  F1 gap = {gap:.4f}")
    print(f"  Wilcoxon p = {p:.6f}")
    if p < 0.05:
        print("  → Gap IS statistically significant ✅")
    else:
        print("  → Gap is NOT significant ⚠️  — soften claim in paper")

# ── Save summary ───────────────────────────────────────────────
with open(OUTPUT_DIR / "ci_summary.txt", "w") as f:
    f.write("PROBE CONFIDENCE INTERVALS — 95% Bootstrap (10,000 resamples)\n")
    f.write("="*65 + "\n\n")
    for _, r in df.iterrows():
        f.write(f"{r['model']}:\n")
        f.write(f"  F1 = {r['f1_mean']:.4f} "
                f"(95% CI: {r['f1_ci_lower']:.4f}–{r['f1_ci_upper']:.4f})\n")
        f.write(f"  HR = {r['hr_mean']:.4f} "
                f"(95% CI: {r['hr_ci_lower']:.4f}–{r['hr_ci_upper']:.4f})\n\n")

print(f"Saved: {OUTPUT_DIR}/ci_summary.txt")