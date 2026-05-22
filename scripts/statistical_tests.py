"""
PROBE: Statistical Testing Script
Cross-Organism Evaluation of General-Purpose LLMs for GO-Based
Protein Function Prediction

Tests performed:
  1. Wilcoxon Signed-Rank Test — pairwise model comparison
  2. Friedman Test — prompt format significance
  3. Kruskal-Wallis Test — organism significance
  4. Kruskal-Wallis Test — namespace significance
  5. Effect size (Cohen's d) for significant pairs

Outputs:
  - statistical_results.xlsx  (all tables)
  - significance_matrix.csv   (model pair p-values)
  - Console summary report

Requirements:
  pip install pandas scipy numpy openpyxl

Usage:
  python statistical_tests.py
  python statistical_tests.py --eval_dir probe_evaluation --output_dir probe_statistics
"""

import argparse
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# ─────────────────────────────────────────────
# CONFIG
# ─────────────────────────────────────────────

EVAL_DIR   = Path("probe_evaluation")
OUTPUT_DIR = Path("probe_statistics")
ALPHA      = 0.05   # significance threshold

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

ORGANISM_ORDER = ["human", "mouse", "zebrafish", "yeast", "ecoli"]

PROMPT_ORDER = [
    "P1_zeroshot",
    "P2_constrained",
    "P3_fewshot",
    "P4_cot",
    "P5_selection",
]

# ─────────────────────────────────────────────
# LOAD DATA
# ─────────────────────────────────────────────

def load_raw(eval_dir: Path) -> pd.DataFrame:
    path = eval_dir / "evaluation_raw.csv"
    if not path.exists():
        print(f"❌ {path} not found. Run evaluation.py first.")
        exit(1)
    df = pd.read_csv(path)
    df["model_label"] = df["model"].map(MODEL_LABELS).fillna(df["model"])
    print(f"Loaded {len(df):,} evaluation records.")
    print(f"Models   : {df['model'].nunique()}")
    print(f"Organisms: {df['organism'].nunique()}")
    print(f"Prompts  : {df['prompt_id'].nunique()}")
    return df


# ─────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────

def cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    """Cohen's d effect size between two arrays."""
    n1, n2   = len(a), len(b)
    mean_diff = np.mean(a) - np.mean(b)
    pooled_std = np.sqrt(
        ((n1 - 1) * np.std(a, ddof=1)**2 + (n2 - 1) * np.std(b, ddof=1)**2)
        / (n1 + n2 - 2)
    )
    return mean_diff / pooled_std if pooled_std > 0 else 0.0


def effect_label(d: float) -> str:
    d = abs(d)
    if d < 0.2:   return "negligible"
    elif d < 0.5: return "small"
    elif d < 0.8: return "medium"
    else:          return "large"


def stars(p: float) -> str:
    if p < 0.001: return "***"
    elif p < 0.01: return "**"
    elif p < 0.05: return "*"
    else:          return "ns"


def bonferroni(p: float, n_tests: int) -> float:
    """Bonferroni-corrected p-value."""
    return min(p * n_tests, 1.0)


# ─────────────────────────────────────────────
# TEST 1 — Pairwise Model Comparison (Wilcoxon)
# ─────────────────────────────────────────────

def test1_pairwise_models(df: pd.DataFrame) -> tuple:
    """
    Wilcoxon signed-rank test for every pair of models.
    Uses per-protein F1 scores (averaged across prompts and namespaces).
    Returns: (pvalue_matrix_df, detail_df)
    """
    print("\n[Test 1] Pairwise Model Comparison (Wilcoxon Signed-Rank) ...")

    # Average F1 per protein per model
    protein_f1 = (df.groupby(["model_label", "accession"])["f1"]
                    .mean()
                    .reset_index())

    models = sorted(protein_f1["model_label"].unique())
    n_models = len(models)
    n_pairs  = n_models * (n_models - 1) // 2

    # Build pivot: rows=accession, cols=model
    pivot = protein_f1.pivot(index="accession", columns="model_label", values="f1").fillna(0)

    rows = []
    p_matrix = pd.DataFrame(np.nan, index=models, columns=models)
    d_matrix = pd.DataFrame(np.nan, index=models, columns=models)

    for m1, m2 in itertools.combinations(models, 2):
        if m1 not in pivot.columns or m2 not in pivot.columns:
            continue
        a = pivot[m1].values
        b = pivot[m2].values

        # Wilcoxon needs matched pairs — use shared accessions
        shared = pivot[[m1, m2]].dropna()
        a = shared[m1].values
        b = shared[m2].values

        if len(a) < 10:
            continue

        stat, p_raw = stats.wilcoxon(a, b, alternative="two-sided")
        p_corrected  = bonferroni(p_raw, n_pairs)
        d            = cohens_d(a, b)

        p_matrix.loc[m1, m2] = round(p_corrected, 6)
        p_matrix.loc[m2, m1] = round(p_corrected, 6)
        d_matrix.loc[m1, m2] = round(d, 4)
        d_matrix.loc[m2, m1] = round(-d, 4)

        rows.append({
            "model_A":          m1,
            "model_B":          m2,
            "mean_f1_A":        round(np.mean(a), 4),
            "mean_f1_B":        round(np.mean(b), 4),
            "wilcoxon_stat":    round(stat, 4),
            "p_raw":            round(p_raw, 6),
            "p_corrected":      round(p_corrected, 6),
            "significant":      p_corrected < ALPHA,
            "cohens_d":         round(d, 4),
            "effect_size":      effect_label(d),
            "better_model":     m1 if np.mean(a) > np.mean(b) else m2,
            "sig_stars":        stars(p_corrected),
        })

    detail_df = pd.DataFrame(rows)
    if not detail_df.empty and "p_corrected" in detail_df.columns:
        detail_df = detail_df.sort_values("p_corrected")
    n_sig = detail_df["significant"].sum()
    print(f"    Model pairs tested  : {len(rows)}")
    print(f"    Significant (p<.05) : {n_sig} / {len(rows)}")

    return p_matrix, d_matrix, detail_df


# ─────────────────────────────────────────────
# TEST 2 — Prompt Format (Friedman)
# ─────────────────────────────────────────────

def test2_prompt_format(df: pd.DataFrame) -> pd.DataFrame:
    """
    Friedman test: do prompt formats produce significantly different F1?
    Run separately for each model.
    """
    print("\n[Test 2] Prompt Format Significance (Friedman Test) ...")

    rows = []
    for model in df["model_label"].unique():
        mdf = df[df["model_label"] == model]

        # F1 per protein per prompt (averaged across namespaces)
        pivot = (mdf.groupby(["accession", "prompt_id"])["f1"]
                    .mean()
                    .reset_index()
                    .pivot(index="accession", columns="prompt_id", values="f1")
                    .dropna())

        prompts = [p for p in PROMPT_ORDER if p in pivot.columns]
        if len(prompts) < 2 or len(pivot) < 5:
            continue

        groups = [pivot[p].values for p in prompts]
        stat, p = stats.friedmanchisquare(*groups)

        # Mean F1 per prompt for this model
        means = {p: round(pivot[p].mean(), 4) for p in prompts}
        best_prompt = max(means, key=means.get)

        rows.append({
            "model":           model,
            "friedman_stat":   round(stat, 4),
            "p_value":         round(p, 6),
            "significant":     p < ALPHA,
            "sig_stars":       stars(p),
            "best_prompt":     best_prompt,
            "best_prompt_f1":  means[best_prompt],
            **{f"f1_{k}": v for k, v in means.items()},
        })

    result_df = pd.DataFrame(rows)
    if not result_df.empty and "p_value" in result_df.columns:
        result_df = result_df.sort_values("p_value")
    n_sig = result_df["significant"].sum() if not result_df.empty else 0
    print(f"    Models tested       : {len(rows)}")
    print(f"    Significant (p<.05) : {n_sig} / {len(rows)}")
    return result_df


# ─────────────────────────────────────────────
# TEST 3 — Organism Effect (Kruskal-Wallis)
# ─────────────────────────────────────────────

def test3_organism_effect(df: pd.DataFrame) -> pd.DataFrame:
    """
    Kruskal-Wallis test: does organism significantly affect F1?
    Run separately for each model.
    """
    print("\n[Test 3] Organism Effect (Kruskal-Wallis Test) ...")

    rows = []
    for model in df["model_label"].unique():
        mdf = df[df["model_label"] == model]

        groups = []
        org_means = {}
        for org in ORGANISM_ORDER:
            vals = mdf[mdf["organism"] == org]["f1"].values
            if len(vals) > 0:
                groups.append(vals)
                org_means[org] = round(vals.mean(), 4)

        if len(groups) < 2:
            continue

        stat, p = stats.kruskal(*groups)

        rows.append({
            "model":         model,
            "kruskal_stat":  round(stat, 4),
            "p_value":       round(p, 6),
            "significant":   p < ALPHA,
            "sig_stars":     stars(p),
            "best_organism": max(org_means, key=org_means.get),
            **{f"f1_{k}": v for k, v in org_means.items()},
        })

    result_df = pd.DataFrame(rows)
    if not result_df.empty and "p_value" in result_df.columns:
        result_df = result_df.sort_values("p_value")
    n_sig = result_df["significant"].sum() if not result_df.empty else 0
    print(f"    Models tested       : {len(rows)}")
    print(f"    Significant (p<.05) : {n_sig} / {len(rows)}")
    return result_df


# ─────────────────────────────────────────────
# TEST 4 — Namespace Effect (Kruskal-Wallis)
# ─────────────────────────────────────────────

def test4_namespace_effect(df: pd.DataFrame) -> pd.DataFrame:
    """
    Kruskal-Wallis test: does GO namespace (MF/BP/CC) affect F1?
    Run separately for each model.
    """
    print("\n[Test 4] Namespace Effect (Kruskal-Wallis Test) ...")

    rows = []
    namespaces = ["molecular_function", "biological_process", "cellular_component"]

    for model in df["model_label"].unique():
        mdf = df[df["model_label"] == model]

        groups = []
        ns_means = {}
        for ns in namespaces:
            vals = mdf[mdf["namespace"] == ns]["f1"].values
            if len(vals) > 0:
                groups.append(vals)
                ns_means[ns[:2].upper()] = round(vals.mean(), 4)

        if len(groups) < 2:
            continue

        stat, p = stats.kruskal(*groups)

        rows.append({
            "model":        model,
            "kruskal_stat": round(stat, 4),
            "p_value":      round(p, 6),
            "significant":  p < ALPHA,
            "sig_stars":    stars(p),
            "best_namespace": max(ns_means, key=ns_means.get),
            **{f"f1_{k}": v for k, v in ns_means.items()},
        })

    result_df = pd.DataFrame(rows)
    if not result_df.empty and "p_value" in result_df.columns:
        result_df = result_df.sort_values("p_value")
    n_sig = result_df["significant"].sum() if not result_df.empty else 0
    print(f"    Models tested       : {len(rows)}")
    print(f"    Significant (p<.05) : {n_sig} / {len(rows)}")
    return result_df


# ─────────────────────────────────────────────
# EXPORT
# ─────────────────────────────────────────────

def export(p_matrix, d_matrix, t1, t2, t3, t4, output_dir: Path):
    output_dir.mkdir(exist_ok=True)
    print(f"\n[Export] Saving to {output_dir}/ ...")

    # CSVs
    p_matrix.to_csv(output_dir / "significance_matrix_pvalues.csv")
    d_matrix.to_csv(output_dir / "effect_size_matrix_cohens_d.csv")
    t1.to_csv(output_dir / "test1_pairwise_models.csv", index=False)
    t2.to_csv(output_dir / "test2_prompt_format.csv",   index=False)
    t3.to_csv(output_dir / "test3_organism_effect.csv", index=False)
    t4.to_csv(output_dir / "test4_namespace_effect.csv",index=False)

    # Single Excel workbook
    excel_path = output_dir / "PROBE_statistical_tests.xlsx"
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        p_matrix.to_excel(writer, sheet_name="pvalue_matrix")
        d_matrix.to_excel(writer, sheet_name="cohens_d_matrix")
        t1.to_excel(writer, sheet_name="test1_pairwise",  index=False)
        t2.to_excel(writer, sheet_name="test2_prompt",    index=False)
        t3.to_excel(writer, sheet_name="test3_organism",  index=False)
        t4.to_excel(writer, sheet_name="test4_namespace", index=False)

    print(f"    Saved: {excel_path}")
    for f in output_dir.glob("*.csv"):
        print(f"    Saved: {f}")


# ─────────────────────────────────────────────
# PRINT SUMMARY
# ─────────────────────────────────────────────

def print_summary(p_matrix, t1, t2, t3, t4):
    print("\n" + "=" * 65)
    print("PROBE — Statistical Testing Summary")
    print("=" * 65)

    # Test 1
    print("\n📊 Test 1: Pairwise Model Comparisons (Wilcoxon, Bonferroni-corrected)")
    sig = t1[t1["significant"]]
    print(f"   Significant pairs: {len(sig)} / {len(t1)}")
    if not sig.empty:
        print("   Top significant pairs:")
        for _, row in sig.head(5).iterrows():
            print(f"   {row['better_model']} > {row['model_A'] if row['better_model']==row['model_B'] else row['model_A']}  "
                  f"p={row['p_corrected']:.4f}{row['sig_stars']}  d={row['cohens_d']:.3f} ({row['effect_size']})")

    # Test 2
    print("\n📊 Test 2: Prompt Format Effect (Friedman Test per model)")
    sig2 = t2[t2["significant"]]
    print(f"   Models where prompt format matters: {len(sig2)} / {len(t2)}")
    if not t2.empty and "best_prompt" in t2.columns:
        best_counts = t2["best_prompt"].value_counts()
        print(f"   Most often best prompt: {best_counts.index[0]} ({best_counts.iloc[0]} models)")

    # Test 3
    print("\n📊 Test 3: Organism Effect (Kruskal-Wallis per model)")
    if not t3.empty and "significant" in t3.columns:
        sig3 = t3[t3["significant"]]
        print(f"   Models where organism matters: {len(sig3)} / {len(t3)}")
        best_org = t3["best_organism"].value_counts()
        print(f"   Most often best organism: {best_org.index[0]}")
    else:
        print("   Organism test skipped (only 1 organism in test data)")

    # Test 4
    print("\n📊 Test 4: Namespace Effect (Kruskal-Wallis per model)")
    sig4 = t4[t4["significant"]]
    print(f"   Models where namespace matters: {len(sig4)} / {len(t4)}")
    if not t4.empty and "best_namespace" in t4.columns:
        best_ns = t4["best_namespace"].value_counts()
        print(f"   Most often best namespace: {best_ns.index[0]}")

    print("\n" + "=" * 65)
    print("✅ Statistical testing complete!")
    print(f"   Results in: {OUTPUT_DIR}/")
    print(f"   Excel file: {OUTPUT_DIR}/PROBE_statistical_tests.xlsx")
    print("\nKey for significance stars:")
    print("   *** p < 0.001   ** p < 0.01   * p < 0.05   ns = not significant")


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PROBE Statistical Testing")
    parser.add_argument("--eval_dir",   default="probe_evaluation")
    parser.add_argument("--output_dir", default="probe_statistics")
    args = parser.parse_args()

    EVAL_DIR   = Path(args.eval_dir)
    OUTPUT_DIR = Path(args.output_dir)

    print("=" * 65)
    print("PROBE — Statistical Testing Pipeline")
    print("=" * 65)

    df = load_raw(EVAL_DIR)

    p_matrix, d_matrix, t1 = test1_pairwise_models(df)
    t2 = test2_prompt_format(df)
    t3 = test3_organism_effect(df)
    t4 = test4_namespace_effect(df)

    export(p_matrix, d_matrix, t1, t2, t3, t4, OUTPUT_DIR)
    print_summary(p_matrix, t1, t2, t3, t4)
