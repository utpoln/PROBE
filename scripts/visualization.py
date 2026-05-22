"""
PROBE: Visualization Script
Generates all paper-ready figures for BMC Bioinformatics submission

Figures produced:
  Fig 1 — Overall F1 leaderboard (bar chart, all models)
  Fig 2 — Model × Namespace heatmap (MF, BP, CC)
  Fig 3 — Cross-organism F1 curve (organism degradation)
  Fig 4 — Prompt sensitivity heatmap (Model × Prompt)
  Fig 5 — Hallucination rate bar chart
  Fig 6 — Precision-Recall scatter plot per model
  Fig 7 — F1 distribution boxplot per model

Requirements:
  pip install matplotlib seaborn pandas numpy

Usage:
  python visualization.py
  python visualization.py --eval_dir probe_evaluation --output_dir probe_figures
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

# ─────────────────────────────────────────────
# CONFIG
# ─────────────────────────────────────────────

EVAL_DIR   = Path("probe_evaluation")
OUTPUT_DIR = Path("probe_figures")

# Clean model name mapping for display
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

# Organism display order (human → most distant)
ORGANISM_ORDER = ["yeast", "human", "ecoli", "mouse", "zebrafish"]
ORGANISM_LABELS = {
    "human":     "H. sapiens",
    "mouse":     "M. musculus",
    "zebrafish": "D. rerio",
    "yeast":     "S. cerevisiae",
    "ecoli":     "E. coli",
}

PROMPT_LABELS = {
    "P1_zeroshot":    "P1: Zero-shot",
    "P2_constrained": "P2: Constrained",
    "P3_fewshot":     "P3: Few-shot",
    "P4_cot":         "P4: Chain-of-Thought",
    "P5_selection":   "P5: Selection",
}

# Paper-quality style
PALETTE    = sns.color_palette("colorblind", 10)
FIG_DPI    = 300
FIG_FORMAT = "pdf"   # PDF for submission; change to "png" for previews

plt.rcParams.update({
    "font.family":      "DejaVu Sans",
    "font.size":        11,
    "axes.titlesize":   13,
    "axes.labelsize":   12,
    "xtick.labelsize":  10,
    "ytick.labelsize":  10,
    "legend.fontsize":  10,
    "figure.dpi":       FIG_DPI,
    "axes.spines.top":  False,
    "axes.spines.right":False,
    "axes.grid":        True,
    "grid.alpha":       0.3,
    "grid.linestyle":   "--",
})


# ─────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────

def load_tables(eval_dir: Path) -> dict:
    tables = {}
    for csv_file in eval_dir.glob("*.csv"):
        tables[csv_file.stem] = pd.read_csv(csv_file)
    if not tables:
        print(f"❌ No CSV files found in {eval_dir}/")
        print("   Run evaluation.py first.")
        exit(1)
    print(f"Loaded {len(tables)} evaluation tables.")
    return tables


def clean_model_name(name: str) -> str:
    return MODEL_LABELS.get(name, name)


def save(fig, name: str, output_dir: Path):
    output_dir.mkdir(exist_ok=True)
    path_pdf = output_dir / f"{name}.pdf"
    path_png = output_dir / f"{name}.png"
    fig.savefig(path_pdf, bbox_inches="tight", dpi=FIG_DPI)
    fig.savefig(path_png, bbox_inches="tight", dpi=FIG_DPI)
    print(f"    Saved: {path_png}  +  {path_pdf}")
    plt.close(fig)


# ─────────────────────────────────────────────
# FIGURE 1 — Overall Leaderboard
# ─────────────────────────────────────────────

def fig1_leaderboard(tables: dict, output_dir: Path):
    print("\n[Fig 1] Overall F1 Leaderboard ...")
    df = tables["leaderboard"].copy()
    df["model_label"] = df["model"].apply(clean_model_name)
    df = df.sort_values("f1", ascending=True)

    fig, ax = plt.subplots(figsize=(9, 5))
    colors = [PALETTE[i % len(PALETTE)] for i in range(len(df))]
    bars = ax.barh(df["model_label"], df["f1"], color=colors,
                   edgecolor="white", linewidth=0.5, height=0.6)

    # Add value labels
    for bar, val in zip(bars, df["f1"]):
        ax.text(bar.get_width() + 0.002, bar.get_y() + bar.get_height() / 2,
                f"{val:.4f}", va="center", ha="left", fontsize=9)

    ax.set_xlabel("Macro-F1 Score")
    # title removed — caption goes in LaTeX
    # ax.set_title("Fig. 1 — Overall LLM Performance on GO Function Prediction (PROBE Benchmark)",      pad=12)
    ax.set_xlim(0, df["f1"].max() * 1.18)
    fig.tight_layout()
    save(fig, "fig1_leaderboard", output_dir)


# ─────────────────────────────────────────────
# FIGURE 2 — Model × Namespace Heatmap
# ─────────────────────────────────────────────

def fig2_namespace_heatmap(tables: dict, output_dir: Path):
    print("\n[Fig 2] Model × Namespace Heatmap ...")
    df = tables["model_x_namespace"].copy()
    df["model_label"] = df["model"].apply(clean_model_name)

    pivot = df.pivot(index="model_label", columns="namespace_short", values="f1")
    # Reorder columns
    cols = [c for c in ["MF", "BP", "CC"] if c in pivot.columns]
    pivot = pivot[cols]
    pivot = pivot.sort_values("MF", ascending=False)

    fig, ax = plt.subplots(figsize=(7, 6))
    sns.heatmap(
        pivot, annot=True, fmt=".4f", cmap="YlOrRd",
        linewidths=0.5, linecolor="white",
        cbar_kws={"label": "F1 Score"},
        ax=ax, vmin=0, vmax=pivot.values.max()
    )
    # title removed — caption goes in LaTeX
    ax.set_xlabel("GO Namespace")
    ax.set_ylabel("Model")
    ax.tick_params(axis="x", rotation=0)
    ax.tick_params(axis="y", rotation=0)
    fig.tight_layout()
    save(fig, "fig2_namespace_heatmap", output_dir)


# ─────────────────────────────────────────────
# FIGURE 3 — Cross-Organism Curve
# ─────────────────────────────────────────────

def fig3_organism_curve(tables: dict, output_dir: Path):
    print("\n[Fig 3] Cross-Organism F1 Curve ...")
    df = tables["model_x_organism"].copy()
    df["model_label"] = df["model"].apply(clean_model_name)

    # Filter to organisms in order
    df = df[df["organism"].isin(ORGANISM_ORDER)]
    df["organism"] = pd.Categorical(df["organism"],
                                     categories=ORGANISM_ORDER, ordered=True)
    df = df.sort_values("organism")

    fig, ax = plt.subplots(figsize=(9, 5))
    models = df["model_label"].unique()

    for i, model in enumerate(models):
        subset = df[df["model_label"] == model].sort_values("organism")
        x = [ORGANISM_LABELS.get(o, o) for o in subset["organism"]]
        y = subset["f1"].values
        ax.plot(x, y, marker="o", label=model,
                color=PALETTE[i % len(PALETTE)], linewidth=1.8,
                markersize=6)

    ax.set_xlabel("Organism (ordered by average F1, reflecting research literature density)")
    ax.set_ylabel("Macro-F1 Score")
    # title removed — caption goes in LaTeX
    ax.legend(loc="upper right", fontsize=8, ncol=2,
              framealpha=0.8, edgecolor="lightgrey")
    fig.tight_layout()
    save(fig, "fig3_organism_curve", output_dir)


# ─────────────────────────────────────────────
# FIGURE 4 — Prompt Sensitivity Heatmap
# ─────────────────────────────────────────────

def fig4_prompt_heatmap(tables: dict, output_dir: Path):
    print("\n[Fig 4] Prompt Sensitivity Heatmap ...")
    df = tables["model_x_prompt"].copy()
    df["model_label"] = df["model"].apply(clean_model_name)
    df["prompt_label"] = df["prompt_id"].map(PROMPT_LABELS).fillna(df["prompt_id"])

    pivot = df.pivot(index="model_label", columns="prompt_label", values="f1")
    # Sort columns by prompt order
    ordered_cols = [PROMPT_LABELS[p] for p in
                    ["P1_zeroshot", "P2_constrained", "P3_fewshot",
                     "P4_cot", "P5_selection"]
                    if PROMPT_LABELS[p] in pivot.columns]
    pivot = pivot[ordered_cols]
    pivot = pivot.sort_values(ordered_cols[0], ascending=False)

    fig, ax = plt.subplots(figsize=(10, 6))
    sns.heatmap(
        pivot, annot=True, fmt=".4f", cmap="Blues",
        linewidths=0.5, linecolor="white",
        cbar_kws={"label": "F1 Score"},
        ax=ax, vmin=0
    )
    # title removed — caption goes in LaTeX
    ax.set_xlabel("Prompt Format")
    ax.set_ylabel("Model")
    ax.tick_params(axis="x", rotation=15)
    ax.tick_params(axis="y", rotation=0)
    fig.tight_layout()
    save(fig, "fig4_prompt_heatmap", output_dir)


# ─────────────────────────────────────────────
# FIGURE 5 — Hallucination Rate
# ─────────────────────────────────────────────

def fig5_hallucination(tables: dict, output_dir: Path):
    print("\n[Fig 5] Hallucination Rate ...")
    df = tables["leaderboard"].copy()
    df["model_label"] = df["model"].apply(clean_model_name)
    df = df.sort_values("hallucination_rate", ascending=True)
    df["halluc_pct"] = df["hallucination_rate"] * 100

    fig, ax = plt.subplots(figsize=(9, 5))
    colors = ["#d62728" if v > 50 else "#ff7f0e" if v > 25 else "#2ca02c"
              for v in df["halluc_pct"]]
    bars = ax.barh(df["model_label"], df["halluc_pct"],
                   color=colors, edgecolor="white", height=0.6)

    for bar, val in zip(bars, df["halluc_pct"]):
        ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height() / 2,
                f"{val:.1f}%", va="center", fontsize=9)

    ax.axvline(50, color="red", linestyle="--", alpha=0.4, linewidth=1)
    ax.set_xlabel("Hallucination Rate (% predicted GO IDs not in ground truth)")
    # title removed — caption goes in LaTeX
    ax.set_xlim(0, 110)

    legend = [
        mpatches.Patch(color="#2ca02c", label="Low (<25%)"),
        mpatches.Patch(color="#ff7f0e", label="Medium (25–50%)"),
        mpatches.Patch(color="#d62728", label="High (>50%)"),
    ]
    ax.legend(handles=legend, loc="lower right", fontsize=9)
    fig.tight_layout()
    save(fig, "fig5_hallucination", output_dir)


# ─────────────────────────────────────────────
# FIGURE 6 — Precision-Recall Scatter
# ─────────────────────────────────────────────

def fig6_precision_recall(tables: dict, output_dir: Path):
    print("\n[Fig 6] Precision-Recall Scatter ...")
    df = tables["leaderboard"].copy()
    df["model_label"] = df["model"].apply(clean_model_name)

    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot points
    for i, (_, row) in enumerate(df.iterrows()):
        ax.scatter(row["recall"], row["precision"],
                   color=PALETTE[i % len(PALETTE)],
                   s=160, zorder=3, edgecolors="white", linewidth=1.2)

    # Smart label placement — alternate above/below to avoid overlap
    offsets = [
        (8, 6), (-75, 8), (8, -14), (-80, -14),
        (8, 6), (-85, 8), (8, -14), (8, 6),
        (-85, -14), (8, 6),
    ]
    for i, (_, row) in enumerate(df.iterrows()):
        ox, oy = offsets[i % len(offsets)]
        ax.annotate(row["model_label"],
                    (row["recall"], row["precision"]),
                    textcoords="offset points", xytext=(ox, oy),
                    fontsize=8, color="black",
                    arrowprops=dict(arrowstyle="-", color="grey",
                                   lw=0.5) if abs(ox) > 10 else None)

    # F1 iso-curves — only relevant range for our data
    for f1_val in [0.05, 0.08, 0.10, 0.12]:
        p = np.linspace(0.01, 0.3, 300)
        r = f1_val * p / (2 * p - f1_val)
        mask = (r > 0) & (r <= 0.3)
        if mask.sum() > 5:
            ax.plot(r[mask], p[mask], "--", color="grey",
                    alpha=0.3, linewidth=0.9)
            mid = len(r[mask]) // 3
            ax.annotate(f"F1={f1_val}",
                        (r[mask][mid], p[mask][mid]),
                        fontsize=7.5, color="grey", alpha=0.7)

    # Tight axes around actual data
    pad = 0.015
    ax.set_xlim(df["recall"].min() - pad, df["recall"].max() + pad + 0.06)
    ax.set_ylim(df["precision"].min() - pad, df["precision"].max() + pad + 0.02)

    ax.set_xlabel("Recall", fontsize=11)
    ax.set_ylabel("Precision", fontsize=11)
    ax.grid(True, alpha=0.3, linestyle="--")
    # No title — caption in LaTeX
    fig.tight_layout()
    save(fig, "fig6_precision_recall", output_dir)


# ─────────────────────────────────────────────
# FIGURE 7 — F1 Boxplot per Model
# ─────────────────────────────────────────────

def fig7_f1_boxplot(tables: dict, output_dir: Path):
    print("\n[Fig 7] F1 Distribution Boxplot ...")
    df = tables["evaluation_raw"].copy()
    df["model_label"] = df["model"].apply(clean_model_name)

    # Sort by mean F1 from leaderboard (per-protein median is 0 for most
    # models due to 68% complete-miss rate, so mean is more informative)
    order = (tables["leaderboard"]
             .assign(model_label=lambda x: x["model"].apply(clean_model_name))
             .sort_values("f1", ascending=False)["model_label"]
             .tolist())

    # Keep only models present in both
    present = df["model_label"].unique()
    order = [m for m in order if m in present]

    fig, ax = plt.subplots(figsize=(12, 5))
    sns.boxplot(data=df, x="model_label", y="f1",
                order=order, palette="colorblind",
                linewidth=0.8, fliersize=1.5, ax=ax)

    ax.set_xlabel("Model (ordered by mean F1, best on left)", fontsize=10)
    ax.set_ylabel("F1 Score (per protein)", fontsize=10)
    ax.set_ylim(-0.02, 1.02)
    ax.tick_params(axis="x", rotation=30)
    ax.grid(True, axis="y", alpha=0.3, linestyle="--")
    # No title — caption in LaTeX
    fig.tight_layout()
    save(fig, "fig7_f1_boxplot", output_dir)


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PROBE Visualization Script")
    parser.add_argument("--eval_dir",   default="probe_evaluation",
                        help="Directory with evaluation CSV files")
    parser.add_argument("--output_dir", default="probe_figures",
                        help="Directory to save figures")
    args = parser.parse_args()

    EVAL_DIR   = Path(args.eval_dir)
    OUTPUT_DIR = Path(args.output_dir)

    print("=" * 55)
    print("PROBE — Visualization Pipeline")
    print("=" * 55)

    tables = load_tables(EVAL_DIR)

    fig1_leaderboard(tables, OUTPUT_DIR)
    fig2_namespace_heatmap(tables, OUTPUT_DIR)
    fig3_organism_curve(tables, OUTPUT_DIR)
    fig4_prompt_heatmap(tables, OUTPUT_DIR)
    fig5_hallucination(tables, OUTPUT_DIR)
    fig6_precision_recall(tables, OUTPUT_DIR)
    fig7_f1_boxplot(tables, OUTPUT_DIR)

    print(f"\n✅ All 7 figures saved to: {OUTPUT_DIR}/")
    print("   Each figure saved as both .pdf (submission) and .png (preview)")
