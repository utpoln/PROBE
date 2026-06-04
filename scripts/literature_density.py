"""
PROBE: Literature Density Analysis
Fetches PubMed paper counts per organism and correlates
with per-organism F1 scores.

Usage:
    python literature_density.py

Output:
    probe_statistics/literature_density.csv
    probe_figures/fig_literature_correlation.pdf
"""

import requests
import time
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
from scipy import stats
from pathlib import Path

OUTPUT_DIR = Path("probe_statistics")
FIG_DIR    = Path("probe_figures")
OUTPUT_DIR.mkdir(exist_ok=True)
FIG_DIR.mkdir(exist_ok=True)

# ── Per-organism F1 from PROBE evaluation ──────────────────────
organism_f1 = {
    "Saccharomyces cerevisiae": 0.086,
    "Homo sapiens":             0.081,
    "Escherichia coli":         0.069,
    "Mus musculus":             0.065,
    "Danio rerio":              0.026,
}

# Pool sizes from benchmark
pool_sizes = {
    "Saccharomyces cerevisiae": 4560,
    "Homo sapiens":             11042,
    "Escherichia coli":         2177,
    "Mus musculus":             8758,
    "Danio rerio":              759,
}

# ── PubMed search queries per organism ─────────────────────────
pubmed_queries = {
    "Homo sapiens":             "Homo sapiens[Organism] AND protein function[Title/Abstract]",
    "Mus musculus":             "Mus musculus[Organism] AND protein function[Title/Abstract]",
    "Danio rerio":              "Danio rerio[Organism] AND protein function[Title/Abstract]",
    "Saccharomyces cerevisiae": "Saccharomyces cerevisiae[Organism] AND protein function[Title/Abstract]",
    "Escherichia coli":         "Escherichia coli[Organism] AND protein function[Title/Abstract]",
}

def fetch_pubmed_count(query):
    """Fetch paper count from PubMed E-utilities API."""
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db":      "pubmed",
        "term":    query,
        "rettype": "count",
        "retmode": "json",
    }
    try:
        r = requests.get(url, params=params, timeout=30)
        r.raise_for_status()
        count = int(r.json()["esearchresult"]["count"])
        return count
    except Exception as e:
        print(f"  Error fetching '{query}': {e}")
        return None

# ── Fetch counts ───────────────────────────────────────────────
print("Fetching PubMed paper counts per organism...")
print("(This may take ~30 seconds due to API rate limits)\n")

pubmed_counts = {}
for organism, query in pubmed_queries.items():
    print(f"  Querying: {organism}...")
    count = fetch_pubmed_count(query)
    pubmed_counts[organism] = count
    print(f"    → {count:,} papers")
    time.sleep(0.4)  # NCBI rate limit: max 3 requests/sec

# ── Build dataframe ────────────────────────────────────────────
rows = []
for org in organism_f1:
    rows.append({
        "organism":      org,
        "f1":            organism_f1[org],
        "pubmed_count":  pubmed_counts.get(org),
        "pool_size":     pool_sizes[org],
        "pubmed_per_protein": pubmed_counts.get(org, 0) / pool_sizes[org]
            if pubmed_counts.get(org) else None,
    })

df = pd.DataFrame(rows)
print("\n── Results ──")
print(df[["organism","f1","pubmed_count","pool_size","pubmed_per_protein"]].to_string(index=False))

# ── Correlations ───────────────────────────────────────────────
print("\n── Correlations ──")

# 1. F1 vs PubMed count
r1, p1 = stats.spearmanr(df["f1"], df["pubmed_count"])
print(f"Spearman r (F1 vs PubMed count):           r={r1:.3f}, p={p1:.4f}")

# 2. F1 vs PubMed per protein
r2, p2 = stats.spearmanr(df["f1"], df["pubmed_per_protein"])
print(f"Spearman r (F1 vs PubMed per protein):     r={r2:.3f}, p={p2:.4f}")

# 3. F1 vs pool size (to test if pool size drives it)
r3, p3 = stats.spearmanr(df["f1"], df["pool_size"])
print(f"Spearman r (F1 vs pool size):              r={r3:.3f}, p={p3:.4f}")

# ── Save CSV ───────────────────────────────────────────────────
df["spearman_r_f1_pubmed"]           = r1
df["spearman_p_f1_pubmed"]           = p1
df["spearman_r_f1_pubmed_per_prot"]  = r2
df["spearman_p_f1_pubmed_per_prot"]  = p2
df["spearman_r_f1_pool_size"]        = r3
df["spearman_p_f1_pool_size"]        = p3
df.to_csv(OUTPUT_DIR / "literature_density.csv", index=False)
print(f"\nSaved: {OUTPUT_DIR}/literature_density.csv")

# ── Figure ─────────────────────────────────────────────────────
short_names = {
    "Saccharomyces cerevisiae": "S. cerevisiae",
    "Homo sapiens":             "H. sapiens",
    "Escherichia coli":         "E. coli",
    "Mus musculus":             "M. musculus",
    "Danio rerio":              "D. rerio",
}

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Plot 1: F1 vs PubMed count
ax = axes[0]
for _, row in df.iterrows():
    ax.scatter(row["pubmed_count"], row["f1"], s=120, zorder=5)
    ax.annotate(short_names[row["organism"]],
                (row["pubmed_count"], row["f1"]),
                textcoords="offset points", xytext=(6, 4), fontsize=9)
ax.set_xlabel("PubMed Papers (protein function)", fontsize=11)
ax.set_ylabel("Average F1 Score", fontsize=11)
ax.set_title(f"F1 vs Literature Density\nSpearman r={r1:.3f}, p={p1:.3f}", fontsize=11)
ax.grid(True, alpha=0.3)

# Add trend line
x = df["pubmed_count"].values
y = df["f1"].values
m, b = np.polyfit(x, y, 1)
xline = np.linspace(x.min(), x.max(), 100)
ax.plot(xline, m*xline+b, 'r--', alpha=0.5, linewidth=1.5)

# Plot 2: F1 vs PubMed per protein
ax = axes[1]
for _, row in df.iterrows():
    ax.scatter(row["pubmed_per_protein"], row["f1"], s=120, zorder=5)
    ax.annotate(short_names[row["organism"]],
                (row["pubmed_per_protein"], row["f1"]),
                textcoords="offset points", xytext=(6, 4), fontsize=9)
ax.set_xlabel("PubMed Papers per Annotated Protein", fontsize=11)
ax.set_ylabel("Average F1 Score", fontsize=11)
ax.set_title(f"F1 vs Literature Density (per protein)\nSpearman r={r2:.3f}, p={p2:.3f}", fontsize=11)
ax.grid(True, alpha=0.3)

x2 = df["pubmed_per_protein"].values
m2, b2 = np.polyfit(x2, y, 1)
xline2 = np.linspace(x2.min(), x2.max(), 100)
ax.plot(xline2, m2*xline2+b2, 'r--', alpha=0.5, linewidth=1.5)

plt.tight_layout()
plt.savefig(FIG_DIR / "fig_literature_correlation.pdf", bbox_inches='tight', dpi=300)
plt.savefig(FIG_DIR / "fig_literature_correlation.png", bbox_inches='tight', dpi=300)
plt.close()
print(f"Saved: {FIG_DIR}/fig_literature_correlation.pdf")

print("\n── FOR PAPER ──")
print(f"Spearman correlation between per-organism F1 and PubMed paper count:")
print(f"  r = {r1:.3f}, p = {p1:.4f}")
print(f"Spearman correlation between per-organism F1 and PubMed papers per protein:")
print(f"  r = {r2:.3f}, p = {p2:.4f}")
print(f"Spearman correlation between per-organism F1 and pool size:")
print(f"  r = {r3:.3f}, p = {p3:.4f}")
print("\nIf r(F1, pubmed) > r(F1, pool_size) → literature density drives the effect")
print("If r(F1, pool_size) > r(F1, pubmed) → pool size drives the effect")
