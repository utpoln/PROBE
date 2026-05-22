# PROBE: A Cross-Organism Benchmark for Evaluating General-Purpose LLMs on Zero-Shot Gene Ontology Protein Function Prediction

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

PROBE (**PR**otein functi**O**n **B**enchmark **E**valuation) is a systematic benchmark for evaluating whether general-purpose open-source Large Language Models (LLMs) can predict Gene Ontology (GO) protein function terms in a zero-shot, inference-only setting.

**Key findings:**
- Best model (Mistral Large 123B) achieves macro-F1 of **9.3%** — far below specialized tools
- Hallucination rates range from **30.6% to 50.4%** (mean 41.3%)
- *S. cerevisiae* (yeast) outperforms *H. sapiens* — reflecting **research literature density**, not evolutionary distance
- Candidate-selection prompting (P5) reduces hallucination by **33.5 percentage points**
- **42/45** model pairs significantly different (Wilcoxon, Bonferroni-corrected)

---

## Repository Structure

```
PROBE/
├── benchmark_dataset/
│   └── benchmark_master.json      # 1,000 proteins, 5 organisms, experimental GO only
├── scripts/
│   ├── dataset_preparation.py     # Swiss-Prot parsing, filtering, stratified sampling
│   ├── llm_prompting.py           # MindRouter API calls, 5 prompt formats, auto-resume
│   ├── evaluation.py              # F1, Precision, Recall, Hallucination Rate
│   ├── statistical_tests.py       # Wilcoxon, Friedman, Kruskal-Wallis, Cohen's d
│   ├── error_analysis.py          # Error categorization (5 categories)
│   └── visualization.py           # All 7 publication figures
├── probe_evaluation/
│   ├── leaderboard.csv            # Overall model rankings
│   ├── model_x_namespace.csv      # F1 by model and GO namespace
│   ├── model_x_organism.csv       # F1 by model and organism
│   ├── model_x_prompt.csv         # F1 by model and prompt format
│   ├── organism_x_namespace.csv   # F1 by organism and namespace
│   ├── full_detail.csv            # Full breakdown
│   └── PROBE_results.xlsx         # All tables in one Excel workbook
├── probe_figures/
│   ├── fig0_pipeline.png          # Pipeline overview diagram
│   ├── fig1_leaderboard.pdf       # Macro-F1 leaderboard
│   ├── fig2_namespace_heatmap.pdf # F1 by model and namespace
│   ├── fig3_organism_curve.pdf    # Cross-organism performance
│   ├── fig4_prompt_heatmap.pdf    # Prompt sensitivity
│   ├── fig5_hallucination.pdf     # Hallucination rates
│   ├── fig6_precision_recall.pdf  # Precision-Recall scatter
│   └── fig7_f1_boxplot.pdf        # Per-protein F1 distributions
├── paper/
│   ├── probe_paper_final.tex      # LaTeX manuscript
│   └── references.bib             # Bibliography
├── requirements.txt
└── README.md
```

---

## Dataset

**Source:** UniProt Swiss-Prot release 2024_03

**Construction:**
- 5 organisms: *H. sapiens*, *M. musculus*, *D. rerio*, *S. cerevisiae*, *E. coli*
- Evidence codes: IDA, IMP, IPI, IGI, IEP, EXP and high-throughput variants only
- IEA (electronic annotation) excluded to prevent circular reasoning
- Minimum 3 experimental GO terms per protein
- 200 proteins per organism = **1,000 total** (stratified random sample, seed=42)

**Statistics:** Mean 7.7 GO terms/protein, max 64, min 3. MF: 68%, BP: 94%, CC: 76%

---

## Models Evaluated

| Model | Family | Size | Architecture |
|---|---|---|---|
| Mistral Large 123B | Mistral | 123B | Dense |
| Llama 3.3 70B | Llama 3 | 70B | Dense |
| Qwen2.5 72B | Qwen 2.5 | 72B | Dense |
| Qwen3 32B | Qwen 3 | 32B | Dense |
| Gemma-4 31B | Gemma | 31B | Dense |
| Llama 3.1 8B | Llama 3 | 8B | Dense |
| Mixtral 8×7B | Mistral | 56B (14B active) | MoE |
| Qwen2.5 7B | Qwen 2.5 | 7B | Dense |
| Mistral 7B | Mistral | 7B | Dense |
| Gemma3 12B | Gemma | 12B | Dense |

All models accessed via [University of Idaho MindRouter API](https://mindrouter.uidaho.edu/v1) at T=0, max 512 tokens.

---

## Prompt Formats

| ID | Name | Description |
|---|---|---|
| P1 | Zero-shot | Direct GO term prediction |
| P2 | Constrained | Namespace-restricted zero-shot |
| P3 | Few-shot | 3 annotated examples (TP53, EGFR, ACT1) |
| P4 | Chain-of-thought | Step-by-step biological reasoning |
| P5 | Candidate selection | Select from 15 pre-validated GO terms |

---

## Reproduce Results

### 1. Install dependencies

```bash
pip install -r requirements.txt
```

### 2. Build dataset

Download Swiss-Prot flat file from https://www.uniprot.org/downloads and place as `uniprot_sprot.dat` in the root folder, then:

```bash
python scripts/dataset_preparation.py
```

Output: `benchmark_dataset/benchmark_master.json`

### 3. Run LLM inference

Requires a [MindRouter API key](https://mindrouter.uidaho.edu). Run all 10 models:

```bash
python scripts/llm_prompting.py --api_key YOUR_KEY --all_models
```

Or run a single model:

```bash
python scripts/llm_prompting.py --api_key YOUR_KEY --model mistral-large:123b
```

The script auto-resumes if interrupted. Results saved to `probe_results/`.

**Note:** 150,000 API calls total. Estimated runtime: 5–20 hours per model depending on size.

### 4. Evaluate

```bash
python scripts/evaluation.py
```

Output: `probe_evaluation/` — all CSV tables and Excel workbook.

### 5. Statistical tests

```bash
python scripts/statistical_tests.py
```

Output: `probe_statistics/` — Wilcoxon, Friedman, Kruskal-Wallis results.

### 6. Error analysis

```bash
python scripts/error_analysis.py
```

Output: `probe_error_analysis/` — error categorization and sampled failure cases.

### 7. Generate figures

```bash
python scripts/visualization.py
```

Output: `probe_figures/` — all 7 publication-ready figures (PDF + PNG).

---

## Results Summary

### Leaderboard

| Rank | Model | F1 | Precision | Recall | HR (%) |
|---|---|---|---|---|---|
| 1 | Mistral Large 123B | 0.093 | 0.089 | 0.140 | 37.1 |
| 2 | Llama 3.3 70B | 0.092 | 0.082 | 0.144 | 38.2 |
| 3 | Qwen2.5 72B | 0.086 | 0.078 | 0.133 | 30.6 |
| 4 | Qwen3 32B | 0.076 | 0.074 | 0.110 | 33.9 |
| 5 | Gemma-4 31B | 0.069 | 0.064 | 0.097 | 41.7 |
| 6 | Llama 3.1 8B | 0.062 | 0.058 | 0.084 | 39.3 |
| 7 | Mixtral 8×7B | 0.059 | 0.049 | 0.109 | 45.7 |
| 8 | Qwen2.5 7B | 0.044 | 0.049 | 0.050 | 47.4 |
| 9 | Mistral 7B | 0.042 | 0.031 | 0.087 | 49.0 |
| 10 | Gemma3 12B | 0.032 | 0.028 | 0.053 | 50.4 |

HR = Hallucination Rate (predicted GO IDs absent from full reference set)

### Key Findings

**Namespace:** CC (0.141) >> MF (0.031) > BP (0.024) — all 10/10 models significant (KW, p<0.001)

**Organism:** Yeast (0.086) > Human (0.081) > E. coli (0.069) > Mouse (0.065) > Zebrafish (0.026) — all 10/10 models significant (KW, p<0.001)

**Prompt:** P5 (0.082) > P3 (0.063) > P4 (0.062) > P1 (0.060) > P2 (0.060) — 10/10 models significant (Friedman, p<0.001); P5 best for 7/10 models

**Error analysis:** Complete miss 68.71% | Complete hallucination 8.65% | Correct 0.34%

---

## Citation

If you use PROBE in your research, please cite:

```bibtex
@article{naha2025probe,
  title={PROBE: A Cross-Organism Benchmark for Evaluating General-Purpose
         Large Language Models on Zero-Shot Gene Ontology Protein Function Prediction},
  author={Naha, Kallol and Jamil, Hasan M.},
  journal={BMC Bioinformatics},
  year={2025},
  institution={University of Idaho}
}
```

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

## Contact

- Kallol Naha — naha7197@vandals.uidaho.edu
- Hasan M. Jamil — jamil@uidaho.edu
- Department of Computer Science, University of Idaho
