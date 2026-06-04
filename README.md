# PROBE: A Cross-Organism Benchmark for Evaluating General-Purpose Large Language Models on Gene Ontology Protein Function Prediction

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![BMC Bioinformatics](https://img.shields.io/badge/journal-BMC%20Bioinformatics-green.svg)](https://bmcbioinformatics.biomedcentral.com/)

> **Kallol Naha, Hasan M. Jamil**
> Department of Computer Science, University of Idaho
> naha7197@vandals.uidaho.edu | jamil@uidaho.edu

---

## Overview

**PROBE** (**PR**otein functi**O**n **B**enchmark **E**valuation) is the first systematic benchmark for evaluating whether general-purpose Large Language Models (LLMs) can predict Gene Ontology (GO) protein function terms in a **zero-shot, inference-only setting** — no fine-tuning, no retrieval augmentation.

We evaluate **10 open-source LLMs** across **5 prompt formats**, **5 organisms**, and **3 GO namespaces**, yielding **150,000 API calls** on a curated benchmark of **1,000 proteins** from UniProt Swiss-Prot with experimental evidence-only annotations.

A focused supplementary comparison of **4 leading closed-source models** (Claude Sonnet 4.6, Gemini 2.5 Pro, GPT-5.4, Gemini 2.5 Flash) confirms that the annotation deficit is a fundamental LLM limitation regardless of model source or scale — the gap between best closed-source and best open-source is only **0.002 F1 points**.

GO semantic similarity validation (Wang/BMA via GOSemSim) confirms that LLM predictions are biologically meaningful beyond chance, exceeding random baselines by **2.1–2.5×** across all namespaces — showing that low exact F1 scores reflect a precision deficit rather than random guessing.

---

## Seven Key Findings

| # | Finding | Result |
|---|---|---|
| 1 | Best open-source model (Mistral Large 123B) | Macro-F1 = **9.3%** — far below specialized tools |
| 2 | Most-frequent-5 heuristic beats 4/10 LLMs | F1 = **0.062** outperforms Gemma3 12B, Mistral 7B, Qwen2.5 7B, Mixtral 8×7B |
| 3 | Hallucination decomposed | Syntactic **6.6%** · Out-of-benchmark **40.2%** · Misattribution **47.7%** |
| 4 | Best organism | *S. cerevisiae* (yeast) > *H. sapiens* — literature density effect |
| 5 | Best prompt format | P5 candidate selection — reduces hallucination by **33.5 pp** |
| 6 | Semantic similarity vs random | LLMs exceed random by **2.1–2.5×** across all namespaces |
| 7 | Closed-source gap | Best closed (Claude Sonnet 4.6, F1=0.127) vs best open (Mistral Large, F1=0.125) — gap of only **0.002** |

---

## Repository Structure

```
PROBE/
├── benchmark_dataset/
│   └── benchmark_master.json              # 1,000 proteins, 5 organisms, experimental GO only
│
├── scripts/
│   ├── dataset_preparation.py             # Swiss-Prot parsing, filtering, stratified sampling
│   ├── llm_prompting.py                   # API calls, 5 prompt formats, auto-resume
│   ├── evaluation.py                      # F1, Precision, Recall, Hallucination Rate
│   ├── statistical_tests.py               # Wilcoxon, Friedman, Kruskal-Wallis, Cohen's d
│   ├── error_analysis.py                  # 5-category error classification
│   ├── visualization.py                   # Publication-ready figures (fig1–fig7)
│   ├── hallucination_decompose.py         # Decompose HR into syntactic/out-of-bench/misattr
│   ├── closed_source_eval.py              # Closed-source model evaluation (P5, 50 proteins)
│   ├── ablation_function_text.py          # FUNCTION text ablation study (A1/A2/A3)
│   ├── confidence_intervals.py            # Bootstrap 95% CIs for all metrics
│   ├── literature_density.py              # PubMed organism correlation analysis
│   ├── semantic_similarity_f1.py          # Jaccard-based semantic F1 (no R required)
│   ├── semantic_similarity.R              # GOSemSim Wang/BMA analysis (top 3 models)
│   ├── semantic_similarity_all10.R        # GOSemSim Wang/BMA analysis (all 10 models)
│   └── random_baseline.R                  # Random baseline per namespace
│
├── probe_evaluation/
│   ├── leaderboard.csv                    # Overall model rankings
│   ├── model_x_namespace.csv              # F1 by model × GO namespace
│   ├── model_x_organism.csv               # F1 by model × organism
│   ├── model_x_prompt.csv                 # F1 by model × prompt format
│   ├── organism_x_namespace.csv           # F1 by organism × namespace
│   ├── full_detail.csv                    # Full per-record breakdown
│   ├── evaluation_raw.csv                 # Raw per-prediction scores (150,000 rows)
│   ├── hallucination_decomposed.csv       # Hallucination decomposition (a/b/c)
│   ├── semantic_similarity_results_all10.csv  # Wang/BMA scores all 10 models
│   └── PROBE_results.xlsx                 # All tables in one Excel workbook
│
├── probe_figures/
│   ├── fig0_pipeline.png                  # Pipeline overview diagram
│   ├── fig1_leaderboard.pdf/png           # Macro-F1 leaderboard
│   ├── fig2_namespace_heatmap.pdf/png     # F1 heatmap by model and namespace
│   ├── fig3_organism_curve.pdf/png        # Cross-organism F1 line plot
│   ├── fig4_prompt_heatmap.pdf/png        # Prompt sensitivity heatmap
│   ├── fig5_hallucination.pdf/png         # Aggregate hallucination rate
│   ├── fig6_precision_recall.pdf/png      # Precision-Recall scatter
│   ├── fig7_f1_boxplot.pdf/png            # Per-protein F1 distributions
│   ├── fig9_MF_semantic_similarity.pdf/png  # Semantic similarity — MF namespace
│   ├── fig9_BP_semantic_similarity.pdf/png  # Semantic similarity — BP namespace
│   └── fig9_CC_semantic_similarity.pdf/png  # Semantic similarity — CC namespace
│
├── probe_error_analysis/
│   ├── error_analysis_raw.csv             # Per-record error classification (150,000 rows)
│   ├── error_cases_sampled.csv            # 50 sampled failure cases per model
│   ├── error_counts_per_model.csv         # Error category counts per model
│   ├── error_pivot.csv                    # Error % pivot table
│   ├── hallucination_by_namespace.csv     # HR by model and namespace
│   ├── error_by_organism.csv              # Error counts by organism
│   ├── error_by_prompt.csv                # Error counts by prompt format
│   └── PROBE_error_analysis.xlsx          # All error tables in one workbook
│
├── probe_statistics/
│   ├── test1_pairwise_models.csv          # Wilcoxon results (45 pairs)
│   ├── test2_prompt_format.csv            # Friedman test results (10 models)
│   ├── test3_organism_effect.csv          # Kruskal-Wallis organism results
│   ├── test4_namespace_effect.csv         # Kruskal-Wallis namespace results
│   ├── significance_matrix_pvalues.csv    # 10×10 p-value matrix
│   ├── effect_size_matrix_cohens_d.csv    # 10×10 Cohen's d matrix
│   ├── confidence_intervals.csv           # 95% bootstrap CIs for all models
│   ├── ci_summary.txt                     # CI summary in plain text
│   ├── semantic_f1_leaderboard.csv        # Jaccard semantic F1 vs exact F1
│   ├── literature_density.csv             # PubMed counts and Spearman correlations
│   └── PROBE_statistical_tests.xlsx       # All statistical results in one workbook
│
├── probe_results/                         # 10 open-source model raw responses
│   ├── llama3.1-8b.jsonl                  # 15,000 raw LLM responses each
│   ├── llama3.3-70b.jsonl
│   ├── mistral-7b.jsonl
│   ├── mistral-large-123b.jsonl
│   ├── mixtral-8x7b.jsonl
│   ├── qwen2.5-7b.jsonl
│   ├── qwen2.5-72b.jsonl
│   ├── Qwen_Qwen3-32B.jsonl
│   ├── gemma3-12b.jsonl
│   └── google_gemma-4-31b.jsonl           # Total: 150,000 API call records
│
├── probe_results_ablation/                # FUNCTION text ablation results
│   ├── ablation_gemma3-12b.jsonl          # 900 records (100 proteins × 3 conditions × 3 ns)
│   ├── ablation_llama3-1-8b.jsonl
│   └── ablation_mistral-large-123b.jsonl
│
├── probe_results_closed/                  # Closed-source results (P5, 50 proteins)
│   ├── claude-sonnet-4-6.jsonl            # 150 records each
│   ├── gemini-2-5-pro.jsonl
│   ├── gemini-2-5-flash.jsonl
│   └── gpt-5-4-2026-03-05.jsonl
│
├── go-basic.obo                           # GO ontology file (release 2026-03-25, 48,291 terms)
├── requirements.txt
├── .gitignore
└── README.md
```

> **Note:** The manuscript and LaTeX source will be released upon journal acceptance.
> All 150,000 raw LLM response records are included in `probe_results/` for full reproducibility.

---

## Dataset

**Source:** UniProt Swiss-Prot release 2024\_03 (574,627 reviewed entries)

**Construction pipeline (`dataset_preparation.py`):**

1. **Organism filter** — 5 target organisms: *H. sapiens*, *M. musculus*, *D. rerio*, *S. cerevisiae*, *E. coli*
2. **Evidence filter** — IDA, IMP, IPI, IGI, IEP, EXP, HDA, HMP, HGI, HEP (IEA excluded to prevent circular reasoning)
3. **Depth filter** — minimum 3 experimental GO terms per protein
4. **Stratified sampling** — 200 proteins per organism (seed=42) = **1,000 total**

| Organism | N | Pool | Avg GO terms |
|---|---|---|---|
| *H. sapiens* | 200 | 11,042 | 7.7 |
| *M. musculus* | 200 | 8,758 | 7.7 |
| *D. rerio* | 200 | 759 | 7.7 |
| *S. cerevisiae* | 200 | 4,560 | 7.7 |
| *E. coli* | 200 | 2,177 | 7.7 |
| **Total** | **1,000** | **27,296** | **7.7** |

Overall: MF 68% · BP 94% · CC 76% · max 64 GO terms · min 3 (depth filter threshold)

---

## Models Evaluated

### Open-Source (primary benchmark — 150,000 API calls)

| Rank | Model | Family | Size | Architecture | Ollama command |
|---|---|---|---|---|---|
| 1 | Mistral Large 123B | Mistral | 123B | Dense | `ollama pull mistral-large` |
| 2 | Llama 3.3 70B | Llama 3 | 70B | Dense | `ollama pull llama3.3:70b` |
| 3 | Qwen2.5 72B | Qwen 2.5 | 72B | Dense | `ollama pull qwen2.5:72b` |
| 4 | Qwen3 32B | Qwen 3 | 32B | Dense | `ollama pull qwen3:32b` |
| 5 | Gemma-4 31B | Gemma | 31B | Dense | `ollama pull gemma4:31b` |
| 6 | Llama 3.1 8B | Llama 3 | 8B | Dense | `ollama pull llama3.1:8b` |
| 7 | Mixtral 8×7B | Mistral | 56B (14B active) | MoE | `ollama pull mixtral:8x7b` |
| 8 | Qwen2.5 7B | Qwen 2.5 | 7B | Dense | `ollama pull qwen2.5:7b` |
| 9 | Mistral 7B | Mistral | 7B | Dense | `ollama pull mistral:7b` |
| 10 | Gemma3 12B | Gemma | 12B | Dense | `ollama pull gemma3:12b` |

All models evaluated at T=0, max 512 tokens, no system prompt, inference only.

### Closed-Source (supplementary — P5 only, 50 proteins)

| Model | Provider | F1 | Prec | Recall |
|---|---|---|---|---|
| Claude Sonnet 4.6 | Anthropic | 0.127 | 0.132 | 0.150 |
| Gemini 2.5 Pro | Google | 0.118 | 0.123 | 0.141 |
| GPT-5.4 | OpenAI | 0.116 | 0.136 | 0.119 |
| Gemini 2.5 Flash | Google | 0.081 | 0.167 | 0.059 |
| **Mistral Large 123B (open)** | — | **0.125** | **0.131** | **0.152** |

Gap between best closed-source and best open-source: **0.002 F1 points**.

---

## Prompt Formats

| ID | Name | Description |
|---|---|---|
| P1 | Zero-shot | Free GO term prediction, no examples or constraints |
| P2 | Constrained | Namespace-restricted zero-shot with strict format rules |
| P3 | Few-shot (3-shot) | 3 examples: TP53, EGFR, ACT1 (verified absent from test set — no contamination) |
| P4 | Chain-of-thought | Step-by-step biological reasoning before prediction |
| P5 | Candidate selection | Select from 15 pre-validated GO terms per namespace |

> **Note on P5:** P5 changes the task from open-ended generation to constrained selection. P5 results should not be interpreted as autonomous GO discovery performance and are not directly comparable to open generation prompts (P1–P4).

Total: 1,000 × 5 × 3 × 10 = **150,000 API calls**

### Representative Prompt Examples (TP53, *H. sapiens*, Molecular Function)

**P1 — Zero-shot direct**
```
You are an expert molecular biologist.

Protein: TP53
Full name: Cellular tumor antigen p53
Organism: Homo sapiens (Human)
Known function: Acts as a tumor suppressor. Induces growth arrest
or apoptosis depending on physiological circumstances and cell type.

List the Gene Ontology (GO) Molecular Function terms that best
describe this protein. Format each as: GO:XXXXXXX | term name
```

**P2 — Constrained zero-shot**
```
[Same context block as P1]

Task: Predict ONLY the Gene Ontology Molecular Function (GO MF) terms.
Rules:
- List only experimentally supported functions
- Format each term as: GO:XXXXXXX | term name
- Do not include BP or CC terms
- Output only the GO terms, nothing else
```

**P3 — Few-shot (3-shot)**
```
Example 1:
Protein: TP53, Organism: Homo sapiens
GO:0003700 | DNA-binding transcription factor activity
GO:0046872 | metal ion binding
GO:0042802 | identical protein binding

Example 2:
Protein: EGFR, Organism: Homo sapiens
GO:0004714 | transmembrane receptor protein tyrosine kinase activity
GO:0005515 | protein binding
GO:0016301 | kinase activity

Example 3:
Protein: ACT1, Organism: Saccharomyces cerevisiae
GO:0005524 | ATP binding
GO:0003779 | actin binding
GO:0005198 | structural molecule activity

Now predict for: [same context block as P1]
List GO Molecular Function terms. Format: GO:XXXXXXX | term name
```

**P4 — Chain-of-thought**
```
[Same context block as P1]

Think step by step about what Molecular Function GO terms apply:
Step 1: What is the protein's primary biochemical activity?
Step 2: What molecular processes does it participate in?
Step 3: Based on your reasoning, list the GO Molecular Function terms.

Format final answer as: GO:XXXXXXX | term name
```

**P5 — Candidate selection**
```
[Same context block as P1]

From the following GO Molecular Function terms, select ALL that apply:
GO:0003700 | DNA-binding transcription factor activity
GO:0004672 | protein kinase activity
GO:0005515 | protein binding
GO:0005524 | ATP binding
GO:0003677 | DNA binding
GO:0046872 | metal ion binding
GO:0004722 | protein serine/threonine phosphatase activity
GO:0016301 | kinase activity
GO:0003723 | RNA binding
GO:0004725 | protein tyrosine phosphatase activity
GO:0003779 | actin binding
GO:0005198 | structural molecule activity
GO:0008270 | zinc ion binding
GO:0016787 | hydrolase activity
GO:0016740 | transferase activity

List only matching terms as: GO:XXXXXXX | term name
If none apply, write: NONE
```

---

## Reproduce Results

### Step 0 — Install dependencies

```bash
pip install -r requirements.txt
```

### Step 1 — Build dataset (optional — pre-built included)

```bash
python scripts/dataset_preparation.py
```

Output: `benchmark_dataset/benchmark_master.json`

### Step 2 — Run LLM Inference

#### Option A — University of Idaho MindRouter (internal access)

```bash
python scripts/llm_prompting.py --api_key YOUR_KEY --all_models
```

#### Option B — Ollama (public, local inference — no institutional access required)

```bash
curl -fsSL https://ollama.ai/install.sh | sh

ollama pull llama3.3:70b
ollama pull llama3.1:8b
ollama pull mistral:7b
ollama pull qwen2.5:7b
ollama pull qwen2.5:72b
ollama pull gemma3:12b
ollama pull mixtral:8x7b

python scripts/llm_prompting.py \
  --api_key ollama \
  --base_url http://localhost:11434/v1 \
  --model llama3.1:8b
```

> Hardware requirements: 7B models need 8GB VRAM · 70B models need 40GB+ VRAM

#### Option C — Any OpenAI-compatible API

```bash
python scripts/llm_prompting.py \
  --api_key YOUR_KEY \
  --base_url https://your-endpoint/v1 \
  --model your-model-name
```

The script auto-resumes if interrupted. Only successful calls are counted.

### Step 3 — Evaluate

```bash
python scripts/evaluation.py
```

Output: `probe_evaluation/` — all CSV tables and Excel workbook.

### Step 4 — Statistical Tests

```bash
python scripts/statistical_tests.py
```

Output: `probe_statistics/` — Wilcoxon, Friedman, Kruskal-Wallis, Cohen's d.

### Step 5 — Error Analysis

```bash
python scripts/error_analysis.py
```

### Step 6 — Hallucination Decomposition

```bash
python scripts/hallucination_decompose.py
```

Requires `go-basic.obo` (included). Decomposes hallucination rate into:
- **(a) Syntactic** — predicted GO ID absent from GO ontology
- **(b) Out-of-benchmark** — valid GO ID not in reference set R
- **(c) Misattribution** — valid GO ID in R but wrong protein

### Step 7 — Generate Figures

```bash
python scripts/visualization.py
```

### Step 8 — Semantic Similarity (Python, no R required)

```bash
python scripts/semantic_similarity_f1.py
```

Uses Jaccard similarity over GO ancestor sets (go-basic.obo included). Output: `probe_statistics/semantic_f1_leaderboard.csv`.

### Step 9 — Semantic Similarity (R/GOSemSim Wang/BMA)

```r
BiocManager::install(c("GOSemSim","org.Hs.eg.db","org.Sc.sgd.db",
                        "org.Mm.eg.db","org.Dr.eg.db","org.EcK12.eg.db",
                        "GO.db","jsonlite"))
source("scripts/semantic_similarity_all10.R")
```

Output: `probe_evaluation/semantic_similarity_results_all10.csv`.

### Step 10 — Confidence Intervals

```bash
python scripts/confidence_intervals.py
```

Output: `probe_statistics/confidence_intervals.csv` — 95% bootstrap CIs (10,000 resamples) for all models.

### Step 11 — Ablation Study

```bash
python scripts/ablation_function_text.py \
  --api_key YOUR_KEY \
  --base_url http://localhost:11434/v1 \
  --model mistral-large:123b \
  --limit 100
```

Tests whether performance reflects pre-training knowledge (A1: name only) vs FUNCTION text extraction (A3: full context).

### Step 12 — Closed-Source Evaluation

```bash
python scripts/closed_source_eval.py \
  --api_key YOUR_KEY \
  --model claude-sonnet-4-6 \
  --provider anthropic

python scripts/closed_source_eval.py \
  --api_key YOUR_KEY \
  --model gpt-5.4-2026-03-05 \
  --provider openai

python scripts/closed_source_eval.py \
  --api_key YOUR_KEY \
  --model gemini-2.5-pro \
  --provider google
```

---

## Results Summary

### Leaderboard

| Rank | Model | F1 | Prec | Recall | HR% | ERR% |
|---|---|---|---|---|---|---|
| 1 | Mistral Large 123B | 0.093 | 0.089 | 0.140 | 37.1 | 1.8 |
| 2 | Llama 3.3 70B | 0.092 | 0.082 | 0.144 | 38.2 | 0.9 |
| 3 | Qwen2.5 72B | 0.086 | 0.078 | 0.133 | 30.6 | 6.2 |
| 4 | Qwen3 32B | 0.076 | 0.074 | 0.110 | 33.9 | 10.1 |
| 5 | Gemma-4 31B | 0.069 | 0.064 | 0.097 | 41.7 | 2.8 |
| 6 | Llama 3.1 8B | 0.062 | 0.058 | 0.084 | 39.3 | 1.0 |
| — | **Most-frequent-5 heuristic** | **0.062** | — | — | **0.0** | — |
| 7 | Mixtral 8×7B | 0.059 | 0.049 | 0.109 | 45.7 | 2.8 |
| 8 | Qwen2.5 7B | 0.044 | 0.049 | 0.050 | 47.4 | 5.8 |
| 9 | Mistral 7B | 0.042 | 0.031 | 0.087 | 49.0 | 0.0 |
| 10 | Gemma3 12B | 0.032 | 0.028 | 0.053 | 50.4 | 0.7 |

HR = Hallucination Rate · ERR = Empty Response Rate

### Hallucination Decomposition

| Model | Correct % | (a) Syntactic % | (b) Out-of-bench % | (c) Misattr. % |
|---|---|---|---|---|
| Mistral Large 123B | 6.6 | 5.7 | 42.4 | 45.3 |
| Llama 3.3 70B | 6.7 | 4.5 | 42.9 | 45.9 |
| Qwen2.5 72B | 7.1 | 5.2 | 34.0 | 53.7 |
| Qwen3 32B | 7.2 | 5.3 | 30.2 | 57.3 |
| Gemma-4 31B | 6.5 | 8.5 | 37.1 | 47.9 |
| Llama 3.1 8B | 6.2 | 5.3 | 36.3 | 52.3 |
| Mixtral 8×7B | 4.2 | 7.4 | 47.9 | 40.6 |
| Qwen2.5 7B | 4.7 | 8.2 | 44.5 | 42.6 |
| Mistral 7B | 3.2 | 7.4 | 44.8 | 44.7 |
| Gemma3 12B | 2.7 | 8.4 | 42.1 | 46.8 |
| **Mean** | **5.7** | **6.6** | **40.2** | **47.7** |

The dominant failure mode is **misattribution (47.7%)** — models predict valid GO terms but assign them to the wrong protein. Syntactic hallucination (predicting non-existent GO IDs) accounts for only **6.6%**.

### Semantic Similarity (Wang/BMA, GOSemSim)

| Model | MF | BP | CC |
|---|---|---|---|
| Mistral Large 123B | 0.452 | 0.377 | 0.658 |
| Llama 3.3 70B | 0.478 | 0.373 | 0.654 |
| Qwen2.5 72B | 0.444 | 0.361 | 0.663 |
| Qwen3 32B | 0.377 | 0.330 | 0.618 |
| Gemma-4 31B | 0.345 | 0.281 | 0.609 |
| Llama 3.1 8B | 0.380 | 0.283 | 0.556 |
| Mixtral 8×7B | 0.383 | 0.300 | 0.558 |
| Mistral 7B | 0.321 | 0.259 | 0.504 |
| Qwen2.5 7B | 0.234 | 0.237 | 0.492 |
| Gemma3 12B | 0.281 | 0.196 | 0.390 |
| **Mean** | **0.370** | **0.299** | **0.571** |
| **Random baseline** | **0.145** | **0.133** | **0.258** |

All models exceed the random baseline in all namespaces (2.5× MF · 2.2× BP · 2.2× CC).

### By GO Namespace (avg across models)

| Namespace | Avg F1 | Avg HR |
|---|---|---|
| Cellular Component (CC) | 0.141 | 26.3% |
| Molecular Function (MF) | 0.031 | 47.8% |
| Biological Process (BP) | 0.024 | 49.9% |

### By Organism (avg across models)

| Organism | Avg F1 |
|---|---|
| *S. cerevisiae* (yeast) | 0.086 |
| *H. sapiens* (human) | 0.081 |
| *E. coli* | 0.069 |
| *M. musculus* (mouse) | 0.065 |
| *D. rerio* (zebrafish) | 0.026 |

### By Prompt Format (avg across models)

| Prompt | Avg F1 | Avg HR |
|---|---|---|
| P5 — Candidate selection | 0.082 | 16.3% |
| P3 — Few-shot | 0.063 | 43.9% |
| P4 — Chain-of-thought | 0.062 | 47.6% |
| P1 — Zero-shot | 0.060 | 49.8% |
| P2 — Constrained | 0.060 | 49.1% |

### FUNCTION Text Ablation

| Condition | Mistral Large 123B | Llama 3.1 8B | Gemma3 12B |
|---|---|---|---|
| A1 — Name + organism only | 0.119 | 0.084 | 0.097 |
| A2 — Name + first sentence | 0.124 | 0.092 | 0.105 |
| A3 — Full FUNCTION text | 0.125 | 0.097 | 0.100 |
| Drop A3→A1 (relative) | 5.0% | 13.4% | 2.9% |

Performance drop when FUNCTION text is removed is small (2.9–13.4%), confirming that PROBE measures **genuine pre-training biological knowledge**, not text extraction ability.

### Statistical Results

- **42/45** model pairs significantly different (Wilcoxon, Bonferroni-corrected, p<0.05)
- Effect sizes: d = 0.20–0.76 (negligible to medium)
- **10/10** models: namespace, organism, and prompt format effects significant (p<0.001)
- P5 is best prompt for **7/10** models
- Mistral Large 123B and Llama 3.3 70B are **statistically indistinguishable** (p=0.511) — both treated as co-best
- P3 exemplars (TP53, EGFR, ACT1) verified absent from test set — **no benchmark contamination**

---

## Citation

If you use PROBE in your research, please cite:

```bibtex
@article{naha2026probe,
  title={{PROBE}: A Cross-Organism Benchmark for Evaluating General-Purpose
         Large Language Models on Gene Ontology Protein Function Prediction},
  author={Naha, Kallol and Jamil, Hasan M.},
  journal={BMC Bioinformatics},
  year={2026},
  publisher={BioMed Central}
}
```

---

## License

MIT License — see [LICENSE](LICENSE) for details.
Benchmark dataset derived from UniProt Swiss-Prot ([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)).

---

## Contact

| | |
|---|---|
| **Kallol Naha** | naha7197@vandals.uidaho.edu |
| **Hasan M. Jamil** | jamil@uidaho.edu |
| **Institution** | Department of Computer Science, University of Idaho |
| **Issues** | Please open a GitHub issue for bugs or questions |
