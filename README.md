# PROBE: A Cross-Organism Benchmark for Evaluating General-Purpose LLMs on Zero-Shot Gene Ontology Protein Function Prediction

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![BMC Bioinformatics](https://img.shields.io/badge/journal-BMC%20Bioinformatics-green.svg)](https://bmcbioinformatics.biomedcentral.com/)

> **Kallol Naha, Hasan M. Jamil**  
> Department of Computer Science, University of Idaho  
> naha7197@vandals.uidaho.edu | jamil@uidaho.edu

---

## Overview

**PROBE** (**PR**otein functi**O**n **B**enchmark **E**valuation) is the first systematic benchmark for evaluating whether general-purpose open-source Large Language Models (LLMs) can predict Gene Ontology (GO) protein function terms in a **zero-shot, inference-only setting** — no fine-tuning, no retrieval augmentation.

We evaluate **10 open-source LLMs** across **5 prompt formats**, **5 organisms**, and **3 GO namespaces**, yielding **150,000 API calls** on a curated benchmark of **1,000 proteins** from UniProt Swiss-Prot with experimental evidence-only annotations.

---

## Key Findings

| Finding | Result |
|---|---|
| Best model (Mistral Large 123B) | Macro-F1 = **9.3%** — far below specialized tools |
| Hallucination rate | **30.6% – 50.4%** (mean 41.3%) |
| Best organism | *S. cerevisiae* (yeast) > *H. sapiens* — literature density effect |
| Worst organism | *D. rerio* (zebrafish) — avg F1 = 2.6% |
| Best prompt format | P5 candidate selection — reduces hallucination by **33.5 pp** |
| Model pairs significantly different | **42 / 45** (Wilcoxon, Bonferroni-corrected) |
| Dominant failure mode | Complete miss **68.71%** — not hallucination (8.65%) |
| Fully correct predictions | Only **0.34%** |

---

## Repository Structure

```
PROBE/
├── benchmark_dataset/
│   └── benchmark_master.json          # 1,000 proteins, 5 organisms, experimental GO only
│
├── scripts/
│   ├── dataset_preparation.py         # Swiss-Prot parsing, filtering, stratified sampling
│   ├── llm_prompting.py               # API calls, 5 prompt formats, auto-resume
│   ├── evaluation.py                  # F1, Precision, Recall, Hallucination Rate
│   ├── statistical_tests.py           # Wilcoxon, Friedman, Kruskal-Wallis, Cohen's d
│   ├── error_analysis.py              # 5-category error classification
│   └── visualization.py              # 7 publication-ready figures
│
├── probe_evaluation/
│   ├── leaderboard.csv                # Overall model rankings
│   ├── model_x_namespace.csv          # F1 by model × GO namespace
│   ├── model_x_organism.csv           # F1 by model × organism
│   ├── model_x_prompt.csv             # F1 by model × prompt format
│   ├── organism_x_namespace.csv       # F1 by organism × namespace
│   ├── full_detail.csv                # Full breakdown
│   └── PROBE_results.xlsx             # All tables in one Excel workbook
│
├── probe_figures/
│   ├── fig0_pipeline.png              # Pipeline overview diagram
│   ├── fig1_leaderboard.pdf           # Macro-F1 leaderboard bar chart
│   ├── fig2_namespace_heatmap.pdf     # F1 heatmap by model and namespace
│   ├── fig3_organism_curve.pdf        # Cross-organism F1 line plot
│   ├── fig4_prompt_heatmap.pdf        # Prompt sensitivity heatmap
│   ├── fig5_hallucination.pdf         # Hallucination rate bar chart
│   ├── fig6_precision_recall.pdf      # Precision-Recall scatter
│   └── fig7_f1_boxplot.pdf            # Per-protein F1 distributions
│
├── probe_error_analysis/              # Error categorization results
├── probe_statistics/                  # Statistical test results
│
├── requirements.txt
├── .gitignore
└── README.md
```
> **Note:** The manuscript and LaTeX source will be released upon journal acceptance.
> Raw LLM response files (`probe_results/`, ~400MB) are not included due to size constraints
> but are available from the corresponding author upon request.

---

## Dataset

**Source:** UniProt Swiss-Prot release 2024\_03 (574,627 reviewed entries)

**Construction pipeline (`dataset_preparation.py`):**

1. **Organism filter** — 5 target organisms:
   - *Homo sapiens*, *Mus musculus*, *Danio rerio*, *Saccharomyces cerevisiae*, *Escherichia coli*

2. **Evidence filter** — experimental codes only:
   - IDA, IMP, IPI, IGI, IEP, EXP, HDA, HMP, HGI, HEP
   - IEA (electronic annotation) **excluded** — prevents circular reasoning

3. **Depth filter** — minimum 3 experimental GO terms per protein

4. **Stratified sampling** — 200 proteins per organism (seed=42) = **1,000 total**

**Statistics:**

| Organism | N | Pool | Avg GO terms |
|---|---|---|---|
| *H. sapiens* | 200 | 11,042 | 7.7 |
| *M. musculus* | 200 | 8,758 | 7.7 |
| *D. rerio* | 200 | 759 | 7.7 |
| *S. cerevisiae* | 200 | 4,560 | 7.7 |
| *E. coli* | 200 | 2,177 | 7.7 |
| **Total** | **1,000** | **27,296** | **7.7** |

Overall: MF 68% · BP 94% · CC 76% · max 64 GO terms · min 3

---

## Models Evaluated

| Rank | Model | Family | Size | Architecture |
|---|---|---|---|---|
| 1 | Mistral Large 123B | Mistral | 123B | Dense |
| 2 | Llama 3.3 70B | Llama 3 | 70B | Dense |
| 3 | Qwen2.5 72B | Qwen 2.5 | 72B | Dense |
| 4 | Qwen3 32B | Qwen 3 | 32B | Dense |
| 5 | Gemma-4 31B | Gemma | 31B | Dense |
| 6 | Llama 3.1 8B | Llama 3 | 8B | Dense |
| 7 | Mixtral 8×7B | Mistral | 56B (14B active) | MoE |
| 8 | Qwen2.5 7B | Qwen 2.5 | 7B | Dense |
| 9 | Mistral 7B | Mistral | 7B | Dense |
| 10 | Gemma3 12B | Gemma | 12B | Dense |

All evaluated at T=0, max 512 tokens, no system prompt, inference only.

---

## Prompt Formats

| ID | Name | Description |
|---|---|---|
| P1 | Zero-shot | Free GO term prediction with no examples or constraints |
| P2 | Constrained | Namespace-restricted zero-shot with strict format rules |
| P3 | Few-shot (3-shot) | 3 annotated examples: TP53, EGFR, ACT1 |
| P4 | Chain-of-thought | Step-by-step biological reasoning before prediction |
| P5 | Candidate selection | Select from 15 pre-validated GO terms per namespace |

Each protein is queried 5 prompts × 3 namespaces = 15 times per model.  
Total: 1,000 × 15 × 10 = **150,000 API calls**.

---

## Reproduce Results

### Step 0 — Install dependencies

```bash
pip install -r requirements.txt
```

### Step 1 — Build dataset

Download Swiss-Prot flat file from https://www.uniprot.org/downloads  
Place as `uniprot_sprot.dat` in the PROBE root folder, then:

```bash
python scripts/dataset_preparation.py
```

Output: `benchmark_dataset/benchmark_master.json`

> **Note:** The pre-built dataset is already included in this repository.  
> You only need this step if you want to rebuild from scratch.

---

### Step 2 — Run LLM Inference

#### Option A — University of Idaho MindRouter (internal access)

```bash
python scripts/llm_prompting.py \
  --api_key YOUR_MINDROUTER_KEY \
  --model mistral-large:123b
```

Run all 10 models:

```bash
python scripts/llm_prompting.py \
  --api_key YOUR_MINDROUTER_KEY \
  --all_models
```

#### Option B — Ollama (public, local inference)

Researchers without MindRouter access can reproduce results using [Ollama](https://ollama.ai), which runs the same models locally with an OpenAI-compatible API.

**Install Ollama:**

```bash
# macOS / Linux
curl -fsSL https://ollama.ai/install.sh | sh
```

**Pull models:**

```bash
ollama pull llama3.3:70b
ollama pull llama3.1:8b
ollama pull mistral:7b
ollama pull qwen2.5:7b
ollama pull qwen2.5:72b
ollama pull gemma3:12b
ollama pull mixtral:8x7b
```

> Note: Qwen3-32B and Gemma-4-31B may be available as `qwen3:32b` and `gemma4:31b` depending on your Ollama version.

**Run PROBE with Ollama:**

```bash
python scripts/llm_prompting.py \
  --api_key ollama \
  --base_url http://localhost:11434/v1 \
  --model llama3.1:8b
```

> **Hardware requirements:** 7B models run on 8GB VRAM. 70B models require 40GB+ VRAM or CPU offloading. 

#### Option C — Any OpenAI-compatible API

```bash
python scripts/llm_prompting.py \
  --api_key YOUR_KEY \
  --base_url https://your-api-endpoint/v1 \
  --model your-model-name
```

#### Auto-resume

The script automatically resumes from where it left off if interrupted. Only successful calls are counted — failed calls (timeouts, connection errors) are retried on the next run.

---

### Step 3 — Evaluate

```bash
python scripts/evaluation.py
```

Output: `probe_evaluation/` — all CSV tables and Excel workbook.

Computes: F1, Precision, Recall, Hallucination Rate, Empty Response Rate  
per model × organism × namespace × prompt format.

---

### Step 4 — Statistical Tests

```bash
python scripts/statistical_tests.py
```

Output: `probe_statistics/`

| Test | Purpose |
|---|---|
| Wilcoxon signed-rank (45 pairs) | Pairwise model comparison |
| Friedman test (per model) | Prompt format significance |
| Kruskal-Wallis (per model) | Organism effect significance |
| Kruskal-Wallis (per model) | Namespace effect significance |
| Cohen's d | Effect size for significant pairs |

---

### Step 5 — Error Analysis

```bash
python scripts/error_analysis.py
```

Output: `probe_error_analysis/`

Error categories:
- **Complete miss** — valid GO IDs predicted but none match ground truth
- **Complete hallucination** — all predicted IDs absent from reference set
- **Partial match** (poor / high-recall / high-precision / good)
- **Empty response** — no GO IDs extracted

---

### Step 6 — Generate Figures

```bash
python scripts/visualization.py
```

Output: `probe_figures/` — 7 figures as PDF + PNG.

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
| 7 | Mixtral 8×7B | 0.059 | 0.049 | 0.109 | 45.7 | 2.8 |
| 8 | Qwen2.5 7B | 0.044 | 0.049 | 0.050 | 47.4 | 5.8 |
| 9 | Mistral 7B | 0.042 | 0.031 | 0.087 | 49.0 | 0.0 |
| 10 | Gemma3 12B | 0.032 | 0.028 | 0.053 | 50.4 | 0.7 |

HR = Hallucination Rate · ERR = Empty Response Rate

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

### Error Distribution (all models combined)

| Error Category | % |
|---|---|
| Complete miss | 68.71% |
| Complete hallucination | 8.65% |
| Partial match (various) | 19.10% |
| Empty response | 3.20% |
| **Correct** | **0.34%** |

---

## Statistical Results

- **42/45** model pairs significantly different (Wilcoxon, Bonferroni-corrected, p<0.05)
- Effect sizes: d = 0.20–0.76 (negligible to medium)
- **10/10** models: namespace effect significant (Kruskal-Wallis, p<0.001)
- **10/10** models: organism effect significant (Kruskal-Wallis, p<0.001)
- **10/10** models: prompt format effect significant (Friedman, p<0.001)
- P5 is best prompt for **7/10** models

---

## Citation

If you use PROBE in your research, please cite:

```bibtex
@article{naha2025probe,
  title={{PROBE}: A Cross-Organism Benchmark for Evaluating
         General-Purpose Large Language Models on Zero-Shot
         Gene Ontology Protein Function Prediction},
  author={Naha, Kallol and Jamil, Hasan M.},
  journal={BMC Bioinformatics},
  year={2025},
  publisher={BioMed Central}
}
```

---

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.

The benchmark dataset is derived from UniProt Swiss-Prot, which is available under the [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/) license.

---

## Contact

| | |
|---|---|
| **Kallol Naha** | naha7197@vandals.uidaho.edu |
| **Hasan M. Jamil** | jamil@uidaho.edu |
| **Institution** | Department of Computer Science, University of Idaho |
| **Issues** | Please open a GitHub issue for bugs or questions |
