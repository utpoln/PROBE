"""
PROBE: Closed-Source Model Evaluation
P5 candidate selection only, 50 proteins, 3 namespaces
Usage:
  python scripts/closed_source_eval.py --api_key YOUR_KEY --model gpt-5.4-2026-03-05 --provider openai
  python scripts/closed_source_eval.py --api_key YOUR_KEY --model claude-sonnet-4-5 --provider anthropic
"""

import json
import re
import time
import argparse
from pathlib import Path
from datetime import datetime
from openai import OpenAI
from tqdm import tqdm

# ─────────────────────────────────────────────
# CONFIG
# ─────────────────────────────────────────────

DATASET    = "benchmark_dataset/benchmark_master.json"
OUTPUT_DIR = Path("probe_results_closed")
SLEEP_SEC  = 0.5
TEMPERATURE = 0.0
MAX_TOKENS  = 512
N_PROTEINS  = 50

NAMESPACES = ["molecular_function", "biological_process", "cellular_component"]
NAMESPACE_LABELS = {
    "molecular_function": "Molecular Function",
    "biological_process": "Biological Process",
    "cellular_component": "Cellular Component",
}

PROVIDER_URLS = {
    "openai":    "https://api.openai.com/v1",
    "anthropic": "https://api.anthropic.com/v1",
    "google":    "https://generativelanguage.googleapis.com/v1beta/openai",
    "xai":       "https://api.x.ai/v1",
}

# ─────────────────────────────────────────────
# P5 CANDIDATE TERMS
# ─────────────────────────────────────────────

CANDIDATES = {
    "molecular_function": (
        "GO:0003700 | DNA-binding transcription factor activity\n"
        "GO:0004672 | protein kinase activity\n"
        "GO:0005515 | protein binding\n"
        "GO:0005524 | ATP binding\n"
        "GO:0003677 | DNA binding\n"
        "GO:0046872 | metal ion binding\n"
        "GO:0004722 | protein serine/threonine phosphatase activity\n"
        "GO:0016301 | kinase activity\n"
        "GO:0003723 | RNA binding\n"
        "GO:0004725 | protein tyrosine phosphatase activity\n"
        "GO:0003779 | actin binding\n"
        "GO:0005198 | structural molecule activity\n"
        "GO:0008270 | zinc ion binding\n"
        "GO:0016787 | hydrolase activity\n"
        "GO:0016740 | transferase activity"
    ),
    "biological_process": (
        "GO:0006915 | apoptotic process\n"
        "GO:0007049 | cell cycle\n"
        "GO:0006281 | DNA repair\n"
        "GO:0006351 | DNA-templated transcription\n"
        "GO:0006468 | protein phosphorylation\n"
        "GO:0007165 | signal transduction\n"
        "GO:0006355 | regulation of DNA-templated transcription\n"
        "GO:0008283 | cell population proliferation\n"
        "GO:0006954 | inflammatory response\n"
        "GO:0045944 | positive regulation of transcription by RNA polymerase II\n"
        "GO:0006412 | translation\n"
        "GO:0006810 | transport\n"
        "GO:0000077 | DNA damage checkpoint signaling\n"
        "GO:0045786 | negative regulation of cell cycle\n"
        "GO:0006986 | response to unfolded protein"
    ),
    "cellular_component": (
        "GO:0005654 | nucleoplasm\n"
        "GO:0005829 | cytosol\n"
        "GO:0005737 | cytoplasm\n"
        "GO:0005634 | nucleus\n"
        "GO:0005886 | plasma membrane\n"
        "GO:0005739 | mitochondrion\n"
        "GO:0005783 | endoplasmic reticulum\n"
        "GO:0005615 | extracellular space\n"
        "GO:0016020 | membrane\n"
        "GO:0005874 | microtubule\n"
        "GO:0015629 | actin cytoskeleton\n"
        "GO:0005925 | focal adhesion\n"
        "GO:0000785 | chromatin\n"
        "GO:0005694 | chromosome\n"
        "GO:0070062 | extracellular exosome"
    ),
}

# ─────────────────────────────────────────────
# PROMPT BUILDER - P5 ONLY
# ─────────────────────────────────────────────

def build_p5_prompt(entry: dict, namespace: str) -> str:
    gene     = entry.get("gene_name", "Unknown")
    protein  = entry.get("protein_name", "").split("{")[0].strip()
    org      = entry.get("organism_full", "")
    func     = entry.get("function_text", "")
    ns_label = NAMESPACE_LABELS[namespace]

    context = f"Protein: {gene}\nFull name: {protein}\nOrganism: {org}"
    if func:
        context += f"\nKnown function: {func}"

    candidates = CANDIDATES[namespace]

    prompt = (
        f"You are an expert molecular biologist.\n\n"
        f"{context}\n\n"
        f"From the following GO {ns_label} terms, select ALL that apply.\n"
        f"Candidate terms:\n{candidates}\n\n"
        f"List only the matching terms as: GO:XXXXXXX | term name\n"
        f"If none apply, write: NONE"
    )
    return prompt

# ─────────────────────────────────────────────
# EVALUATION METRICS
# ─────────────────────────────────────────────

def compute_f1(predicted: set, ground_truth: set):
    if not predicted or not ground_truth:
        return 0.0, 0.0, 0.0
    tp = len(predicted & ground_truth)
    prec = tp / len(predicted) if predicted else 0.0
    rec  = tp / len(ground_truth) if ground_truth else 0.0
    f1   = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0.0
    return round(prec, 4), round(rec, 4), round(f1, 4)

def extract_go_terms(text: str) -> set:
    return set(re.findall(r"GO:\d{7}", text))

# ─────────────────────────────────────────────
# API CALLER
# ─────────────────────────────────────────────

def call_model(client: OpenAI, model: str, prompt: str) -> dict:
    try:
        # Gemini 2.5 Pro is a thinking model and needs more tokens
        if "gemini-2.5-pro" in model:
            max_tok = 8192
        else:
            max_tok = MAX_TOKENS

        kwargs = {
            "model": model,
            "messages": [{"role": "user", "content": prompt}],
            "max_completion_tokens": max_tok,
        }
        if "gpt-5.5" not in model:
            kwargs["temperature"] = TEMPERATURE

        response = client.chat.completions.create(**kwargs)
        content = response.choices[0].message.content or ""
        return {
            "success":  True,
            "response": content,
            "usage": {
                "prompt_tokens":     response.usage.prompt_tokens,
                "completion_tokens": response.usage.completion_tokens,
            } if response.usage else {},
        }
    except Exception as e:
        return {
            "success":  False,
            "response": "",
            "error":    str(e),
        }
# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

def run(api_key: str, model: str, provider: str):
    OUTPUT_DIR.mkdir(exist_ok=True)

    base_url = PROVIDER_URLS.get(provider, PROVIDER_URLS["openai"])
    client   = OpenAI(base_url=base_url, api_key=api_key)

    print(f"Loading dataset...")
    with open(DATASET) as f:
        entries = json.load(f)

    # Take first 50 proteins
    entries = entries[:N_PROTEINS]
    print(f"Proteins: {len(entries)}")
    print(f"Model: {model}")
    print(f"Provider: {provider}")
    print(f"Base URL: {base_url}")
    print(f"Total calls: {len(entries) * len(NAMESPACES)} (P5 only)")

    model_slug  = model.replace("/", "_").replace(":", "-").replace(".", "-")
    output_file = OUTPUT_DIR / f"{model_slug}.jsonl"

    # Load completed
    completed = set()
    if output_file.exists():
        with open(output_file) as f:
            for line in f:
                try:
                    rec = json.loads(line)
                    if rec.get("success"):
                        completed.add((rec["accession"], rec["namespace"]))
                except:
                    continue
    print(f"Already completed: {len(completed)}")

    # Summary accumulators
    total_f1    = []
    total_prec  = []
    total_rec   = []

    with open(output_file, "a") as out_f:
        for entry in tqdm(entries, desc=model_slug):
            accession = entry["accession"]
            for namespace in NAMESPACES:
                key = (accession, namespace)
                if key in completed:
                    continue

                prompt = build_p5_prompt(entry, namespace)
                result = call_model(client, model, prompt)

                # Get ground truth
                gt_raw    = entry.get(f"go_{namespace}", "")
                gt_terms  = extract_go_terms(gt_raw)
                pred_terms = extract_go_terms(result["response"])

                prec, rec, f1 = compute_f1(pred_terms, gt_terms)
                total_f1.append(f1)
                total_prec.append(prec)
                total_rec.append(rec)

                record = {
                    "accession":    accession,
                    "gene_name":    entry.get("gene_name", ""),
                    "organism":     entry.get("organism", ""),
                    "namespace":    namespace,
                    "prompt_id":    "P5_selection",
                    "model":        model,
                    "response":     result["response"],
                    "success":      result["success"],
                    "error":        result.get("error", ""),
                    "usage":        result.get("usage", {}),
                    "ground_truth": gt_raw,
                    "predicted_go": list(pred_terms),
                    "gt_go":        list(gt_terms),
                    "precision":    prec,
                    "recall":       rec,
                    "f1":           f1,
                    "timestamp":    datetime.utcnow().isoformat(),
                }

                out_f.write(json.dumps(record) + "\n")
                out_f.flush()
                time.sleep(SLEEP_SEC)

    # Print summary
    if total_f1:
        print("\n" + "="*50)
        print(f"Model: {model}")
        print(f"Proteins evaluated: {len(entries)}")
        print(f"Prompt format: P5 (candidate selection)")
        print(f"Mean Precision: {sum(total_prec)/len(total_prec):.4f}")
        print(f"Mean Recall:    {sum(total_rec)/len(total_rec):.4f}")
        print(f"Mean F1:        {sum(total_f1)/len(total_f1):.4f}")
        print(f"Results saved:  {output_file}")
        print("="*50)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--api_key",  required=True)
    parser.add_argument("--model",    required=True)
    parser.add_argument("--provider", default="openai",
                        choices=["openai", "anthropic", "google", "xai"],
                        help="API provider")
    args = parser.parse_args()
    run(api_key=args.api_key, model=args.model, provider=args.provider)