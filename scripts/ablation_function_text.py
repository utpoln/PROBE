"""
PROBE: FUNCTION Text Ablation Study
Tests whether performance comes from biological reasoning or FUNCTION text extraction.

Three input conditions per protein:
  A1 - Protein name + organism only (no FUNCTION text)
  A2 - Protein name + organism + partial hint (first sentence of FUNCTION only)
  A3 - Full context (current PROBE setup — gene name + full FUNCTION text)

Usage:
  python scripts/ablation_function_text.py \
    --api_key YOUR_KEY \
    --base_url http://localhost:11434/v1 \
    --model mistral-large:123b \
    --limit 100
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

DATASET     = "benchmark_dataset/benchmark_master.json"
OUTPUT_DIR  = Path("probe_results_ablation")
SLEEP_SEC   = 0.5
TEMPERATURE = 0.0
MAX_TOKENS  = 512
N_PROTEINS  = 100  # 100 proteins is enough for ablation

NAMESPACES = ["molecular_function", "biological_process", "cellular_component"]
NAMESPACE_LABELS = {
    "molecular_function": "Molecular Function",
    "biological_process": "Biological Process",
    "cellular_component": "Cellular Component",
}

# Use P1 zero-shot only for clean comparison
# Same prompt structure, only input changes
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
# PROMPT BUILDERS — 3 ablation conditions
# ─────────────────────────────────────────────

def build_prompts(entry: dict, namespace: str) -> dict:
    gene        = entry.get("gene_name", "Unknown")
    protein     = entry.get("protein_name", "").split("{")[0].strip()
    org         = entry.get("organism_full", "")
    func        = entry.get("function_text", "")
    ns_label    = NAMESPACE_LABELS[namespace]
    candidates  = CANDIDATES[namespace]

    # First sentence of FUNCTION text only
    first_sentence = func.split(".")[0].strip() + "." if func else ""

    # ── A1: Name + organism only (no FUNCTION text) ──────────────
    context_A1 = (
        f"Protein: {gene}\n"
        f"Full name: {protein}\n"
        f"Organism: {org}"
    )

    # ── A2: Name + organism + first sentence only ─────────────────
    context_A2 = (
        f"Protein: {gene}\n"
        f"Full name: {protein}\n"
        f"Organism: {org}\n"
        f"Known function: {first_sentence}"
    )

    # ── A3: Full context (current PROBE setup) ────────────────────
    context_A3 = (
        f"Protein: {gene}\n"
        f"Full name: {protein}\n"
        f"Organism: {org}\n"
        f"Known function: {func}"
    )

    prompts = {}

    # P5 format for all three conditions
    for cond, context in [("A1_name_only", context_A1),
                           ("A2_first_sentence", context_A2),
                           ("A3_full_function", context_A3)]:
        prompts[cond] = (
            f"You are an expert molecular biologist.\n\n"
            f"{context}\n\n"
            f"From the following GO {ns_label} terms, select ALL that apply.\n"
            f"Candidate terms:\n{candidates}\n\n"
            f"List only the matching terms as: GO:XXXXXXX | term name\n"
            f"If none apply, write: NONE"
        )

    return prompts

# ─────────────────────────────────────────────
# EVALUATION
# ─────────────────────────────────────────────

def extract_go_terms(text: str) -> set:
    return set(re.findall(r"GO:\d{7}", text))

def compute_f1(predicted: set, ground_truth: set):
    if not predicted or not ground_truth:
        return 0.0, 0.0, 0.0
    tp   = len(predicted & ground_truth)
    prec = tp / len(predicted) if predicted else 0.0
    rec  = tp / len(ground_truth) if ground_truth else 0.0
    f1   = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0.0
    return round(prec, 4), round(rec, 4), round(f1, 4)

# ─────────────────────────────────────────────
# API CALLER
# ─────────────────────────────────────────────

def call_model(client: OpenAI, model: str, prompt: str) -> dict:
    try:
        response = client.chat.completions.create(
            model=model,
            messages=[{"role": "user", "content": prompt}],
            max_tokens=MAX_TOKENS,
            temperature=TEMPERATURE,
        )
        content = response.choices[0].message.content or ""
        return {"success": True, "response": content}
    except Exception as e:
        return {"success": False, "response": "", "error": str(e)}

# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

def run(api_key: str, model: str, base_url: str, limit: int):
    OUTPUT_DIR.mkdir(exist_ok=True)

    client = OpenAI(base_url=base_url, api_key=api_key)

    print(f"Loading dataset...")
    with open(DATASET) as f:
        entries = json.load(f)

    entries = entries[:limit]
    model_slug  = model.replace("/", "_").replace(":", "-").replace(".", "-")
    output_file = OUTPUT_DIR / f"ablation_{model_slug}.jsonl"

    # Load completed
    completed = set()
    if output_file.exists():
        with open(output_file) as f:
            for line in f:
                try:
                    rec = json.loads(line)
                    if rec.get("success"):
                        completed.add((rec["accession"], rec["namespace"], rec["condition"]))
                except:
                    continue

    print(f"Model      : {model}")
    print(f"Proteins   : {len(entries)}")
    print(f"Conditions : A1 (name only), A2 (first sentence), A3 (full function)")
    print(f"Namespaces : {len(NAMESPACES)}")
    print(f"Total calls: {len(entries) * 3 * len(NAMESPACES)}")
    print(f"Already done: {len(completed)}")
    print()

    # Accumulators per condition
    results = {
        "A1_name_only":      {"prec": [], "rec": [], "f1": []},
        "A2_first_sentence": {"prec": [], "rec": [], "f1": []},
        "A3_full_function":  {"prec": [], "rec": [], "f1": []},
    }

    with open(output_file, "a") as out_f:
        for entry in tqdm(entries, desc="Ablation"):
            accession = entry["accession"]
            for namespace in NAMESPACES:
                prompts = build_prompts(entry, namespace)
                gt_raw  = entry.get(f"go_{namespace}", "")
                gt_terms = extract_go_terms(gt_raw)

                for condition, prompt_text in prompts.items():
                    key = (accession, namespace, condition)
                    if key in completed:
                        # Still accumulate from file
                        continue

                    result = call_model(client, model, prompt_text)
                    pred_terms = extract_go_terms(result["response"])
                    prec, rec, f1 = compute_f1(pred_terms, gt_terms)

                    results[condition]["prec"].append(prec)
                    results[condition]["rec"].append(rec)
                    results[condition]["f1"].append(f1)

                    record = {
                        "accession":   accession,
                        "gene_name":   entry.get("gene_name", ""),
                        "organism":    entry.get("organism", ""),
                        "namespace":   namespace,
                        "condition":   condition,
                        "model":       model,
                        "response":    result["response"],
                        "success":     result["success"],
                        "error":       result.get("error", ""),
                        "ground_truth": gt_raw,
                        "predicted_go": list(pred_terms),
                        "gt_go":        list(gt_terms),
                        "precision":   prec,
                        "recall":      rec,
                        "f1":          f1,
                        "timestamp":   datetime.utcnow().isoformat(),
                    }

                    out_f.write(json.dumps(record) + "\n")
                    out_f.flush()
                    time.sleep(SLEEP_SEC)

    # ── Print summary ─────────────────────────────────────────────
    print("\n" + "="*65)
    print(f"FUNCTION TEXT ABLATION RESULTS — {model}")
    print("="*65)
    print(f"{'Condition':<25} {'Prec':>8} {'Recall':>8} {'F1':>8}")
    print("-"*65)

    condition_labels = {
        "A1_name_only":      "A1 — Name + organism only",
        "A2_first_sentence": "A2 — Name + first sentence",
        "A3_full_function":  "A3 — Full FUNCTION text",
    }

    for cond, label in condition_labels.items():
        r = results[cond]
        if r["f1"]:
            p = sum(r["prec"]) / len(r["prec"])
            rc = sum(r["rec"]) / len(r["rec"])
            f = sum(r["f1"]) / len(r["f1"])
            print(f"{label:<25} {p:>8.4f} {rc:>8.4f} {f:>8.4f}")
        else:
            print(f"{label:<25}  (loaded from file — rerun to recompute)")

    print("="*65)
    print(f"\nResults saved: {output_file}")
    print("\nKey question: How much does F1 drop from A3 → A1?")
    print("If drop is small  → LLMs have genuine biological knowledge")
    print("If drop is large  → performance depends on FUNCTION text extraction")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--api_key",  required=True)
    parser.add_argument("--model",    required=True)
    parser.add_argument("--base_url", default="https://mindrouter.uidaho.edu/v1")
    parser.add_argument("--limit",    type=int, default=100)
    args = parser.parse_args()

    run(
        api_key=args.api_key,
        model=args.model,
        base_url=args.base_url,
        limit=args.limit,
    )