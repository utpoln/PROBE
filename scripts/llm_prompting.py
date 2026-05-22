"""
PROBE: LLM Prompting Pipeline
Cross-Organism Evaluation of General-Purpose LLMs for GO-Based
Protein Function Prediction

Usage:
  python llm_prompting.py --api_key YOUR_KEY --model phi4-reasoning:14b
  python llm_prompting.py --api_key YOUR_KEY --all_models
  python llm_prompting.py --api_key YOUR_KEY --all_models --limit 50
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

BASE_URL    = "https://mindrouter.uidaho.edu/v1"
DATASET     = "benchmark_dataset/benchmark_master.json"
OUTPUT_DIR  = Path("probe_results")
SLEEP_SEC   = 0.4
TEMPERATURE = 0.0

# Default max tokens — increased for reasoning models
MAX_TOKENS_DEFAULT  = 512
MAX_TOKENS_THINKING = 4096   # phi4-reasoning needs full thinking chain to complete

# Models that use thinking/reasoning mode — need special handling
THINKING_MODELS = {
    "phi4-reasoning:14b",
    "Qwen/Qwen3-32B",
    "qwen3:32b",
}

# 10 benchmark models
MODELS = [
    "llama3.3:70b",
    "llama3.1:8b",
    "mistral-large:123b",
    "mistral:7b",
    "qwen2.5:72b",
    "qwen2.5:7b",
    "Qwen/Qwen3-32B",
    "gemma3:12b",
    "phi4-reasoning:14b",
    "mixtral:8x7b",
]

NAMESPACES = ["molecular_function", "biological_process", "cellular_component"]
NAMESPACE_LABELS = {
    "molecular_function": "Molecular Function",
    "biological_process": "Biological Process",
    "cellular_component": "Cellular Component",
}

# ─────────────────────────────────────────────
# PROMPT BUILDERS
# ─────────────────────────────────────────────

def build_prompts(entry: dict, namespace: str, thinking_model: bool = False) -> dict:
    """Build all 5 prompt formats. Adds /no_think for reasoning models."""
    gene         = entry.get("gene_name", "Unknown")
    protein      = entry.get("protein_name", "")
    org          = entry.get("organism_full", "")
    func         = entry.get("function_text", "")
    ns_label     = NAMESPACE_LABELS[namespace]
    protein_short = protein.split("{")[0].strip()

    context = f"Protein: {gene}\nFull name: {protein_short}\nOrganism: {org}"
    if func:
        context += f"\nKnown function: {func}"

    # No prefix needed — /no_think does not work for phi4-reasoning
    # Instead we use max_tokens=4096 to let thinking complete,
    # then strip <think>...</think> in extract_response_content()
    prefix = ""

    prompts = {}

    prompts["P1_zeroshot"] = (
        f"{prefix}"
        f"You are an expert molecular biologist.\n\n"
        f"{context}\n\n"
        f"List the Gene Ontology (GO) {ns_label} terms that best describe "
        f"this protein. Provide GO IDs and term names. "
        f"Format each as: GO:XXXXXXX | term name"
    )

    prompts["P2_constrained"] = (
        f"{prefix}"
        f"You are an expert molecular biologist specializing in protein annotation.\n\n"
        f"{context}\n\n"
        f"Task: Predict ONLY the Gene Ontology {ns_label} "
        f"(GO {namespace.upper()[:2]}) terms for this protein.\n"
        f"Rules:\n"
        f"- List only experimentally supported functions\n"
        f"- Format each term as: GO:XXXXXXX | term name\n"
        f"- Do not include Biological Process or Cellular Component terms\n"
        f"- Output only the GO terms, nothing else"
    )

    few_shot_examples = _get_fewshot_examples(namespace)
    prompts["P3_fewshot"] = (
        f"{prefix}"
        f"You are an expert molecular biologist. "
        f"Predict Gene Ontology {ns_label} terms for proteins.\n\n"
        f"Here are three examples:\n\n"
        f"{few_shot_examples}\n\n"
        f"Now predict for:\n"
        f"{context}\n\n"
        f"List GO {ns_label} terms. Format: GO:XXXXXXX | term name"
    )

    prompts["P4_cot"] = (
        f"{prefix}"
        f"You are an expert molecular biologist.\n\n"
        f"{context}\n\n"
        f"Think step by step about what {ns_label} GO terms apply to this protein:\n"
        f"Step 1: What is the protein's primary biochemical activity?\n"
        f"Step 2: What molecular processes does it participate in?\n"
        f"Step 3: Based on your reasoning, list the GO {ns_label} terms.\n\n"
        f"Format final answer as: GO:XXXXXXX | term name"
    )

    candidates = _get_candidate_terms(namespace)
    prompts["P5_selection"] = (
        f"{prefix}"
        f"You are an expert molecular biologist.\n\n"
        f"{context}\n\n"
        f"From the following GO {ns_label} terms, select ALL that apply.\n"
        f"Candidate terms:\n{candidates}\n\n"
        f"List only the matching terms as: GO:XXXXXXX | term name\n"
        f"If none apply, write: NONE"
    )

    return prompts


def _get_fewshot_examples(namespace: str) -> str:
    examples = {
        "molecular_function": (
            "Example 1:\nProtein: TP53, Organism: Homo sapiens\n"
            "GO Molecular Function terms:\n"
            "GO:0003700 | DNA-binding transcription factor activity\n"
            "GO:0046872 | metal ion binding\n"
            "GO:0042802 | identical protein binding\n\n"
            "Example 2:\nProtein: EGFR, Organism: Homo sapiens\n"
            "GO Molecular Function terms:\n"
            "GO:0004714 | transmembrane receptor protein tyrosine kinase activity\n"
            "GO:0005515 | protein binding\n"
            "GO:0016301 | kinase activity\n\n"
            "Example 3:\nProtein: ACT1, Organism: Saccharomyces cerevisiae\n"
            "GO Molecular Function terms:\n"
            "GO:0005524 | ATP binding\n"
            "GO:0003779 | actin binding\n"
            "GO:0005198 | structural molecule activity"
        ),
        "biological_process": (
            "Example 1:\nProtein: TP53, Organism: Homo sapiens\n"
            "GO Biological Process terms:\n"
            "GO:0006915 | apoptotic process\n"
            "GO:0006351 | DNA-templated transcription\n"
            "GO:0007050 | cell cycle arrest\n\n"
            "Example 2:\nProtein: BRCA1, Organism: Homo sapiens\n"
            "GO Biological Process terms:\n"
            "GO:0006281 | DNA repair\n"
            "GO:0007050 | cell cycle arrest\n"
            "GO:0045786 | negative regulation of cell cycle\n\n"
            "Example 3:\nProtein: CDC28, Organism: Saccharomyces cerevisiae\n"
            "GO Biological Process terms:\n"
            "GO:0007049 | cell cycle\n"
            "GO:0006468 | protein phosphorylation\n"
            "GO:0045736 | negative regulation of cyclin-dependent protein kinase activity"
        ),
        "cellular_component": (
            "Example 1:\nProtein: TP53, Organism: Homo sapiens\n"
            "GO Cellular Component terms:\n"
            "GO:0005654 | nucleoplasm\n"
            "GO:0005829 | cytosol\n"
            "GO:0000785 | chromatin\n\n"
            "Example 2:\nProtein: MYH9, Organism: Homo sapiens\n"
            "GO Cellular Component terms:\n"
            "GO:0005925 | focal adhesion\n"
            "GO:0015629 | actin cytoskeleton\n"
            "GO:0032587 | ruffle membrane\n\n"
            "Example 3:\nProtein: TUB1, Organism: Saccharomyces cerevisiae\n"
            "GO Cellular Component terms:\n"
            "GO:0005737 | cytoplasm\n"
            "GO:0005874 | microtubule\n"
            "GO:0000779 | condensed chromosome, centromeric region"
        ),
    }
    return examples.get(namespace, "")


def _get_candidate_terms(namespace: str) -> str:
    candidates = {
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
    return candidates.get(namespace, "")


# ─────────────────────────────────────────────
# RESPONSE PARSER
# ─────────────────────────────────────────────

def extract_response_content(choice, model: str) -> str:
    """
    Extract actual answer content from model response.
    Handles three cases:
      1. Normal models: content is the answer
      2. Qwen3: content=None, answer in reasoning_content
      3. phi4-reasoning: content has <think>...</think> wrapping the answer
         Strip <think> block and return only the final answer after </think>
    """
    content = choice.message.content

    # Case 2: Qwen3-style — content is None
    if content is None:
        content = getattr(choice.message, "reasoning_content", None) or ""

    # Case 3: phi4-reasoning — strip <think>...</think> block
    # The actual GO term answer comes AFTER </think>
    if content and "<think>" in content:
        # Try to extract content after </think>
        after_think = re.split(r"</think>", content, flags=re.IGNORECASE)
        if len(after_think) > 1:
            # Take everything after the last </think>
            content = after_think[-1].strip()
        else:
            # </think> not found — model was cut off inside thinking
            # Fall back: look for GO IDs anywhere in the response
            # including inside <think> block (better than nothing)
            content = content  # keep full content, GO parser will find IDs

    return content or ""


# ─────────────────────────────────────────────
# API CALLER
# ─────────────────────────────────────────────

def call_model(client: OpenAI, model: str, prompt: str,
               is_thinking: bool = False) -> dict:
    """Call MindRouter API and return response dict."""
    max_tokens = MAX_TOKENS_THINKING if is_thinking else MAX_TOKENS_DEFAULT
    try:
        response = client.chat.completions.create(
            model=model,
            messages=[{"role": "user", "content": prompt}],
            max_tokens=max_tokens,
            temperature=TEMPERATURE,
        )
        content = extract_response_content(response.choices[0], model)
        return {
            "success":  True,
            "response": content,
            "model":    model,
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
            "model":    model,
        }


# ─────────────────────────────────────────────
# RESUME LOGIC
# ─────────────────────────────────────────────

def load_completed(output_file: Path) -> set:
    """Return only SUCCESSFUL calls — failed records will be retried."""
    completed = set()
    if not output_file.exists():
        return completed
    with open(output_file) as f:
        for line in f:
            try:
                rec = json.loads(line)
                # Only skip if the call actually succeeded
                # Failed records (connection errors, timeouts) get retried
                if rec.get("success", False):
                    completed.add((rec["accession"], rec["namespace"], rec["prompt_id"]))
            except Exception:
                continue
    return completed


# ─────────────────────────────────────────────
# MAIN RUNNER
# ─────────────────────────────────────────────

def run(api_key: str, models: list, limit: int = None):
    OUTPUT_DIR.mkdir(exist_ok=True)

    print(f"Loading dataset from {DATASET} ...")
    with open(DATASET) as f:
        entries = json.load(f)

    if limit:
        entries = entries[:limit]
        print(f"Running in LIMIT mode: {limit} proteins only")

    print(f"Total proteins     : {len(entries):,}")
    print(f"Models to run      : {len(models)}")
    print(f"Prompts per protein: 5")
    print(f"Namespaces         : {len(NAMESPACES)}")
    print(f"Total API calls    : {len(entries)*5*len(NAMESPACES)*len(models):,}\n")

    client = OpenAI(base_url=BASE_URL, api_key=api_key)

    for model in models:
        is_thinking = model in THINKING_MODELS
        model_slug  = model.replace("/", "_").replace(":", "-")
        output_file = OUTPUT_DIR / f"{model_slug}.jsonl"
        completed   = load_completed(output_file)

        print(f"\n{'='*60}")
        print(f"Model: {model}")
        print(f"Output: {output_file}")
        print(f"Thinking mode: {'YES (max_tokens=4096, </think> stripper active)' if is_thinking else 'NO (max_tokens=512)'}")
        print(f"Already completed: {len(completed):,} calls")
        print(f"{'='*60}")

        with open(output_file, "a") as out_f:
            for entry in tqdm(entries, desc=model_slug):
                accession = entry["accession"]
                for namespace in NAMESPACES:
                    prompts = build_prompts(entry, namespace,
                                           thinking_model=is_thinking)
                    for prompt_id, prompt_text in prompts.items():
                        key = (accession, namespace, prompt_id)
                        if key in completed:
                            continue

                        result = call_model(client, model, prompt_text,
                                           is_thinking=is_thinking)

                        record = {
                            "accession":    accession,
                            "gene_name":    entry.get("gene_name", ""),
                            "organism":     entry.get("organism", ""),
                            "namespace":    namespace,
                            "prompt_id":    prompt_id,
                            "prompt_text":  prompt_text,
                            "model":        model,
                            "response":     result["response"],
                            "success":      result["success"],
                            "error":        result.get("error", ""),
                            "usage":        result.get("usage", {}),
                            "ground_truth": entry.get(f"go_{namespace}", ""),
                            "timestamp":    datetime.utcnow().isoformat(),
                        }

                        out_f.write(json.dumps(record) + "\n")
                        out_f.flush()
                        time.sleep(SLEEP_SEC)

    print("\n✅ Prompting complete!")
    print(f"   Results saved in: {OUTPUT_DIR}/")


# ─────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--api_key",    required=True)
    parser.add_argument("--model",      default=None)
    parser.add_argument("--all_models", action="store_true")
    parser.add_argument("--limit",      type=int, default=None)
    args = parser.parse_args()

    if args.all_models:
        selected = MODELS
    elif args.model:
        selected = [args.model]
    else:
        print("Please specify --model MODEL_NAME or --all_models")
        exit(1)

    run(api_key=args.api_key, models=selected, limit=args.limit)
