#!/usr/bin/env python3
"""
LLM-as-a-Judge evaluation script for Arena data.
Uses GPT-5.2 and GPT-4o-mini via OpenRouter to judge model battles.

Usage:
    source ~/.bash_profile  # Load OPENROUTER_API_KEY
    python3 scripts/llm_judge_eval.py

Progress is saved after each sample, so you can interrupt and resume.
"""

import os
import json
import time
import pandas as pd
import numpy as np
from pathlib import Path
import requests
from tqdm import tqdm

# OpenRouter config
OPENROUTER_API_KEY = os.environ.get('OPENROUTER_API_KEY')
OPENROUTER_URL = "https://openrouter.ai/api/v1/chat/completions"

MODELS = {
    "gpt-5.2": "openai/gpt-5.2",
    "gpt-4o-mini": "openai/gpt-4o-mini"
}

JUDGE_PROMPT = """You are an impartial judge evaluating two AI assistant responses to a user prompt.

## User Prompt:
{prompt}

## Response A:
{response_a}

## Response B:
{response_b}

## Task:
Evaluate which response is better overall considering helpfulness, accuracy, and quality.
Respond with JSON: {{"winner": "A" or "B"}}"""


def extract_text_from_content(content):
    """Extract text from nested content structure."""
    if content is None:
        return ""
    if isinstance(content, np.ndarray):
        content = content.tolist()
    if isinstance(content, str):
        return content
    if isinstance(content, list):
        texts = []
        for item in content:
            if isinstance(item, dict):
                if 'text' in item and item['text']:
                    texts.append(str(item['text']))
            elif isinstance(item, str):
                texts.append(item)
        return ' '.join(texts)
    if isinstance(content, dict) and 'text' in content:
        return str(content['text'])
    return ""


def get_prompt_and_responses(row):
    """Extract user prompt and both model responses from a row."""
    conv_a = row['conversation_a']
    conv_b = row['conversation_b']
    if isinstance(conv_a, np.ndarray):
        conv_a = conv_a.tolist()
    if isinstance(conv_b, np.ndarray):
        conv_b = conv_b.tolist()

    prompt = ""
    response_a = ""
    response_b = ""

    for msg in conv_a:
        if isinstance(msg, dict):
            role = msg.get('role', '')
            content = msg.get('content', '')
            text = extract_text_from_content(content)
            if role == 'user' and not prompt:
                prompt = text
            elif role == 'assistant' and not response_a:
                response_a = text

    for msg in conv_b:
        if isinstance(msg, dict):
            role = msg.get('role', '')
            content = msg.get('content', '')
            text = extract_text_from_content(content)
            if role == 'assistant' and not response_b:
                response_b = text

    return prompt, response_a, response_b


def call_judge(prompt, response_a, response_b, model_id, max_retries=3):
    """Call the LLM judge via OpenRouter with JSON structured output."""
    headers = {
        "Authorization": f"Bearer {OPENROUTER_API_KEY}",
        "Content-Type": "application/json",
        "HTTP-Referer": "https://llm-judge-debias.example",
        "X-Title": "LLM Judge Debias"
    }

    judge_input = JUDGE_PROMPT.format(
        prompt=prompt[:8000],
        response_a=response_a[:8000],
        response_b=response_b[:8000]
    )

    payload = {
        "model": model_id,
        "messages": [{"role": "user", "content": judge_input}],
        "response_format": {"type": "json_object"},
        "temperature": 0
    }

    for attempt in range(max_retries):
        try:
            response = requests.post(OPENROUTER_URL, headers=headers, json=payload, timeout=60)
            response.raise_for_status()
            result = response.json()
            content = result['choices'][0]['message']['content'].strip()
            # Parse JSON response
            parsed = json.loads(content)
            winner = parsed.get('winner', '').upper()
            if winner in ['A', 'B']:
                return winner
            return None
        except json.JSONDecodeError:
            # Fallback: try to extract A or B from raw text
            if 'A' in content and 'B' not in content:
                return 'A'
            elif 'B' in content and 'A' not in content:
                return 'B'
            return None
        except Exception as e:
            if attempt < max_retries - 1:
                time.sleep(2 ** attempt)
            else:
                print(f"\nError calling {model_id}: {e}")
                return None
    return None


def load_json_results(path):
    """Load results from JSON file."""
    if path.exists():
        with open(path, 'r') as f:
            return json.load(f)
    return []


def save_json_results(results, path):
    """Save results to JSON file."""
    with open(path, 'w') as f:
        json.dump(results, f, indent=2)


def evaluate_dataset(df, model_name, model_id, output_path, checkpoint_path):
    """Evaluate all rows in a dataset with a specific judge model."""
    # Load checkpoint if exists
    results = load_json_results(checkpoint_path)
    done_ids = {r['id'] for r in results}

    if checkpoint_path.exists():
        print(f"  Resuming from checkpoint: {len(done_ids)} already done")

    remaining = df[~df['id'].isin(done_ids)]

    if len(remaining) == 0:
        print(f"  All {len(df)} samples already evaluated!")
        save_json_results(results, output_path)
        return results

    pbar = tqdm(remaining.iterrows(), total=len(remaining),
                desc=f"  {model_name}", unit="sample",
                bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]')

    for idx, row in pbar:
        prompt, response_a, response_b = get_prompt_and_responses(row)

        if not prompt or not response_a or not response_b:
            results.append({
                'id': row['id'],
                'model_a': row['model_a'],
                'model_b': row['model_b'],
                'human_winner': row['winner'],
                'judge_pick': None,
                'error': 'missing_content'
            })
        else:
            judge_pick = call_judge(prompt, response_a, response_b, model_id)
            results.append({
                'id': row['id'],
                'model_a': row['model_a'],
                'model_b': row['model_b'],
                'human_winner': row['winner'],
                'judge_pick': judge_pick,
            })

        # Save checkpoint every 10 samples
        if len(results) % 10 == 0:
            save_json_results(results, checkpoint_path)

        time.sleep(0.05)  # Rate limiting

    # Final save
    save_json_results(results, output_path)
    if checkpoint_path.exists():
        checkpoint_path.unlink()  # Remove checkpoint

    return results


def main():
    if not OPENROUTER_API_KEY:
        print("ERROR: OPENROUTER_API_KEY not set!")
        print("Run: source ~/.bash_profile")
        return

    data_dir = Path("data/arena_subsets")
    output_dir = Path("data/judge_results")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load all datasets
    datasets = {}
    total_samples = 0
    for f in sorted(data_dir.glob("*.parquet")):
        name = f.stem
        datasets[name] = pd.read_parquet(f)
        total_samples += len(datasets[name])
        print(f"Loaded {name}: {len(datasets[name])} rows")

    print(f"\nTotal: {total_samples} samples x 2 judges = {total_samples * 2} API calls")
    print(f"Estimated time: ~{total_samples * 2 * 4 / 60:.0f} minutes\n")

    # Evaluate with each judge model
    for model_name, model_id in MODELS.items():
        print(f"{'='*60}")
        print(f"Judge: {model_name} ({model_id})")
        print('='*60)

        for dataset_name, df in datasets.items():
            output_path = output_dir / f"{dataset_name}_{model_name}.json"
            checkpoint_path = output_dir / f"{dataset_name}_{model_name}.checkpoint.json"

            print(f"\n{dataset_name} ({len(df)} samples)")

            if output_path.exists() and not checkpoint_path.exists():
                existing = load_json_results(output_path)
                valid = [r for r in existing if r.get('judge_pick')]
                print(f"  Already complete: {len(valid)}/{len(existing)} valid")
                continue

            results = evaluate_dataset(df, model_name, model_id, output_path, checkpoint_path)
            valid = [r for r in results if r.get('judge_pick')]
            print(f"  Done: {len(valid)}/{len(results)} valid picks")

    print("\n" + "="*60)
    print("COMPLETE! Results saved to data/judge_results/")
    print("="*60)

    # Quick summary
    print("\nSummary:")
    for f in sorted(output_dir.glob("*.json")):
        if "checkpoint" in f.name:
            continue
        results = load_json_results(f)
        valid = [r for r in results if r.get('judge_pick')]
        print(f"  {f.name}: {len(valid)}/{len(results)} valid")


if __name__ == "__main__":
    main()
