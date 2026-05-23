import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv('/Users/kallolnaha/Documents/mix_projects/mindrouter/PROBE/semantic_similarity_results.csv')

df['namespace'] = df['namespace'].map({
    'molecular_function': 'MF',
    'biological_process': 'BP',
    'cellular_component': 'CC'
})

summary = df.groupby(['model', 'namespace'])[['llm_sim']].mean().reset_index()

# Correctly computed random baselines using proper GO ancestor databases
random_baselines = {
    'MF': 0.145,
    'BP': 0.133,
    'CC': 0.258
}

print("Random baselines being used:")
print(random_baselines)

models = ['Mistral Large 123B', 'Llama 3.3 70B', 'Qwen2.5 72B']
namespaces = ['MF', 'BP', 'CC']
color_llm = '#2166AC'
color_random = '#BDBDBD'

fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharey=False)

for idx, ns in enumerate(namespaces):
    ax = axes[idx]
    ns_data = summary[summary['namespace'] == ns]
    x = np.arange(len(models))
    width = 0.35
    llm_vals = []
    rand_val = random_baselines[ns]
    print(f"NS: {ns}, Random: {rand_val}")
    for m in models:
        row = ns_data[ns_data['model'] == m]
        llm_vals.append(row['llm_sim'].values[0] if len(row) > 0 else 0)

    bars1 = ax.bar(x - width/2, llm_vals, width,
                   label='LLM Prediction', color=color_llm,
                   alpha=0.85, edgecolor='black', linewidth=0.5)
    bars2 = ax.bar(x + width/2, [rand_val]*3, width,
                   label='Random Baseline', color=color_random,
                   alpha=0.85, edgecolor='black', linewidth=0.5)

    for bar in bars1:
        ax.text(bar.get_x() + bar.get_width()/2.,
                bar.get_height() + 0.01,
                f'{bar.get_height():.3f}',
                ha='center', va='bottom', fontsize=8, fontweight='bold')
    for bar in bars2:
        ax.text(bar.get_x() + bar.get_width()/2.,
                bar.get_height() + 0.01,
                f'{bar.get_height():.3f}',
                ha='center', va='bottom', fontsize=8)

    ax.set_title(f'GO Namespace: {ns}', fontsize=12, fontweight='bold', pad=10)
    ax.set_xticks(x)
    ax.set_xticklabels(['Mistral\nLarge', 'Llama\n3.3 70B', 'Qwen2.5\n72B'], fontsize=9)
    ax.set_ylim(0, 0.85)
    ax.set_ylabel('Mean Semantic Similarity (Wang/BMA)', fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(axis='y', alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('/Users/kallolnaha/Documents/mix_projects/mindrouter/PROBE/fig8_semantic_similarity.pdf',
            bbox_inches='tight', dpi=300)
plt.savefig('/Users/kallolnaha/Documents/mix_projects/mindrouter/PROBE/fig8_semantic_similarity.png',
            bbox_inches='tight', dpi=300)
print('Done!')