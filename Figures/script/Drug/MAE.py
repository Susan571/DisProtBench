import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors

# File paths and dataset labels
files = [
    ('Full', '/home/susan2000/Benchmark-final/Drug/full/test/Test_results.csv'),
    ('Disorder_30', '/home/susan2000/Benchmark-final/Drug/disorder_30/test/Test_results.csv'),
    ('Disorder_50', '/home/susan2000/Benchmark-final/Drug/disorder_50/test/Test_results.csv')
]

fig, axes = plt.subplots(1, 3, figsize=(28, 7), sharex=True)

for idx, (label, path) in enumerate(files):
    df = pd.read_csv(path)
    models = df['Model']
    mae = df['MAE']
    r = df['R']
    x = np.arange(len(models))

    ax1 = axes[idx]
    # Define start and end colors
    start_color = '#335372'
    end_color = '#CEDCE9'

    # Convert hex to RGB
    start_rgb = np.array(mcolors.to_rgb(start_color))
    end_rgb = np.array(mcolors.to_rgb(end_color))

    # Interpolate colors
    n_models = len(models)
    bar_colors = [mcolors.to_hex(start_rgb + (end_rgb - start_rgb) * i / (n_models - 1)) for i in range(n_models)]

    bars = ax1.bar(x, mae, color=bar_colors, width=0.6)
    ax1.set_ylim(0, max(mae)*1.3)
    # Only leftmost subplot gets MAE y-label
    if idx == 0:
        ax1.set_ylabel('MAE', fontsize=20)
        ax1.tick_params(axis='y', labelsize=20)
    else:
        ax1.set_ylabel('')
        ax1.tick_params(axis='y', labelsize=0)
    ax1.set_xticks(x)
    # Remove x tick labels
    ax1.set_xticklabels([''] * len(models))
    #ax1.set_xlabel('Model', fontsize=20)
    # Remove subplot title, will add label below
    # ax1.set_title(label, fontsize=20)
    #labelcolor='#4C72B0'

    # Add MAE value labels
    for bar in bars:
        height = bar.get_height()
        ax1.annotate(f'{height:.3f}', xy=(bar.get_x() + bar.get_width() / 2, height),
                     xytext=(0, 8), textcoords="offset points",
                     ha='center', va='bottom', fontsize=15, color='#4C72B0')

    # Add legend handles for later
    if idx == 0:
        legend_handles = [plt.Line2D([0], [0], color=bar_colors[i], lw=8, label=models.iloc[i]) for i in range(len(models))]

    # Add dataset label under each subplot
    if label == 'Disorder_30':
        display_label = r'plDDT $\geq$ 30'
    elif label == 'Disorder_50':
        display_label = r'plDDT $\geq$ 50'
    else:
        display_label = label
    ax1.text(0.5, -0.06, display_label, ha='center', va='center', fontsize=20, transform=ax1.transAxes)

    # Second y-axis for R
    ax2 = ax1.twinx()
    ax2.plot(x, r, color='#2e7ebb', marker='o', linewidth=2, label='R')
    ax2.set_ylim(0, 1.1)
    if idx == len(files) - 1:
        ax2.set_ylabel('R', fontsize=20, color='#2e7ebb', labelpad=20)
        ax2.tick_params(axis='y', labelcolor='#2e7ebb', labelsize=20)
    else:
        ax2.set_ylabel('')
        ax2.tick_params(axis='y', labelcolor='#2e7ebb', labelsize=0)

    # Add R value labels
    for xi, ri in zip(x, r):
        ax2.annotate(f'{ri:.3f}', xy=(xi, ri), xytext=(0, 8), textcoords="offset points",
                     ha='center', va='bottom', fontsize=15, color='#2e7ebb')

    # No legends for R

# Remove the global x label
# fig.text(0.5, 0.01, 'Model', ha='center', va='center', fontsize=20)

# Add a two-row legend for all subplots at the bottom
num_first_row = 6
num_second_row = 5
handles_row1 = legend_handles[:num_first_row]
labels_row1 = [h.get_label() for h in handles_row1]
handles_row2 = legend_handles[num_first_row:num_first_row+num_second_row]
labels_row2 = [h.get_label() for h in handles_row2]

legend1 = fig.legend(handles=handles_row1, labels=labels_row1, fontsize=20, loc='lower right', bbox_to_anchor=(0.98, -0.03), ncol=num_first_row, frameon=False)
legend2 = fig.legend(handles=handles_row2, labels=labels_row2, fontsize=20, loc='lower right', bbox_to_anchor=(0.98, -0.08), ncol=num_second_row, frameon=False)
fig.add_artist(legend1)
fig.add_artist(legend2)

plt.subplots_adjust(left=0.15, right=0.95, wspace=0.03, bottom=0.13)
plt.savefig('/home/susan2000/Benchmark-final/Drug/MAE_R_overall.png', bbox_inches='tight')
plt.close()
print("Combined MAE and R dual y-axis plots saved as MAE_R_overall.png")
