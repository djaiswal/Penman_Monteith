import matplotlib.pyplot as plt
import numpy as np

# --- 1. Data Setup ---
# Define the main categories for the x-axis
periods = ['Near-Term (2021-2030)', 'Mid-Term (2051-2060)', 'Long-Term (2091-100)']

# Store all the data in a dictionary for easy access.
# I've included the data from your SSP585 image and created plausible data for SSP126.
plot_data = {
    'SSP1-2.6': {
        'without CO2': [4.339826, 4.51581, 4.540816],
        'with CO2':    [3.801329, 3.889381, 3.948705]
    },
    'SSP5-8.5': {
        'without CO2': [4.224385, 4.511875, 4.950793],
        'with CO2':    [3.78733, 3.701205, 3.532085]
    }
}

# --- 2. Plotting Configuration ---
bar_width = 0.35  # The width of the bars
#change to pastel color
bar_color_without_co2 = "#c78746"  # Pastel orange
bar_color_with_co2 = '#a3c6f1'  # Pastel blue
panel_titles = ['(a) SSP1-2.6', '(b) SSP5-8.5'] # Titles for each subplot
y_axis_label = r"(mm day$^{-1}$)" # Example Y-axis Label
y_axis_limits = [3.0, 5.2]

# --- 3. Create the Figure and Subplots ---
# Create a figure with 1 row and 2 columns of subplots.
# `sharey=True` makes both plots share the same y-axis scale.
fig, axes = plt.subplots(1, 2, figsize=(15, 7), sharey=True)

# Flatten axes array for easy iteration, in case you want more rows later
axes = axes.flatten()

# --- 4. Loop Through Each Subplot to Draw the Bars ---
ssp_scenarios = list(plot_data.keys())

for i, ax in enumerate(axes):
    ssp_name = ssp_scenarios[i]
    data = plot_data[ssp_name]
    without_co2_vals = data['without CO2']
    with_co2_vals = data['with CO2']

    # Set up the x-axis positions for the bars
    x = np.arange(len(periods))  # the label locations [0, 1, 2]
    
    # Draw the bars for "without CO2"
    rects1 = ax.bar(x - bar_width/2, without_co2_vals, bar_width, 
                    label=r'$ET_o^{original}$', color=bar_color_without_co2, 
                     edgecolor='black')

    # Draw the bars for "with CO2"
    rects2 = ax.bar(x + bar_width/2, with_co2_vals, bar_width, 
                    label=r"$ET_o^{CO_2}$", color=bar_color_with_co2,
                     edgecolor='black')

    # --- 5. Add Labels, Titles, and Ticks ---
    ax.set_title(panel_titles[i], fontsize=16, fontweight='bold', pad=15)
    
    # Set the x-ticks to be in the middle of the two bars for each group
    ax.set_xticks(x)
    ax.set_xticklabels(periods, fontsize=12)
    
    # Set y-axis properties (only show label on the first plot)
    if i == 0:
        ax.set_ylabel(y_axis_label, fontsize=14)
    ax.tick_params(axis='y', labelsize=12)
    
    # Add a clean gridline
    ax.yaxis.grid(True, linestyle='--', which='major', color='grey', alpha=.25)
    ax.set_axisbelow(True) # Send gridlines to the back

    # Set Y-axis limits
    ax.set_ylim(y_axis_limits)

    # Add data labels on top of each bar
    ax.bar_label(rects1, padding=3, fmt='%.2f', fontsize=11)
    ax.bar_label(rects2, padding=3, fmt='%.2f', fontsize=11)

# --- 6. Final Touches ---
# Add a single, shared legend for the entire figure
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.95), ncol=2, fontsize=12, frameon=False)

# Adjust layout to prevent labels from overlapping
plt.tight_layout(rect=[0, 0, 1, 0.9]) # Adjust rect to make space for the figure legend

# Save and show the plot
plt.savefig("ssp_comparison_barchartv6.png", dpi=300)
plt.show()