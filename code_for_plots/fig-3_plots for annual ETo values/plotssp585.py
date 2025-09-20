import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# --- Configuration ---
excel_file = r'D:\CO2 Paper\Results_units_corrected\ssp585\annual_avg_eto_ssp585\ssp585_all_gcms_annual_avg.xlsx'
sheet_name = 'Sheet1'
output_filename = r'Fig4\avg_annual_eto_ssp585_large_font.png'

# --- Font and Style Configuration ---
# NOTE: A font size of 20 is very large and requires a larger figure size and layout adjustments.
main_font_size = 20

# Define the GCMs
gcms_in_plot = [
    'GFDL-ESM4', 'IPSL-CM6A-LR', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'UKESM1-0-LL'
]

# Define consistent colors and markers for GCMs
gcm_colors = {
    'GFDL-ESM4': '#2ca02c', 'IPSL-CM6A-LR': '#d62728', 'MPI-ESM1-2-HR': '#17becf',
    'MRI-ESM2-0': '#ff7f0e', 'UKESM1-0-LL': '#e377c2', 'Average (five GCMs)': 'black'
}
gcm_markers = {gcm: 'o' for gcm in gcms_in_plot}
gcm_markers['Average (five GCMs)'] = 'o'

# Define subplot configurations
subplot_configs = [
    {'co2': 'without_CO2', 'years': (2021, 2030), 'title': '(i)'},
    {'co2': 'with_CO2',    'years': (2021, 2030), 'title': '(ii)'},
    {'co2': 'without_CO2', 'years': (2051, 2060), 'title': '(iii)'},
    {'co2': 'with_CO2',    'years': (2051, 2060), 'title': '(iv)'},
    {'co2': 'without_CO2', 'years': (2091, 2100), 'title': '(v)'},
    {'co2': 'with_CO2',    'years': (2091, 2100), 'title': '(vi)'},
]

# Axis label
y_axis_label = r'Mean ET$_o$ (mm day$^{-1}$)'

# --- Data Loading and Processing ---
try:
    df = pd.read_excel(excel_file, sheet_name=sheet_name)
except FileNotFoundError:
    print(f"Error: File not found at '{excel_file}'"); exit()
except Exception as e:
    print(f"Error reading Excel file: {e}"); exit()

required_cols = ['GCM', 'CO2', 'Year', 'Mean Ref. ETo']
if not all(col in df.columns for col in required_cols):
    print(f"Error: Missing required columns. Found: {list(df.columns)}. Required: {required_cols}"); exit()

df['Year'] = df['Year'].astype(int)

# --- Plotting ---
# <<< CHANGE 1: Increased figsize and switched to constrained_layout=True >>>
fig, axes = plt.subplots(3, 2, figsize=(15, 18), sharey=True, constrained_layout=True)
axes = axes.flatten()

plot_handles, plot_labels = [], []

# Loop through each subplot configuration
for i, config in enumerate(subplot_configs):
    ax = axes[i]
    co2_scenario, (start_year, end_year) = config['co2'], config['years']

    df_sub = df[(df['CO2'] == co2_scenario) & (df['Year'].between(start_year, end_year))].copy()

    if df_sub.empty:
        # <<< CHANGE 2: Using main_font_size variable >>>
        ax.set_title(config['title'], loc='left', fontsize=main_font_size, weight='bold')
        ax.text(0.5, 0.5, 'No Data', ha='center', va='center', transform=ax.transAxes, fontsize=main_font_size)
        continue

    # Calculate average
    df_sub_avg = df_sub[df_sub['GCM'].isin(gcms_in_plot)]
    average_eto = df_sub_avg.groupby('Year')['Mean Ref. ETo'].mean().reset_index()

    # Plot individual GCM lines
    for gcm in gcms_in_plot:
        gcm_data = df_sub[df_sub['GCM'] == gcm]
        if not gcm_data.empty:
            line, = ax.plot(gcm_data['Year'], gcm_data['Mean Ref. ETo'], marker=gcm_markers.get(gcm, 'o'),
                            linestyle='-', color=gcm_colors.get(gcm, 'gray'), label=gcm,
                            linewidth=2.0, markersize=7) # Slightly bigger lines/markers
            if gcm not in plot_labels:
                plot_handles.append(line); plot_labels.append(gcm)

    # Plot the average line
    if not average_eto.empty:
        avg_label = 'Average (five GCMs)'
        line, = ax.plot(average_eto['Year'], average_eto['Mean Ref. ETo'], marker=gcm_markers.get(avg_label, 'o'),
                        linestyle='--', color=gcm_colors.get(avg_label, 'black'), label=avg_label,
                        linewidth=2.2, markersize=7) # Slightly bigger lines/markers
        if avg_label not in plot_labels:
            plot_handles.append(line); plot_labels.append(avg_label)

    # --- Subplot Customization ---
    # <<< CHANGE 3: Using main_font_size for all text elements >>>
    ax.set_title(config['title'], loc='left', fontsize=12, weight='bold')
    ax.grid(True, linestyle=':', linewidth=0.7)
    ax.set_xlim(start_year - 0.5, end_year + 0.5)
    ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True, nbins=len(range(start_year, end_year + 1))))
    
    # Set tick label sizes for both axes at once
    ax.tick_params(axis='both', labelsize=14)

    ax.set_ylim(3.0, 6.0) # Adjusted y-limits for better visibility

    # Add x-axis labels only to bottom row plots
    if i >= 4: # Bottom row
        ax.set_xlabel('Year', fontsize=main_font_size)

# --- Final Figure Adjustments ---

# Add ONE centered Y-axis label for the whole figure
# <<< CHANGE 4: Using main_font_size for the shared Y-axis label >>>
try:
    fig.supylabel(y_axis_label, fontsize=main_font_size)
except AttributeError:
    fig.text(0.02, 0.5, y_axis_label, va='center', rotation='vertical', fontsize=main_font_size)

# Create a single legend below the subplots
ordered_handles, ordered_labels = [], []
avg_label = 'Average (five GCMs)'
if avg_label in plot_labels:
    avg_index = plot_labels.index(avg_label)
    ordered_handles.append(plot_handles.pop(avg_index)); ordered_labels.append(plot_labels.pop(avg_index))
for gcm in gcms_in_plot:
    if gcm in plot_labels:
        gcm_index = plot_labels.index(gcm)
        ordered_handles.append(plot_handles.pop(gcm_index)); ordered_labels.append(plot_labels.pop(gcm_index))
ordered_handles.extend(plot_handles); ordered_labels.extend(plot_labels)

# <<< CHANGE 5: Using main_font_size and adjusting position for the legend >>>
fig.legend(ordered_handles, ordered_labels,
           loc='lower center',
           bbox_to_anchor=(0.5, -0.05), # Moved further down to avoid overlap
           ncol=3, # Fewer columns may look better with large font
           fontsize=main_font_size,
           frameon=False)

# <<< CHANGE 6: Removed plt.tight_layout(), as constrained_layout is now active >>>

# --- Save and Show ---
fig.text(0.03, 0.99, '(b)', ha='center', va='top', fontsize=main_font_size, fontweight='bold')
fig.suptitle('Scenario : SSP5-8.5', fontsize=main_font_size, fontweight='bold')
plt.savefig(output_filename, dpi=300, bbox_inches='tight')
print(f"Plot saved as {output_filename}")
plt.show()

# r""