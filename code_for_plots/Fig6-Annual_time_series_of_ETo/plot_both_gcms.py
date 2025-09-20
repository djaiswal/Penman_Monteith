import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D # Needed for creating a custom legend

# --- Overall Configuration ---
# --- DEFINE THE CONFIGURATIONS FOR EACH SCENARIO PANEL ---
scenario_configs = {
    'ssp126': {
        # !!! IMPORTANT: VERIFY THIS IS THE CORRECT FILE FOR SSP1-2.6 !!!
        'excel_file': r'D:\CO2 Paper\Results_units_corrected\ssp126\annual_avg_ssp126\combined_annual_avg.xlsx',
        'sheet_name': 'Sheet1',
        'panel_title': 'SSP1-2.6',
        'index': 'a'  # Used for subplot titles
    },
    'ssp585': {
        'excel_file': r"D:\CO2 Paper\Results_units_corrected\ssp585\annual_avg_eto_ssp585\ssp585_all_gcms_annual_avg.xlsx",
        'sheet_name': 'Sheet1',
        'panel_title': 'SSP5-8.5',
        'index': 'b'  # Used for subplot titles
    }
}
# Name for the saved plot file
output_filename = r'Fig4\ssp126_vs_ssp585_timeseries_comparison_font_v7.png'

# --- Shared Plotting Styles & Definitions ---
gcms_in_plot = [
    'GFDL-ESM4', 'IPSL-CM6A-LR', 'MPI-ESM1-2-HR',
    'MRI-ESM2-0', 'UKESM1-0-LL'
]
gcm_colors = {
    'GFDL-ESM4': '#2ca02c', 'IPSL-CM6A-LR': '#d62728',
    'MPI-ESM1-2-HR': '#17becf', 'MRI-ESM2-0': '#ff7f0e',
    'UKESM1-0-LL': '#e377c2', 'Average (five GCMs)': 'black'
}
gcm_markers = {gcm: 'o' for gcm in gcms_in_plot}
gcm_markers['Average (five GCMs)'] = 'o'
gcm_linestyles = {gcm: '-' for gcm in gcms_in_plot}
gcm_linestyles['Average (five GCMs)'] = '--'

# Define the 3x2 grid layout once
subplot_configs = [
    {'co2': 'without_CO2', 'years': (2021, 2030), 'title': '(i)'},
    {'co2': 'with_CO2',    'years': (2021, 2030), 'title': '(ii)'},
    {'co2': 'without_CO2', 'years': (2051, 2060), 'title': '(iii)'},
    {'co2': 'with_CO2',    'years': (2051, 2060), 'title': '(iv)'},
    {'co2': 'without_CO2', 'years': (2091, 2100), 'title': '(v)'},
    {'co2': 'with_CO2',    'years': (2091, 2100), 'title': '(vi)'},
]

y_axis_label = r'Mean ET$_o$ (mm day$^{-1}$)'

#=============================================================================
# REUSABLE PLOTTING FUNCTION
# This function draws the 3x2 grid for one scenario.
#=============================================================================
def plot_scenario_panel(subfig, config):
    """
    Reads data and draws a 3x2 grid of time-series plots on a given subfigure.
    """
    panel_title = config['panel_title']
    excel_file = config['excel_file']
    sheet_name = config['sheet_name']
    
    subfig.text(0.005, 0.99, f'({config["index"]})', fontsize=20, fontweight='bold', va='top', ha='left')
    subfig.suptitle(panel_title, fontsize=20, fontweight='bold')
    # subfig.supylabel(y_axis_label, fontsize=20)

    # --- Load Data for this Panel ---
    try:
        df = pd.read_excel(excel_file, sheet_name=sheet_name)
        df['Year'] = df['Year'].astype(int)
    except Exception as e:
        print(f"FATAL: Could not load data for {panel_title} from {excel_file}. Error: {e}")
        # Draw error message on the panel
        ax = subfig.subplots()
        ax.text(0.5, 0.5, f'Error loading data for\n{panel_title}', ha='center', va='center', color='red')
        ax.set_axis_off()
        return

    # Create the 3x2 axes grid within the subfigure
    axes = subfig.subplots(3, 2, sharey=True)
    axes = axes.flatten()

    # Loop through each of the 6 subplots
    for i, plot_conf in enumerate(subplot_configs):
        ax = axes[i]
        co2_scenario = plot_conf['co2']
        start_year, end_year = plot_conf['years']

        # Filter data for the current small plot
        df_sub = df[(df['CO2'] == co2_scenario) & (df['Year'].between(start_year, end_year))].copy()

        if df_sub.empty:
            print(f"Warning: No data for {panel_title}, {co2_scenario}, {start_year}-{end_year}")
            ax.text(0.5, 0.5, 'No Data', ha='center', va='center', transform=ax.transAxes)
        else:
            # Calculate average
            df_sub_avg = df_sub[df_sub['GCM'].isin(gcms_in_plot)]
            average_eto = df_sub_avg.groupby('Year')['Mean Ref. ETo'].mean().reset_index()

            # Plot individual GCM lines
            for gcm in gcms_in_plot:
                gcm_data = df_sub[df_sub['GCM'] == gcm]
                if not gcm_data.empty:
                    ax.plot(gcm_data['Year'], gcm_data['Mean Ref. ETo'],
                            marker=gcm_markers[gcm], linestyle=gcm_linestyles[gcm],
                            color=gcm_colors[gcm], linewidth=1.5, markersize=5)
            
            # Plot average line
            if not average_eto.empty:
                ax.plot(average_eto['Year'], average_eto['Mean Ref. ETo'],
                        marker=gcm_markers['Average (five GCMs)'],
                        linestyle=gcm_linestyles['Average (five GCMs)'],
                        color=gcm_colors['Average (five GCMs)'],
                        linewidth=1.8, markersize=5)

        # --- Subplot Customization ---
        ax.set_title(plot_conf['title'], loc='left', fontsize=20, weight='bold')
        ax.grid(True, linestyle=':', linewidth=0.5)
        ax.set_xlim(start_year - 0.5, end_year + 0.5)
        ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True, nbins=11))
        # Make all x-axis ticks visible (both bottom and top, major and minor)
        # ax.tick_params(axis='x', which='both', direction='in', length=6, width=1.5, colors='black',
        #            bottom=True, top=True, labelbottom=True)
        ax.tick_params(axis='x', labelsize=12, rotation=45)
        ax.tick_params(axis='y', labelsize=18)
        ax.set_ylim(3.0, 6.0) # Consistent Y-axis limits
        if i >= 4: # Bottom row plots
            ax.set_xlabel('Year', fontsize=20)
        

#=============================================================================
# MAIN SCRIPT EXECUTION
#=============================================================================

# 1. Create the main figure and two subfigures (panels)
fig = plt.figure(figsize=(20, 17), constrained_layout=True)
subfigs = fig.subfigures(nrows=1, ncols=2, wspace=0.07)

# 2. Add ONE shared Y-axis label for the whole figure


# 3. Call the plotting function for each panel
print("--- Plotting SSP1-2.6 Panel ---")
plot_scenario_panel(subfigs[0], scenario_configs['ssp126'])

print("\n--- Plotting SSP5-8.5 Panel ---")
plot_scenario_panel(subfigs[1], scenario_configs['ssp585'])

# 4. Create a single, shared legend for the whole figure
legend_handles = []
avg_label = 'Average (five GCMs)'

# Add the handle for the Average line first
legend_handles.append(Line2D([0], [0], color=gcm_colors[avg_label], lw=2,
                             linestyle=gcm_linestyles[avg_label],
                             marker=gcm_markers[avg_label], label=avg_label))
# Add handles for each GCM
for gcm in gcms_in_plot:
    legend_handles.append(Line2D([0], [0], color=gcm_colors[gcm], lw=2,
                                 linestyle=gcm_linestyles[gcm],
                                 marker=gcm_markers[gcm], label=gcm))

fig.legend(handles=legend_handles,
           loc='lower center',      # Position below plots
           bbox_to_anchor=(0.5, -0.03), # Fine-tune position (y<0 is below axes)
           ncol=6,                  # Arrange in 6 columns
           fontsize=20,
           frameon=False)

# 5. Save and show the final combined plot
# A final adjustment for the suptitle and legend might be needed
fig.supylabel(y_axis_label, fontsize=20)
fig.set_constrained_layout_pads(w_pad=0.04, h_pad=0.04, wspace=0.02, hspace=0.02)
plt.savefig(output_filename, dpi=300, bbox_inches='tight')
print(f"\nPlot saved as {output_filename}")
plt.show()