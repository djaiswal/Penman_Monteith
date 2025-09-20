import matplotlib.pyplot as plt
import matplotlib as mpl
import xarray as xr
import numpy as np
import os
import geopandas as gpd
import fiona

# --- Overall Configuration ---
output_filename = r'D:\CO2 Paper\plot_codes\Fig3\ssp126_vs_ssp585_comparisonV10.png'
output_dpi = 300

# --- Shared Plotting Configuration ---
lon_name = 'lon'
lat_name = 'lat'
variable_units = "mm day" + r'$^{-1}$'
data_threshold = 200.0
num_contour_levels = 6
row_titles = ['2021-2030', '2051-2060', '2091-2100']
column_titles = ["GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL"]
gcms = ["GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL"]
periods = ["2021-2030", "2051-2060", "2091-2100"]
row_cmaps = ["Oranges", "Greens", "Blues"]
n_rows, n_cols = len(row_titles), len(column_titles)

# --- SSP-Specific Configurations ---
ssp_configs = {
    'SSP1-2.6': {
        'base_dir': r"D:\CO2 Paper\Results_units_corrected\ssp126\decadal_average",
        'panel_title': 'SSP1-2.6',
        'variable_name': 'decadal_avg'
    },
    'SSP5-8.5': {
        'base_dir': r"D:\CO2 Paper\Results_units_corrected\ssp585\decadal_average",
        'panel_title': 'SSP5-8.5',
        'variable_name': "decadal_avg"
    }
}

# <<< MODIFIED: Manual Colorbar Limits >>>
# Define the min and max for each row of each panel.
# The structure is {panel_title: [(row0_min, row0_max), (row1_min, row1_max), ...]}
manual_clims = {
    'SSP1-2.6': [
       (0, 1.5),  # Row 0: 2021-2030
        (0, 2),  # Row 1: 2051-2060
        (0, 3.5)  
    ],
    'SSP5-8.5': [
        (0, 1.5),  # Row 0: 2021-2030
        (0, 2),  # Row 1: 2051-2060
        (0, 3.5)   # Row 2: 2091-2100
    ]
}


# --- Shapefile Configuration ---
india_shapefile_path = r"D:\CO2 Paper\Latest_frontiers_in_water\Shapefiles\india_state_boundary_projected.shp"
shapefile_outline_color = 'black'
attribute_column_name = 'State_Name'
island_names_to_exclude = ['Lakshadweep', 'Andaman & Nicobar']

#=============================================================================
# REUSABLE PLOTTING FUNCTION
#=============================================================================
# <<< MODIFIED: Function signature updated to accept manual limits >>>
def plot_ssp_panel(subfig, ssp_config, manual_clims_for_panel, india_map, index, show_row_titles, show_colorbars):
    """
    Creates a 3x5 grid of plots for a single SSP scenario on a given subfigure.
    
    Args:
        ...
        manual_clims_for_panel (list): A list of (min, max) tuples for the colorbars.
        ...
    """
    base_netcdf_dir = ssp_config['base_dir']
    panel_title = ssp_config['panel_title']
    
    if panel_title == 'SSP1-2.6':
        subfig.text(0.05, 0.45, "Latitude (°N)", fontsize=18 , rotation=90, va='center')
    subfig.suptitle(panel_title, fontsize=16, fontweight='bold')
    
    # 1. Construct file paths for this SSP
    netcdf_files = []
    for period in periods:
        if panel_title == 'SSP1-2.6':
            row_files = [os.path.join(base_netcdf_dir, f"difference_average_{gcm}_{period}_daily_difference.nc") for gcm in gcms]
        elif panel_title == 'SSP5-8.5':
            row_files = [os.path.join(base_netcdf_dir, f"difference_average_{gcm}_{period}_daily_difference.nc") for gcm in gcms]
        else:
            print(f"Warning: Unrecognized SSP panel title '{panel_title}'. Exiting.")
            exit()
        netcdf_files.append(row_files)

    # <<< MODIFIED: Use manual color limits instead of calculating them >>>
    print(f"\n--- Using manual color limits for {panel_title} ---")
    row_min = [lim[0] for lim in manual_clims_for_panel]
    row_max = [lim[1] for lim in manual_clims_for_panel]
    for i, (vmin, vmax) in enumerate(manual_clims_for_panel):
        print(f"Row {i} ('{row_titles[i]}') using manual range: min={vmin:.2f}, max={vmax:.2f}")


    # 3. Create axes grid within the subfigure
    axes = subfig.subplots(nrows=n_rows, ncols=n_cols, sharex=True, sharey=True)
    row_mappables = [None] * n_rows
    cbar_label_fontsize = 16
    cbar_tick_fontsize = 16
    axis_tick_fontsize = 16

    # 4. Loop through files and plot on each axis
    print(f"--- Plotting grid for {panel_title} ---")
    for i in range(n_rows):
        for j in range(n_cols):
            ax = axes[i, j]
            nc_file = netcdf_files[i][j]
            if not os.path.exists(nc_file):
                ax.text(0.5, 0.5, 'File Not Found', ha='center', va='center', transform=ax.transAxes, color='red')
                ax.set_facecolor('lightgrey')
                continue

            try:
                with xr.open_dataset(nc_file) as ds:
                    if ssp_config['variable_name'] not in ds.variables: 
                        ax.text(0.5, 0.5, 'Var Not Found', ha='center', va='center', transform=ax.transAxes, color='orange')
                        ax.set_facecolor('lightgrey'); continue

                    data = ds[ssp_config['variable_name']].isel(time=0) if 'time' in ds[ssp_config['variable_name']].dims else ds[ssp_config['variable_name']]
                    lons, lats = ds[lon_name].values, ds[lat_name].values
                    if lons.ndim == 1: lons, lats = np.meshgrid(lons, lats)
                    
                    # NOTE: data_threshold is less relevant with manual limits, but can still be used for masking extreme outliers if needed
                    data_filtered = data.where(data <= data_threshold).values
                    
                    levels = np.linspace(row_min[i], row_max[i], num_contour_levels + 1)
                    norm = mpl.colors.Normalize(vmin=row_min[i], vmax=row_max[i])
                    
                    contour_plot = ax.contourf(lons, lats, data_filtered, levels=levels, cmap=row_cmaps[i], norm=norm, extend='both')
                    row_mappables[i] = contour_plot
                    india_map.plot(ax=ax, facecolor='none', edgecolor=shapefile_outline_color, linewidth=0.3)
                    ax.tick_params(axis='both', labelsize=axis_tick_fontsize)
            except Exception as e:
                print(f"ERROR plotting {os.path.basename(nc_file)}: {e}")
                ax.set_facecolor('salmon')

            # Ticks, labels, and titles
            if show_row_titles and j == 0:
                ax.text(-0.50, 0.5, row_titles[i], transform=ax.transAxes, ha='right', va='center', rotation=90, fontsize=20)
            
            if i == 0:
                ax.set_title(column_titles[j], fontsize=14, pad=10)
            
            ax.tick_params(axis='x', labelbottom=(i == n_rows - 1))
            ax.tick_params(axis='y', labelleft=(j == 0))
            ax.grid(True, linestyle='--', alpha=0.5)

    # 5. Add colorbars for the panel
    if show_colorbars:
        for row_idx, mappable in enumerate(row_mappables):
            if mappable:
                cbar_ticks = np.linspace(row_min[row_idx], row_max[row_idx], 6)
                cbar = subfig.colorbar(mappable, ax=axes[row_idx, :].tolist(), location='right',
                                       orientation='vertical', shrink=0.8, aspect=20, pad=0.03, ticks=cbar_ticks)
                
                if row_idx == n_rows // 2:
                    cbar.set_label(f'({variable_units})', size=cbar_label_fontsize, labelpad=15)
                
                cbar.ax.tick_params(labelsize=cbar_tick_fontsize)
                
    subfig.text(0.02, 0.98, f"({index})", fontsize=20, fontweight='bold', va='top', ha='left')


#=============================================================================
# MAIN SCRIPT EXECUTION
#=============================================================================

# 1. Load shapefile once
print("Loading and preparing shapefile...")
with fiona.Env(SHAPE_RESTORE_SHX='YES'):
    india_map = gpd.read_file(india_shapefile_path)
india_map = india_map.to_crs(epsg=4326)
india_map = india_map[~india_map[attribute_column_name].isin(island_names_to_exclude)]
print("Shapefile ready.")

# 2. Create the main figure and subfigures
fig = plt.figure(figsize=(24, 8), constrained_layout=True)
subfigs = fig.subfigures(nrows=1, ncols=2, wspace=0.05)

# 3. Call the plotting function for each subfigure/SSP
ssp_names = list(ssp_configs.keys())

# <<< MODIFIED: Pass the manual color limits for the specific panel to the function >>>
ssp1_name = ssp_names[0]
plot_ssp_panel(subfigs[0], ssp_configs[ssp1_name], manual_clims[ssp1_name], india_map, index='a',
               show_row_titles=True, show_colorbars=False)

ssp2_name = ssp_names[1]
plot_ssp_panel(subfigs[1], ssp_configs[ssp2_name], manual_clims[ssp2_name], india_map, index='b',
               show_row_titles=False, show_colorbars=True)

# Add single, shared labels for the entire figure
fig.supxlabel("Longitude (°E)", fontsize=22)
# fig.supylabel("Latitude (°N)", fontsize=16)


# 4. Save the final figure
print(f"\nSaving combined figure to {output_filename}...")
plt.savefig(output_filename, dpi=output_dpi, bbox_inches='tight')
print("Figure saved successfully.")

# 5. Show plot
plt.show()