import matplotlib.pyplot as plt
import matplotlib as mpl
import xarray as xr
import numpy as np
import os
import geopandas as gpd
import fiona

# --- Overall Configuration ---
output_filename = r'D:\CO2 Paper\plot_codes\Fig3\ssp126_vs_ssp585_comparison_final_v3.png'
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
row_cmaps = ["Oranges", "Greens", "Blues"]  # One colormap for each row
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

# --- Shapefile Configuration ---
india_shapefile_path = r"D:\CO2 Paper\Latest_frontiers_in_water\Shapefiles\india_state_boundary_projected.shp"
shapefile_outline_color = 'black'
attribute_column_name = 'State_Name'
island_names_to_exclude = ['Lakshadweep', 'Andaman & Nicobar']

#=============================================================================
# REUSABLE PLOTTING FUNCTION
# This function will draw one entire 3x5 panel for a given SSP.
#=============================================================================
def plot_ssp_panel(subfig, ssp_config, india_map, index, show_row_titles, show_colorbars):
    """
    Creates a 3x5 grid of plots for a single SSP scenario on a given subfigure.
    
    Args:
        subfig (matplotlib.figure.SubFigure): The subfigure object to draw on.
        ssp_config (dict): A dictionary containing 'base_dir' and 'panel_title'.
        india_map (GeoDataFrame): The loaded and processed GeoPandas shapefile.
        index (str): The index label for the panel (e.g., 'a', 'b').
        show_row_titles (bool): If True, display titles for each row.
        show_colorbars (bool): If True, display colorbars for each row.
    """
    base_netcdf_dir = ssp_config['base_dir']
    panel_title = ssp_config['panel_title']
    
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

    # 2. Calculate row-specific min/max values for this SSP's data
    print(f"\n--- Calculating min/max for {panel_title} ---")
    row_min = [np.inf] * n_rows
    row_max = [-np.inf] * n_rows
    for row_idx in range(n_rows):
        found_data_in_row = False
        for nc_file in netcdf_files[row_idx]:
            if not os.path.exists(nc_file):
                print(f"Warning: File not found, skipping for min/max: {nc_file}")
                continue
            try:
                with xr.open_dataset(nc_file) as ds:
                    if ssp_config['variable_name'] not in ds.variables: continue
                    data = ds[ssp_config['variable_name']].isel(time=0) if 'time' in ds[ssp_config['variable_name']].dims else ds[ssp_config['variable_name']]
                    data_filtered = data.where(data <= data_threshold)
                    current_min, current_max = data_filtered.min().item(), data_filtered.max().item()
                    if not np.isnan(current_min): row_min[row_idx] = min(row_min[row_idx], current_min); found_data_in_row = True
                    if not np.isnan(current_max): row_max[row_idx] = max(row_max[row_idx], current_max)
            except Exception as e:
                print(f"Warning: Could not read {os.path.basename(nc_file)} for min/max: {e}")
        if not found_data_in_row: row_min[row_idx], row_max[row_idx] = 0, 1
        if row_min[row_idx] == row_max[row_idx]: row_max[row_idx] += 1e-6
        print(f"Row {row_idx} ('{row_titles[row_idx]}') range: min={row_min[row_idx]:.2f}, max={row_max[row_idx]:.2f}")

    # 3. Create axes grid within the subfigure
    axes = subfig.subplots(nrows=n_rows, ncols=n_cols, sharex=True, sharey=True)
    row_mappables = [None] * n_rows
    cbar_label_fontsize = 16
    cbar_tick_fontsize = 16
    axis_tick_fontsize = 16 # Adjusted for better fit

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
            # <<< MODIFIED: Only show row titles if requested >>>
            if show_row_titles and j == 0:
                ax.text(-0.35, 0.5, row_titles[i], transform=ax.transAxes, ha='right', va='center', rotation=90, fontsize=20)
            
            if i == 0:
                ax.set_title(column_titles[j], fontsize=14, pad=10) # Adjusted fontsize
            
            ax.tick_params(axis='x', labelbottom=(i == n_rows - 1))
            ax.tick_params(axis='y', labelleft=(j == 0))
            ax.grid(True, linestyle='--', alpha=0.5)

    # 5. Add colorbars for the panel
    # <<< MODIFIED: Only show colorbars if requested >>>
    if show_colorbars:
        for row_idx, mappable in enumerate(row_mappables):
            if mappable:
                cbar_ticks = np.linspace(row_min[row_idx], row_max[row_idx], 6)
                cbar = subfig.colorbar(mappable, ax=axes[row_idx, :].tolist(), location='right',
                                       orientation='vertical', shrink=0.8, aspect=20, pad=0.03, ticks=cbar_ticks)
                
                # <<< MODIFIED: Add unit label to the middle colorbar only >>>
                if row_idx == n_rows // 2:
                    cbar.set_label(f'({variable_units})', size=cbar_label_fontsize, labelpad=15)
                
                cbar.ax.tick_params(labelsize=cbar_tick_fontsize)
                
    subfig.text(0.02, 0.98, f"({index})", fontsize=20, fontweight='bold', va='top', ha='left')
    subfig.supxlabel("Longitude (°E)", fontsize=18)


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
fig = plt.figure(figsize=(24, 8), constrained_layout=True) # Adjusted size for better layout
subfigs = fig.subfigures(nrows=1, ncols=2, wspace=0.05)

# 3. Call the plotting function for each subfigure/SSP
ssp_names = list(ssp_configs.keys())

# <<< MODIFIED: Call left panel with row titles, but NO colorbars >>>
plot_ssp_panel(subfigs[0], ssp_configs[ssp_names[0]], india_map, index='a',
               show_row_titles=True, show_colorbars=False)

# <<< MODIFIED: Call right panel with colorbars, but NO row titles >>>
plot_ssp_panel(subfigs[1], ssp_configs[ssp_names[1]], india_map, index='b',
               show_row_titles=False, show_colorbars=True)

# <<< MODIFIED: Add single, shared labels for the entire figure >>>
fig.supylabel("Latitude (°N)", fontsize=18)


# 4. Save the final figure
print(f"\nSaving combined figure to {output_filename}...")
plt.savefig(output_filename, dpi=output_dpi, bbox_inches='tight')
print("Figure saved successfully.")

# 5. Show plot
plt.show()