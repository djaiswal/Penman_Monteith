import matplotlib.pyplot as plt
import matplotlib as mpl
import xarray as xr
import numpy as np
import os
import geopandas as gpd
import fiona
import matplotlib.ticker as mticker

# --- Overall Configuration ---
output_filename = 'Fig5/ssp126_vs_ssp585_seasonal_comparison_font1_v10.png'
output_dpi = 300

# --- Shared Plotting Configuration ---
shared_config = {
    'variable_name': 'difference_decadal_avg',
    'lon_name': 'lon',
    'lat_name': 'lat',
    'variable_units': "mm day" + r'$^{-1}$',
    'data_threshold': 200.0,
    'num_contour_levels': 6,
    # Grid Layout
    'row_titles': ['2021-2030', '2051-2060', '2091-2100'],
    'column_titles': ["Monsoon\n(Jun-Sep)", "Post-Monsoon\n(Oct-Dec)", "Winter\n(Jan-Feb)", "Pre-Monsoon\n(Mar-May)"],
    'seasons_in_order': ["monsoon", "post_monsoon", "winter", "pre_monsoon"],
    'periods_in_order': ["2021-2030", "2051-2060", "2091-2100"],
    'row_cmaps': ["viridis_r", "PuBu", "cividis_r"],
    # Font Sizes
    'panel_title_fontsize': 20,
    'row_col_title_fontsize': 20,
    'axis_label_fontsize': 16, # Increased for the new text label
    'axis_tick_fontsize': 16,
    'cbar_label_fontsize': 20,
    'cbar_tick_fontsize': 16,
}
n_rows, n_cols = len(shared_config['row_titles']), len(shared_config['column_titles'])

# --- SSP-Specific Configurations ---
ssp_configs = {
    'ssp126': {
        'base_dir': r"D:\CO2 Paper\Results_units_corrected\ssp126\seasonal_difference",
        'panel_title': 'SSP1-2.6',
        'index': 'a',
        'filename_template': "{season}_{period}_difference.nc",
        'manual_row_clim': [
            (0, 1.75),
            (0, 2.5),
            (0, 4.25),
        ]
    },
    'ssp585': {
        'base_dir': r"D:\CO2 Paper\Results_units_corrected\ssp585\seasonal_difference",
        'panel_title': 'SSP5-8.5',
        'index': 'b',
        'filename_template': "{season}_{period}_difference.nc",
        'manual_row_clim': [
            (0, 1.75),
            (0, 2.5),
            (0, 4.25),
        ]
    }
}

# --- Shapefile Configuration ---
india_shapefile_path = r"D:\CO2 Paper\Latest_frontiers_in_water\Shapefiles\india_state_boundary_projected.shp"
shapefile_outline_color = 'black'
attribute_column_name = 'State_Name'
island_names_to_exclude = ['Lakshadweep', 'Andaman & Nicobar']

#=============================================================================
# REUSABLE PLOTTING FUNCTION
#=============================================================================
def plot_ssp_panel(subfig, ssp_config, india_map, shared_config, index, show_colorbar=True, show_yticklabels=True):
    """
    Creates a 3x4 grid of seasonal difference plots on a given subfigure.

    Args:
        ... (standard arguments) ...
        show_colorbar (bool): If True, adds colorbars to the right of the panel.
        show_yticklabels (bool): If True, shows y-axis tick labels and the Latitude label.
    """
    panel_title = ssp_config['panel_title']
    variable_name = shared_config['variable_name']

    subfig.suptitle(panel_title, fontsize=shared_config['panel_title_fontsize'], fontweight='bold')

    # 1. Construct file paths
    netcdf_files = []
    for period in shared_config['periods_in_order']:
        row_files = []
        for season in shared_config['seasons_in_order']:
            fpath = os.path.join(ssp_config['base_dir'], ssp_config['filename_template'].format(season=season, period=period))
            if os.path.exists(fpath):
                row_files.append(fpath)
            else:
                print(f"ERROR: Could not find file: {fpath}")
                row_files.append(None)
        netcdf_files.append(row_files)

    # 2. Calculate or define row-specific min/max values
    print(f"\n--- Calculating/Defining min/max for {panel_title} ---")
    manual_clim = ssp_config.get('manual_row_clim')
    row_min = [np.inf] * n_rows
    row_max = [-np.inf] * n_rows
    for row_idx in range(n_rows):
        if manual_clim and row_idx < len(manual_clim) and manual_clim[row_idx] is not None:
            row_min[row_idx], row_max[row_idx] = manual_clim[row_idx]
            print(f"Row {row_idx} ('{shared_config['periods_in_order'][row_idx]}') using manual range: min={row_min[row_idx]:.2f}, max={row_max[row_idx]:.2f}")
        else: # Auto-scaling fallback
            found_data_in_row = False
            for nc_file in netcdf_files[row_idx]:
                if nc_file is None: continue
                try:
                    with xr.open_dataset(nc_file) as ds:
                        if variable_name not in ds.variables:
                            print(f"Warning: Variable '{variable_name}' not found in {os.path.basename(nc_file)}")
                            continue
                        data = ds[variable_name].isel(time=0) if 'time' in ds[variable_name].dims else ds[variable_name]
                        data_filtered = data.where(data <= shared_config['data_threshold'])
                        current_min, current_max = data_filtered.min().item(), data_filtered.max().item()
                        print("gotchaaa")
                        print(f"Current min/max for {os.path.basename(nc_file)}: min={current_min:.2f}, max={current_max:.2f}")
                        if not np.isnan(current_min): row_min[row_idx] = min(row_min[row_idx], current_min); found_data_in_row = True
                        if not np.isnan(current_max): row_max[row_idx] = max(row_max[row_idx], current_max)
                except Exception as e:
                    print(f"Warning: Could not read {os.path.basename(nc_file)} for min/max: {e}")
            if not found_data_in_row: row_min[row_idx], row_max[row_idx] = 0, 1
            if row_min[row_idx] == row_max[row_idx]: row_max[row_idx] += 1e-6
            print(f"Row {row_idx} ('{shared_config['periods_in_order'][row_idx]}') calculated range: min={row_min[row_idx]:.2f}, max={row_max[row_idx]:.2f}")

    # 3. Create axes grid and plot
    axes = subfig.subplots(nrows=n_rows, ncols=n_cols, sharex=True, sharey=True)
    row_mappables = [None] * n_rows
    print(f"--- Plotting grid for {panel_title} ---")
    for i in range(n_rows):
        for j in range(n_cols):
            ax = axes[i, j]
            nc_file = netcdf_files[i][j]
            if nc_file is None:
                ax.set_facecolor('lightgrey')
            else:
                try:
                    with xr.open_dataset(nc_file) as ds:
                        if variable_name not in ds.variables: ax.set_facecolor('lightgrey'); continue
                        data = ds[variable_name].isel(time=0) if 'time' in ds[variable_name].dims else ds[variable_name]
                        lons, lats = ds[shared_config['lon_name']].values, ds[shared_config['lat_name']].values
                        if lons.ndim == 1: lons, lats = np.meshgrid(lons, lats)
                        
                        data_filtered = data.where(data <= shared_config['data_threshold']).values
                        current_min, current_max = np.nanmin(data_filtered).item(), np.nanmax(data_filtered).item()
                        # print("gotchaaa")
                        print(f"Current min/max for {os.path.basename(nc_file)}: min={current_min:.2f}, max={current_max:.2f}")

                        levels = np.linspace(row_min[i], row_max[i], shared_config['num_contour_levels'] + 1)
                        norm = mpl.colors.Normalize(vmin=row_min[i], vmax=row_max[i])
                        
                        contour_plot = ax.contourf(lons, lats, data_filtered, levels=levels, cmap=shared_config['row_cmaps'][i], norm=norm, extend='both')
                        row_mappables[i] = contour_plot
                        india_map.plot(ax=ax, facecolor='none', edgecolor=shapefile_outline_color, linewidth=0.3)
                except Exception as e:
                    print(f"ERROR plotting {os.path.basename(nc_file)}: {e}")
                    ax.set_facecolor('salmon')

            # Ticks, labels, and titles
            if j == 0:
                if index == 'a':
                    ax.text(-0.45, 0.5, shared_config['row_titles'][i], transform=ax.transAxes, ha='right', va='center', rotation=90, fontsize=shared_config['row_col_title_fontsize'])
            if i == 0:
                ax.set_title(shared_config['column_titles'][j], fontsize=shared_config['row_col_title_fontsize'], pad=10)
            
            # <<< MODIFIED: Control tick visibility and format based on function argument
            ax.tick_params(axis='x', labelbottom=(i == n_rows - 1), labelsize=shared_config['axis_tick_fontsize'])
            ax.tick_params(axis='y', labelleft=(j == 0 and show_yticklabels), labelsize=shared_config['axis_tick_fontsize'])
            
            if j == 0 and show_yticklabels:
                ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%d'))
            if i == n_rows - 1:
                ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%d'))
            ax.grid(True, linestyle='--', alpha=0.5)

    # 5. Add colorbars and shared labels
    # <<< MODIFIED: Only show colorbars if requested
    if show_colorbar:
        for row_idx, mappable in enumerate(row_mappables):
            if mappable:
                cbar_ticks = np.linspace(row_min[row_idx], row_max[row_idx], 6)
                cbar = subfig.colorbar(mappable, ax=axes[row_idx, :].tolist(), location='right',
                                       orientation='vertical', shrink=0.8, aspect=20, pad=0.03, ticks=cbar_ticks)
                cbar.ax.tick_params(labelsize=shared_config['cbar_tick_fontsize'])
                if row_idx == n_rows // 2: # Center the label on the middle colorbar
                    cbar.set_label(f"({shared_config['variable_units']})", size=shared_config['cbar_label_fontsize'])
    
    subfig.text(0.05, 0.98, f"({index})", fontsize=20, fontweight='bold', ha='center', va='center')
    
    # <<< MODIFIED: Add a single, shared x-label for the subfigure
    subfig.supxlabel("Longitude (°E)", fontsize=shared_config['axis_label_fontsize'])
    
    # <<< MODIFIED: Add Latitude label as rotated text on the left-middle axis if requested
    # if show_yticklabels:
    if index == 'a':
            # Add a vertical label on the left side of the middle row
        mid_row_ax = axes[n_rows // 2, 0]
        mid_row_ax.text(-0.28, 0.5, "Latitude (°N)",
                        transform=mid_row_ax.transAxes,
                        ha='center', va='center', rotation=90,
                        fontsize=shared_config['axis_label_fontsize'])


#=============================================================================
# MAIN SCRIPT EXECUTION
#=============================================================================
if __name__ == "__main__":
    print("Loading and preparing shapefile...")
    with fiona.Env(SHAPE_RESTORE_SHX='YES'):
        india_map = gpd.read_file(india_shapefile_path)
    india_map = india_map.to_crs(epsg=4326)
    india_map = india_map[~india_map[attribute_column_name].isin(island_names_to_exclude)]
    print("Shapefile ready.")

    fig = plt.figure(figsize=(20, 8), constrained_layout=True)
    subfigs = fig.subfigures(nrows=1, ncols=2, wspace=0.02)

    ssp_names = list(ssp_configs.keys())
    
    # --- MODIFIED: Control what is shown on each panel ---
    # Plot left panel (SSP1-2.6): Show Y-ticks and Latitude label, but NO colorbars
    plot_ssp_panel(subfigs[0], ssp_configs[ssp_names[0]], india_map, shared_config, 'a', 
                   show_colorbar=False, show_yticklabels=True)
    
    # Plot right panel (SSP5-8.5): Show colorbars, but NO Y-ticks or Latitude label
    plot_ssp_panel(subfigs[1], ssp_configs[ssp_names[1]], india_map, shared_config, 'b', 
                   show_colorbar=True, show_yticklabels=True)

    print(f"\nSaving combined figure to {output_filename}...")
    # Ensure the output directory exists
    os.makedirs(os.path.dirname(output_filename), exist_ok=True)
    plt.savefig(output_filename, dpi=output_dpi, bbox_inches='tight')
    print("Figure saved successfully.")

    plt.show()