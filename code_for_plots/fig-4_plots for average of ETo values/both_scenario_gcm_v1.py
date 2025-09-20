import matplotlib.pyplot as plt
import matplotlib as mpl
import xarray as xr
import numpy as np
import os
import geopandas as gpd
from matplotlib.patches import PathPatch # For clipping
from matplotlib.path import Path
import fiona
import matplotlib.ticker as mticker

# --- Reusable Plotting Function (Now with new features) ---
def create_grid_plot(subfig, netcdf_files, row_titles, column_titles, row_cmaps, subfigure_title, common_config, index, show_colorbar, show_row_titles):
    """
    Generates a 6x5 grid of plots on a given Matplotlib subfigure object.

    Args:
        ... (previous args)
        show_colorbar (bool): If True, display colorbars for each row.
        show_row_titles (bool): If True, display titles for each row.
    """
    n_rows = common_config['n_rows']
    n_cols = common_config['n_cols']

    # --- 1. Calculate or Define Min/Max for each row ---
    print(f"\n--- Calculating/Defining Min/Max for Subfigure: {subfigure_title} ---")
    data_threshold = common_config['data_threshold']
    variable_name = common_config['variable_name']

    manual_row_clim = common_config.get('manual_row_clim')
    row_min = [np.inf] * n_rows
    row_max = [-np.inf] * n_rows

    for row_idx in range(n_rows):
        if manual_row_clim and row_idx < len(manual_row_clim) and manual_row_clim[row_idx] is not None:
            row_min[row_idx], row_max[row_idx] = manual_row_clim[row_idx]
            print(f"Row {row_idx} using manual range: min={row_min[row_idx]:.2f}, max={row_max[row_idx]:.2f}")
        else:
            found_data_in_row = False
            for j in range(len(netcdf_files[row_idx])):
                nc_file = netcdf_files[row_idx][j]
                if not os.path.exists(nc_file):
                    print(f"Warning: File not found, skipping for min/max: {nc_file}")
                    continue
                try:
                    with xr.open_dataset(nc_file) as ds:
                        if variable_name not in ds.variables: continue
                        data = ds[variable_name].isel(time=0) if 'time' in ds[variable_name].dims else ds[variable_name]
                        data_filtered = data.where(data <= data_threshold)
                        current_min, current_max = data_filtered.min().item(), data_filtered.max().item()
                        if not np.isnan(current_min):
                            row_min[row_idx] = min(row_min[row_idx], current_min)
                            found_data_in_row = True
                        if not np.isnan(current_max):
                            row_max[row_idx] = max(row_max[row_idx], current_max)
                except Exception as e:
                    print(f"Error reading {nc_file} for min/max: {e}. Skipping.")

            if not found_data_in_row:
                row_min[row_idx], row_max[row_idx] = 0, 1
            if row_min[row_idx] == row_max[row_idx]: row_max[row_idx] += 1e-6
            print(f"Row {row_idx} calculated range: min={row_min[row_idx]:.2f}, max={row_max[row_idx]:.2f}")


    # --- 2. Create Axes and Plot ---
    axes = subfig.subplots(nrows=n_rows, ncols=n_cols, sharex=True, sharey=True)
    subfig.suptitle(subfigure_title, fontsize=20, fontweight='bold')
    subfig.supxlabel("Longitude (°E)", fontsize=common_config['axis_label_fontsize'])

    # <<< MODIFIED: Removed individual supxlabel from here. It will be added once for the whole figure.

    row_mappables = [None] * n_rows
    india_map = common_config['india_map']

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
                    data = ds[variable_name].isel(time=0) if 'time' in ds[variable_name].dims else ds[variable_name]
                    lons, lats = np.meshgrid(ds[common_config['lon_name']], ds[common_config['lat_name']])
                    data_filtered = data.where(data <= data_threshold).values

                    levels = np.linspace(row_min[i], row_max[i], common_config['num_contour_levels'])
                    norm = mpl.colors.Normalize(vmin=row_min[i], vmax=row_max[i])
                    current_cmap = row_cmaps[i]

                    india_map.plot(ax=ax, facecolor='none', edgecolor=common_config['shapefile_outline_color'], linewidth=0.5)
                    ax.tick_params(axis='both', labelsize=common_config['axis_tick_fontsize'])
                    contour_plot = ax.contourf(lons, lats, data_filtered, levels=levels, cmap=current_cmap, norm=norm, extend='both')

                    row_mappables[i] = contour_plot
                    india_map.plot(ax=ax, facecolor='none', edgecolor=common_config['shapefile_outline_color'], linewidth=0.25)

                    ax.grid(True, linestyle='--', alpha=0.5)

                    if i < n_rows - 1: ax.tick_params(axis='x', labelbottom=False)
                    if j > 0: ax.tick_params(axis='y', labelleft=False)

                    if j == 0:  # Leftmost column
                        ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%d°')) # Added degree symbol
                    if i == n_rows - 1:  # Bottom row
                        ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%d°')) # Added degree symbol

            except Exception as e:
                print(f"!!! Major Error plotting {os.path.basename(nc_file)}: {e}")
                ax.set_facecolor('salmon')

    # --- 3. Add Titles and Colorbars (Now Conditional) ---
    for j, title in enumerate(column_titles):
        axes[0, j].set_title(title, fontsize=common_config['col_title_fontsize'], pad=10)

    # <<< MODIFIED: Added specific axis label placement here
    if show_row_titles:
        for i, title in enumerate(row_titles):
            axes[i, 0].text(common_config['row_title_xpos'], 0.5, title, transform=axes[i, 0].transAxes,
                            ha='right', va='center', rotation=0, fontsize=common_config['row_title_fontsize'])
        
        # Add the shared Y-axis label to the right of the row titles
        # We anchor it to the middle-left axis for good vertical centering
        label_ax = axes[n_rows // 2, 0]
        label_ax.text(
            -0.4,  # X-position in Axes coords. -0.6 is row title, 0.0 is the plot. This is in between.
            0.5,    # Y-position in Axes coords (vertically centered)
            "Latitude (°N)",
            transform=label_ax.transAxes,
            ha='center',
            va='center',
            rotation='vertical',
            fontsize=common_config['axis_label_fontsize']
        )


    if show_colorbar:
        for row_idx, mappable in enumerate(row_mappables):
            if mappable:
                cbar = subfig.colorbar(mappable, ax=axes[row_idx, :].tolist(), location='right',
                                       orientation='vertical', shrink=0.8, aspect=20, pad=0.02,
                                       ticks=np.linspace(row_min[row_idx], row_max[row_idx], 6))
                
                if row_idx == n_rows // 2:
                    cbar.set_label(f"({common_config['variable_units']})", fontsize=common_config['cbar_label_fontsize'], labelpad=15)

                cbar.ax.tick_params(labelsize=common_config['cbar_tick_fontsize'])

    subfig.text(0.03, 0.98, f'({index})', ha='right', va='top', fontsize=20, fontweight='bold')


if __name__ == "__main__":
    # --- Main Configuration ---
    common_config = {
        'variable_name': 'decadal_avg', 'lon_name': 'lon', 'lat_name': 'lat',
        'variable_units': "mm day" + r'$^{-1}  $', 'output_filename': r'Fig2\figure_ssp126_ssp585_finalv8.png',
        'output_dpi': 300, 'india_shapefile_path': "D:\CO2 Paper\Latest_frontiers_in_water\Shapefiles\india_state_boundary_projected.shp",
        'shapefile_outline_color': 'black', 'attribute_column_name': 'State_Name',
        'island_names_to_exclude': ['Lakshadweep', 'Andaman & Nicobar'],
        'n_rows': 6, 'n_cols': 5, 'row_cmaps': ["BuPu", "Greens", "Blues", "Oranges", "Reds", "YlOrBr"],
        'row_title_xpos': -0.6, 'data_threshold': 200.0, 'num_contour_levels': 6,
        'row_title_fontsize': 20, 'col_title_fontsize': 20, 'cbar_label_fontsize': 16,
        'cbar_tick_fontsize': 16, 'axis_tick_fontsize': 20,
        'axis_label_fontsize': 18,

        'manual_row_clim': [
            (0, 7.5), (0, 6.0),
            (0, 7.5), (0, 6.0),
            (0, 8.5), (0, 6.0),
        ]
    }

    # --- Data and Titles Specific to each Subplot ---
    row_titles = [
        r'$ET_o^{original}  $' + '\n' + '2021-2030', r"$ET_o^{CO_2}     $" + "\n2021-2030",
        r'$ET_o^{original}  $' + '\n' + '2051-2060', r"$ET_o^{CO_2}     $" + "\n2051-2060",
        r'$ET_o^{original}  $' + '\n' + '2091-2100', r"$ET_o^{CO_2}    $" + "\n2091-2100"
    ]
    column_titles = ["GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL"]

    base_path = "D:\CO2 Paper\plot_codes\Fig2"

    # --- FILE PATHS FOR SUBPLOT (a) ---
    subfigure_title_a = "Scenario: SSP1-2.6"
    netcdf_files_a = []
    for period in ["2021-2030", "2051-2060", "2091-2100"]:
        netcdf_files_a.append([os.path.join(base_path, f"eto_files_ssp126\decadal_average_{gcm}_eto_{period}_with_FAO-PM_eq.nc") for gcm in column_titles])
        netcdf_files_a.append([os.path.join(base_path, f"eto_files_ssp126\decadal_average_{gcm}_eto_{period}_with_CO2.nc") for gcm in column_titles])

    # --- FILE PATHS FOR SUBPLOT (b) ---
    subfigure_title_b = "Scenario: SSP5-8.5"
    netcdf_files_b = []
    for period in ["2021-2030", "2051-2060", "2091-2100"]:
        netcdf_files_b.append([os.path.join(base_path, f"eto_files_ssp585\decadal_average_{gcm}_eto_{period}_with_FAO-PM_eq.nc") for gcm in column_titles])
        netcdf_files_b.append([os.path.join(base_path, f"eto_files_ssp585\decadal_average_{gcm}_eto_{period}_with_CO2.nc") for gcm in column_titles])

    # --- Load Shapefile Once ---
    with fiona.Env(SHAPE_RESTORE_SHX='YES'):
        india_map = gpd.read_file(common_config['india_shapefile_path'])
    india_map = india_map.to_crs(epsg=4326)
    india_map = india_map[~india_map[common_config['attribute_column_name']].isin(common_config['island_names_to_exclude'])]
    common_config['india_map'] = india_map

    # --- Create Figure and Subfigures Layout ---
    fig = plt.figure(figsize=(28, 14), constrained_layout=True)
    subfigs = fig.subfigures(nrows=1, ncols=2, wspace=0.07)

    # --- Generate the Plots ---
    print("--- Generating Subplot (a) [LEFT] ---")
    create_grid_plot(subfigs[0], netcdf_files_a, row_titles, column_titles,
                     common_config['row_cmaps'], subfigure_title_a, common_config,
                     index='a', show_colorbar=False, show_row_titles=True)

    print("\n--- Generating Subplot (b) [RIGHT] ---")
    create_grid_plot(subfigs[1], netcdf_files_b, row_titles, column_titles,
                     common_config['row_cmaps'], subfigure_title_b, common_config,
                     index='b', show_colorbar=True, show_row_titles=False)

    # <<< MODIFIED: Add a single, shared x label for the ENTIRE figure
    


    # --- Save the Final Figure ---
    print(f"\nSaving final figure to {common_config['output_filename']}...")
    plt.savefig(common_config['output_filename'], dpi=common_config['output_dpi'], bbox_inches='tight')
    print("Figure saved.")

    plt.show()