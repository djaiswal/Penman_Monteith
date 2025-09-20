import matplotlib.pyplot as plt
import matplotlib as mpl
import xarray as xr
import numpy as np
import os
import geopandas as gpd
from matplotlib.patches import PathPatch # For clipping
from matplotlib.path import Path
import fiona
from scipy.interpolate import griddata

# --- Main Configuration ---
variable_name = 'decadal_avg'
lon_name = 'lon'
lat_name = 'lat'
variable_units ="mm day" + r'$^{-1}$'
# <<< CHANGE 1: Using a raw string (r'...') for the Windows path to avoid errors >>>
output_filename = r'D:\CO2 Paper\plot_codes\Fig3\difference_over_decade_ssp126_font.png'
output_dpi = 300
data_threshold = 200.0 # Values above this will be ignored for min/max calculation

# --- Shapefile Configuration ---
india_shapefile_path = r"D:\CO2 Paper\Latest_frontiers_in_water\Shapefiles\india_state_boundary_projected.shp"
shapefile_outline_color = 'black'
attribute_column_name = 'State_Name'  # Column name in shapefile for state names
island_names_to_exclude = ['Lakshadweep', 'Andaman & Nicobar'] # Islands to exclude from mainland plot

# --- Plot Layout Configuration ---
n_rows = 3
n_cols = 5
n_plots = n_rows * n_cols
row_titles = ['2021-2030', '2051-2060', '2091-2100']
column_titles = ["GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL"]
row_cmaps = ["Oranges", "Greens", "Blues"]  # Colormaps for each row
if len(row_cmaps) < n_rows:
    raise ValueError(f"Need at least {n_rows} colormaps defined in `row_cmaps`.")

# --- File Paths ---
base_netcdf_dir = r"D:\CO2 Paper\Results_units_corrected\ssp126\decadal_average"
gcms = ["GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL"]
periods = ["2021-2030", "2051-2060", "2091-2100"]

netcdf_files = []
for period in periods:
    row_files = [os.path.join(base_netcdf_dir, f"difference_average_{gcm}_{period}_daily_difference.nc") for gcm in gcms]
    netcdf_files.append(row_files)

# --- Plotting Style Configuration ---
row_title_fontsize = 20
col_title_fontsize = 20
colorbar_label = f'({variable_units})'
num_contour_levels = 6
row_title_xpos = -0.35 # Adjust as needed

# --- SCRIPT LOGIC ---

# 1. Calculate row-specific min/max values
print(f"Calculating row-specific min/max values (ignoring values > {data_threshold})...")
row_min = [np.inf] * n_rows
row_max = [-np.inf] * n_rows

for row_idx in range(n_rows):
    found_data_in_row = False
    for nc_file in netcdf_files[row_idx]:
        try:
            with xr.open_dataset(nc_file) as ds:
                if variable_name not in ds.variables:
                    print(f"Warning: Variable '{variable_name}' not found in {os.path.basename(nc_file)}. Skipping.")
                    continue
                data = ds[variable_name]
                if 'time' in data.dims:
                    data = data.isel(time=0)

                data_filtered = data.where(data <= data_threshold)
                current_min = data_filtered.min().item()
                current_max = data_filtered.max().item()

                if not np.isnan(current_min):
                    row_min[row_idx] = min(row_min[row_idx], current_min)
                    found_data_in_row = True
                if not np.isnan(current_max):
                    row_max[row_idx] = max(row_max[row_idx], current_max)
                    found_data_in_row = True
        except Exception as e:
            print(f"Error reading {nc_file} for min/max: {e}. Skipping.")

    if not found_data_in_row:
        print(f"Warning: No valid data found for row {row_idx}. Using default range [0, 1].")
        row_min[row_idx], row_max[row_idx] = 0, 1
    if row_min[row_idx] == row_max[row_idx]: # Avoid zero range
        row_max[row_idx] += 1e-6
    print(f"Row {row_idx} ('{row_titles[row_idx]}') filtered data range: min={row_min[row_idx]:.2f}, max={row_max[row_idx]:.2f}")


# 2. Load and prepare shapefile
with fiona.Env(SHAPE_RESTORE_SHX='YES'):
    india_map = gpd.read_file(india_shapefile_path)
print(f"Shapefile loaded. Original CRS: {india_map.crs}")
india_map = india_map.to_crs(epsg=4326) # Ensure shapefile is in lat/lon CRS
india_map = india_map[~india_map[attribute_column_name].isin(island_names_to_exclude)] # Filter out islands
print("Shapefile reprojected to EPSG:4326 and islands excluded.")

# 3. Create the figure and axes array
# <<< CHANGE 2: Added `constrained_layout=True` for automatic layout management >>>
# and adjusted figsize for a better aspect ratio.
fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols,
                         figsize=(n_cols * 2.5 + 2.5, n_rows * 2.5 + 1),
                         sharex=True, sharey=True,
                         constrained_layout=True)

row_mappables = [None] * n_rows
cbar_label_fontsize = 16
cbar_tick_fontsize = 16
axis_tick_fontsize = 20
# 4. Loop through files and plot on each axis
for i in range(n_rows):
    for j in range(n_cols):
        nc_file = netcdf_files[i][j]
        ax = axes[i, j]
        print(f"Plotting {os.path.basename(nc_file)} on subplot ({i},{j})...")

        try:
            with xr.open_dataset(nc_file) as ds:
                if variable_name not in ds.variables:
                    ax.set_facecolor('lightgrey'); continue

                data = ds[variable_name]
                lats = ds[lat_name]
                lons = ds[lon_name]

                if lons.ndim == 1 and lats.ndim == 1:
                    lons, lats = np.meshgrid(lons.values, lats.values)
                elif lons.ndim == 2 and lats.ndim == 2:
                    lons, lats = lons.values, lats.values
                else: raise ValueError("Unsupported coordinate dimensions")

                if 'time' in data.dims: data = data.isel(time=0)
                data_filtered = data.where(data <= data_threshold).values

                levels = np.linspace(row_min[i], row_max[i], num_contour_levels + 1)
                norm = mpl.colors.Normalize(vmin=row_min[i], vmax=row_max[i])
                current_cmap = row_cmaps[i]

                contour_plot = ax.contourf(lons, lats, data_filtered,
                                           levels=levels, cmap=current_cmap,
                                           norm=norm, extend='both')
                row_mappables[i] = contour_plot
                india_map.plot(ax=ax, facecolor='none',
                               edgecolor=shapefile_outline_color, linewidth=0.25)
                ax.tick_params(axis='both', labelsize=axis_tick_fontsize)
                if i < n_rows - 1: ax.tick_params(axis='x', labelbottom=False)
                if j > 0: ax.tick_params(axis='y', labelleft=False)
                ax.grid(True, linestyle='--', alpha=0.5)

        except Exception as e:
            print(f"Error plotting {nc_file}: {e}")
            import traceback; print(traceback.format_exc())
            ax.set_facecolor('salmon')

# 5. Hide unused axes (if any)
num_files_plotted = sum(len(row) for row in netcdf_files)
for k in range(num_files_plotted, n_plots):
     axes.flatten()[k].set_visible(False)

# 6. Add Row and Column Titles
for j, title in enumerate(column_titles):
    if j < n_cols and axes[0, j].get_visible():
        axes[0, j].set_title(title, fontsize=col_title_fontsize, pad=10)
for i, title in enumerate(row_titles):
     if i < n_rows and axes[i, 0].get_visible():
        axes[i, 0].text(row_title_xpos, 0.5, title, transform=axes[i, 0].transAxes,
                        ha='right', va='center', rotation=90, fontsize=row_title_fontsize)

# 7. Add a Colorbar for Each Row
for row_idx in range(n_rows):
    mappable = row_mappables[row_idx]
    if mappable is not None:
        cbar_ticks = np.linspace(row_min[row_idx], row_max[row_idx], 6)
        # <<< CHANGE 3: Reduced `pad` from 0.2 to a much more reasonable 0.03 >>>
        cbar = fig.colorbar(mappable,
                            ax=axes[row_idx, :].tolist(),
                            label=colorbar_label,
                            location='right',
                            orientation='vertical',
                            shrink=0.8,
                            aspect=20,
                            pad=0.03, # Reduced padding
                            ticks=cbar_ticks)
        cbar.set_label(colorbar_label, fontsize=cbar_label_fontsize) # <<< MODIFY THIS LINE
        cbar.ax.tick_params(labelsize=cbar_tick_fontsize)
        if row_idx != n_rows // 2: # Only label the middle colorbar
             cbar.set_label('')

# 8. Adjust layout and save the figure
# <<< CHANGE 4: Removed the `plt.subplots_adjust` line completely. `constrained_layout` handles it. >>>
fig.suptitle('Scenario : SSP1-2.6', fontsize=20, fontweight='bold')
fig.text(0.03, 0.98, '(a)', ha='right', va='top', fontsize=20, fontweight='bold')
print(f"Saving figure to {output_filename}...")
# Note: `bbox_inches='tight'` works well with constrained_layout
plt.savefig(output_filename, dpi=output_dpi, bbox_inches='tight')
print("Figure saved.")

# 10. Show plot
plt.show()