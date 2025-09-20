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

# --- Configuration ---

netcdf_dir = "D:\CO2 Paper\plot_codes\Fig2\eto_files"
variable_name = 'decadal_avg'
lon_name = 'lon'
lat_name = 'lat'
variable_units = "mm day" + r'$^{-1}  $'

output_filename = 'D:\CO2 Paper\plot_codes\Fig2\decadal_avg_plot_ssp585_font.png' # <<< MODIFIED: New filename
output_dpi = 300
india_shapefile_path = "D:\CO2 Paper\plot_codes\Shapefiles\india_state_boundary_projected.shp"
shapefile_outline_color = 'black'
attribute_column_name = 'State_Name'
island_names_to_exclude = ['Lakshadweep', 'Andaman & Nicobar']

n_rows = 6 # <<< MODIFIED
n_cols = 5
n_plots = n_rows * n_cols

# Row-Specific Colormaps
# <<< MODIFIED: Added two more colormaps for the new rows
row_cmaps = ["Purples", "Reds", "viridis_r", "Blues", "Greens",  "YlOrBr"]
if len(row_cmaps) < n_rows:
    raise ValueError(f"Need at least {n_rows} colormaps defined in `row_cmaps`.")

figure_left_margin = 0.10
row_title_xpos = -0.35
data_threshold = 200.0

# Contourf specific configuration
num_contour_levels = 6

# <<< MODIFIED: Added titles for the 2051-2060 period in chronological order
row_titles = [
    r'$ET_o^{original}  $' + '\n' + '2021-2030',
    r"$ET_o^{CO_2}     $" + "\n2021-2030",
    r'$ET_o^{original}  $' + '\n' + '2051-2060', # New Row
    r"$ET_o^{CO_2}     $" + "\n2051-2060", # New Row
    r'$ET_o^{original}  $' + '\n' + '2091-2100',
    r"$ET_o^{CO_2}    $" + "\n2091-2100"
]
column_titles = ["GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL"]
row_title_fontsize = 20
col_title_fontsize = 20
colorbar_label = f'({variable_units})'

# <<< MODIFIED: Added two new rows of file paths for the 2051-2060 period
netcdf_files = [
    # Row 1: 2021-2030 Original
    ["D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_GFDL-ESM4_eto_2021-2030_with_FAO-PM_eq.nc",
     "D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_IPSL-CM6A-LR_eto_2021-2030_with_FAO-PM_eq.nc",
     "D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_MPI-ESM1-2-HR_eto_2021-2030_with_FAO-PM_eq.nc",
     "D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_MRI-ESM2-0_eto_2021-2030_with_FAO-PM_eq.nc",
     "D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_UKESM1-0-LL_eto_2021-2030_with_FAO-PM_eq.nc"],
    # Row 2: 2021-2030 CO2
    ["D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_GFDL-ESM4_eto_2021-2030_with_CO2.nc",
     "D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_IPSL-CM6A-LR_eto_2021-2030_with_CO2.nc",
     "D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_MPI-ESM1-2-HR_eto_2021-2030_with_CO2.nc",
     "D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_MRI-ESM2-0_eto_2021-2030_with_CO2.nc",
     "D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_UKESM1-0-LL_eto_2021-2030_with_CO2.nc"],
    # Row 3: 2051-2060 Original (NEW)
    ["D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_GFDL-ESM4_eto_2051-2060_with_FAO-PM_eq.nc",
     "D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_IPSL-CM6A-LR_eto_2051-2060_with_FAO-PM_eq.nc",
     "D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_MPI-ESM1-2-HR_eto_2051-2060_with_FAO-PM_eq.nc",
     "D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_MRI-ESM2-0_eto_2051-2060_with_FAO-PM_eq.nc",
     "D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_UKESM1-0-LL_eto_2051-2060_with_FAO-PM_eq.nc"],
    # Row 4: 2051-2060 CO2 (NEW)
    ["D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_GFDL-ESM4_eto_2051-2060_with_CO2.nc",
     "D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_IPSL-CM6A-LR_eto_2051-2060_with_CO2.nc",
     "D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_MPI-ESM1-2-HR_eto_2051-2060_with_CO2.nc",
     "D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_MRI-ESM2-0_eto_2051-2060_with_CO2.nc",
     "D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_UKESM1-0-LL_eto_2051-2060_with_CO2.nc"],
    # Row 5: 2091-2100 Original
    ["D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_GFDL-ESM4_eto_2091-2100_with_FAO-PM_eq.nc",
     "D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_IPSL-CM6A-LR_eto_2091-2100_with_FAO-PM_eq.nc",
     "D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_MPI-ESM1-2-HR_eto_2091-2100_with_FAO-PM_eq.nc",
     "D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_MRI-ESM2-0_eto_2091-2100_with_FAO-PM_eq.nc",
     "D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_UKESM1-0-LL_eto_2091-2100_with_FAO-PM_eq.nc"],
    # Row 6: 2091-2100 CO2
    ["D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_GFDL-ESM4_eto_2091-2100_with_CO2.nc",
     "D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_IPSL-CM6A-LR_eto_2091-2100_with_CO2.nc",
     "D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_MPI-ESM1-2-HR_eto_2091-2100_with_CO2.nc",
     "D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_MRI-ESM2-0_eto_2091-2100_with_CO2.nc",
     "D:\CO2 Paper\plot_codes\Fig2\eto_files_ssp585\decadal_average_UKESM1-0-LL_eto_2091-2100_with_CO2.nc"]
]

# --- Main Script (No changes needed below this line) ---

print(f"Calculating global min/max values (ignoring values > {data_threshold})...")
row_min = [np.inf for _ in range(n_rows)] # Use np.inf for safer initialization
row_max = [-np.inf for _ in range(n_rows)] # Use -np.inf
for row_idx in range(n_rows):
    found_data_in_row = False
    for j in range(len(netcdf_files[row_idx])):
        nc_file = netcdf_files[row_idx][j]
        if not os.path.exists(nc_file):
            print(f"Warning: File not found, skipping for min/max calculation: {nc_file}")
            continue
        try:
            with xr.open_dataset(nc_file) as ds:
                if variable_name not in ds.variables: continue
                data = ds[variable_name]
                if 'time' in data.dims: data = data.isel(time=0)
                data_filtered = data.where(data <= data_threshold)
                current_min = data_filtered.min().item()
                current_max = data_filtered.max().item()

                if not np.isnan(current_min):
                    row_min[row_idx] = min(row_min[row_idx], current_min)
                    found_data_in_row = True
                if not np.isnan(current_max):
                    row_max[row_idx] = max(row_max[row_idx], current_max)
        except Exception as e:
            print(f"Error reading {nc_file} for min/max: {e}. Skipping.")
    if not found_data_in_row:
        print(f"Warning: No valid data found for row {row_idx}. Setting default range 0-1.")
        row_min[row_idx] = 0
        row_max[row_idx] = 1
    if row_min[row_idx] == row_max[row_idx]: row_max[row_idx] += 1e-6 # Avoid zero range
    print(f"Row {row_idx} filtered data range: min={row_min[row_idx]:.2f}, max={row_max[row_idx]:.2f}")


with fiona.Env(SHAPE_RESTORE_SHX='YES'):
    india_map = gpd.read_file(india_shapefile_path)
print(f"Shapefile loaded successfully. Original CRS: {india_map.crs}")

if india_map.crs and not india_map.crs.is_geographic:
    print("Warning: Shapefile CRS does not appear to be geographic (lat/lon).")
    print("Overlay might be inaccurate if it doesn't match data coordinates.")
india_map = india_map.to_crs(epsg=4326)
india_map = india_map[~india_map[attribute_column_name].isin(island_names_to_exclude)]

fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols,
                         figsize=(n_cols * 2.5 + 2.5, n_rows * 2.5 + 1), # Adjusted figsize for 6 rows
                         sharex=True, sharey=True)

row_mappables = [None] * n_rows

cbar_label_fontsize = 16
cbar_tick_fontsize = 16
axis_tick_fontsize = 20

for i in range(n_rows):
    for j in range(n_cols):
        ax = axes[i, j]
        if i >= len(netcdf_files) or j >= len(netcdf_files[i]):
            ax.set_visible(False)
            continue

        nc_file = netcdf_files[i][j]
        print(f"Processing {os.path.basename(nc_file)}...")
        
        if not os.path.exists(nc_file):
            print(f"Error: File not found -> {nc_file}. Skipping plot.")
            ax.text(0.5, 0.5, 'File Not Found', ha='center', va='center', transform=ax.transAxes, color='red')
            ax.set_facecolor('lightgrey')
            continue

        try:
            with xr.open_dataset(nc_file) as ds:
                if variable_name not in ds.variables:
                    ax.set_facecolor('lightgrey'); continue

                data = ds[variable_name]
                ax.xaxis.grid(True, linestyle='--', alpha=0.5)
                ax.yaxis.grid(True, linestyle='--', alpha=0.5)

                current_cmap = row_cmaps[i]
                try:
                    lats = ds[lat_name]; lons = ds[lon_name]
                    if lons.ndim == 1 and lats.ndim == 1: lons, lats = np.meshgrid(lons.values, lats.values)
                    elif lons.ndim == 2 and lats.ndim == 2: lons = lons.values; lats = lats.values
                    else: raise ValueError("Unsupported coordinate dimensions")
                except (KeyError, ValueError) as e:
                    ax.set_facecolor('lightgrey'); print(f"Coord error in {nc_file}: {e}"); continue

                if 'time' in data.dims: data = data.isel(time=0)
                data_filtered = data.where(data <= data_threshold).values

                levels = np.linspace(row_min[i], row_max[i], num_contour_levels)
                norm = mpl.colors.Normalize(vmin=row_min[i], vmax=row_max[i])

                middle_color = mpl.colors.to_hex(plt.cm.get_cmap(current_cmap)(0.5))
                india_map.plot(ax=ax,
                               facecolor=middle_color,
                               edgecolor=shapefile_outline_color,
                               linewidth=0.25)
                contour_plot = ax.contourf(lons, lats, data_filtered,
                                        levels=levels,
                                        cmap=current_cmap,
                                        norm=norm,
                                        extend='both')
                row_mappables[i] = contour_plot
                india_map.plot(ax=ax,
                               facecolor='none',
                               edgecolor=shapefile_outline_color,
                               linewidth=0.25)
                ax.tick_params(axis='both', labelsize=axis_tick_fontsize)
                if i < n_rows - 1: ax.tick_params(axis='x', labelbottom=False)
                if j > 0: ax.tick_params(axis='y', labelleft=False)

        except Exception as e:
            print(f"Error plotting {nc_file}: {e}")
            import traceback; print(traceback.format_exc())
            ax.set_facecolor('salmon')

# Column Titles
for j, title in enumerate(column_titles):
    if j < n_cols and axes[0, j].get_visible():
        axes[0, j].set_title(title, fontsize=col_title_fontsize, pad=10)
# Row Titles
for i, title in enumerate(row_titles):
     if i < n_rows and axes[i, 0].get_visible():
        axes[i, 0].text(row_title_xpos, 0.5, title, transform=axes[i, 0].transAxes,
                        ha='right', va='center', rotation=0, fontsize=row_title_fontsize)


# Add Colorbar for Each Row
for row_idx in range(n_rows):
    mappable = row_mappables[row_idx]
    if mappable is not None:
        cbar = fig.colorbar(mappable,
                    ax=axes[row_idx, :].tolist(),
                    label=colorbar_label,
                    location='right',
                    orientation='vertical',
                    shrink=0.8,
                    aspect=20,
                    pad=0.02,
                    ticks=np.linspace(row_min[row_idx], row_max[row_idx], 6))
        cbar.set_label(colorbar_label, fontsize=cbar_label_fontsize) # <<< MODIFY THIS LINE
        cbar.ax.tick_params(labelsize=cbar_tick_fontsize)
fig.suptitle('Scenario : SSP5-8.5', fontsize=20, fontweight='bold', y=0.93, x=0.45)
fig.text(0.03, 0.93, '(b)', ha='right', va='top', fontsize=20, fontweight='bold')
print(f"Saving figure to {output_filename}...")
plt.savefig(output_filename, dpi=output_dpi, bbox_inches='tight')
print("Figure saved.")

plt.show()