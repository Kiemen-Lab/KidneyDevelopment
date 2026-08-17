# Figure 2 – Zones

This folder contains MATLAB scripts used to divide reconstructed kidneys into spatial zones and quantify cellular composition across these zones, as presented in Figure 2.

## Requirements

* MATLAB R2023b
* Image Processing Toolbox

## Scripts

| Script                           | Description                                                                                     |
| -------------------------------- | ----------------------------------------------------------------------------------------------- |
| `r_make_zones_volumes.m`         | Generates volumetric zone masks from the reconstructed tissue volumes.                          |
| `r_compute_compositions_zones.m` | Calculates the percentage of each of 16 cellular/tissue classes within each of the three zones. |
| `r_compute_percentage_zone.m`    | Calculates the relative volume occupied by each zone within the reconstructed tissue.           |
| `combine_z_projections.m`        | Combines zone projections for visualization.                                                    |
| `create_stacked_plot.m`          | Generates stacked plots of cellular composition across zones.                                   |
| `make_cmap_legend.m`             | Generates the color legend used for zone/composition visualizations.                            |

## Input data

The analysis requires the processed 3D reconstruction volumes and corresponding zone masks. The scripts expect sample-specific `.mat` files containing the reconstructed volume (`volTA`) and zone-volume files containing `volzone`.

Before running the scripts, update the input and output paths and sample names to match your local file structure.

## Usage

1. Generate the zone volumes using `r_make_zones_volumes.m`.
2. Run `r_compute_compositions_zones.m` to calculate cellular composition within each zone.
3. Run `r_compute_percentage_zone.m` to calculate the relative size of the zones.
4. Use the plotting scripts to generate the corresponding visualizations.

## Output

The analysis generates zone masks, cellular composition measurements, zone percentages, and stacked plots used in the spatial-zone analyses in Figure 2.
