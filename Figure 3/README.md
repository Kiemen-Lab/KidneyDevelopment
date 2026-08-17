# Figure 3 – Glomerular and Vascular Architecture

This folder contains MATLAB scripts used to quantify and visualize glomerular density, glomerular volume, and the spatial relationship between glomeruli and the renal vascular network, as presented in Figure 3.

## Requirements

* MATLAB R2023b
* Image Processing Toolbox
* Statistics and Machine Learning Toolbox

## Scripts

| Script                                          | Description                                                                                                         |
| ----------------------------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| `r_make_volgloms.m`                             | Generates volumetric representations of the glomeruli.                                                              |
| `r_compute_glomeruli_density.m`                 | Calculates glomerular density across the reconstructed kidney.                                                      |
| `r_get_dist_centroid.m`                         | Calculates distances between glomerular centroids and the kidney surface.                                           |
| `r_glom_vess_fraction_distance_surface.m`       | Quantifies the vascular fraction in relation to glomerular position and distance from the kidney surface.           |
| `r_make_glom_heatmap_density_good.m`            | Generates heatmaps showing the spatial distribution and density of glomeruli.                                       |
| `r_make_heatmap_volume.m`                       | Generates heatmaps of glomerular volume.                                                                            |
| `r_plot_gloms_arteries.m`                       | Visualizes the spatial relationship between glomeruli and arteries.                                                 |
| `r_plot_arteries_arterioles.m`                  | Visualizes the reconstructed arterial and arteriolar networks.                                                      |
| `r_mega_plot_as_function_of_distance.m`         | Generates combined plots of glomerular and vascular measurements as a function of distance from the kidney surface. |
| `r_combined_plots_distribs_densities_volumes.m` | Combines distributions, density, and volume measurements for visualization.                                         |

## Usage

1. Ensure the required 3D reconstruction, glomerular segmentation, and vascular data are available.
2. Update the input and output paths in the scripts to match your local file structure.
3. Run the scripts required to generate glomerular volumes and spatial measurements.
4. Use the plotting scripts to generate the corresponding visualizations.

## Output

The scripts generate glomerular density and volume measurements, distance-based spatial analyses, heatmaps, and visualizations of the glomerular and vascular networks used in Figure 3.
