# Figure 4 – Tubular Architecture

This folder contains MATLAB scripts used to quantify and visualize the volume and spatial distribution of renal tubules across developmental stages and species, as presented in Figure 4.

## Requirements

* MATLAB R2023b
* Image Processing Toolbox

## Scripts

| Script                                   | Description                                                                           |
| ---------------------------------------- | ------------------------------------------------------------------------------------- |
| `r_save_cleaned_tubules_volumes.m`       | Processes and saves cleaned 3D tubule volumes for the reconstructed samples.          |
| `r_compute_vf_tubules_distance.m`        | Calculates tubular volume fraction as a function of distance from the kidney surface. |
| `r_plot_tubules_fractions_per_species.m` | Plots tubular volume fractions across species.                                        |
| `r_plot_tubules_fractions_per_type.m`    | Plots tubular volume fractions by tubule type.                                        |
| `r_tubule_plot_E17.m`                    | Generates 3D visualizations of the E17 tubular reconstruction.                        |
| `r_tubule_plot_Macc.m`                   | Generates 3D visualizations of the macaque tubular reconstruction.                    |
| `r_make_legends.m`                       | Generates the legends used for the tubule visualizations.                             |

## Input data

The scripts require processed 3D reconstruction data and segmented/annotated tubule volumes. Before running the scripts, update the input and output paths in the scripts to match your local file structure.

## Usage

1. Generate the cleaned tubule volumes using `r_save_cleaned_tubules_volumes.m`.
2. Run `r_compute_vf_tubules_distance.m` to calculate tubular volume fraction by distance from the kidney surface.
3. Use the `r_plot_*` scripts to generate the quantitative plots.
4. Use the `r_tubule_plot_*` scripts to generate the 3D visualizations.

## Output

The scripts generate cleaned tubule volumes, distance-dependent tubular volume fractions, quantitative comparisons across species and tubule types, and 3D visualizations used in Figure 4.
