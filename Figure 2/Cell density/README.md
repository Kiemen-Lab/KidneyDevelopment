# Figure 2 – Cell Density

This folder contains MATLAB scripts used to calculate and visualize cell density as a function of distance from the kidney outer surface, as presented in Figure 2.

## Requirements

* MATLAB R2023b
* Image Processing Toolbox

## Scripts

| Script                                             | Description                                                                            |
| -------------------------------------------------- | -------------------------------------------------------------------------------------- |
| `r_celldensity_as_function_of_distance_blastema.m` | Calculates cell density across distance bins from the kidney surface for each sample.  |
| `r_make_correct_outershells_human1.m`              | Generates the corrected outer shell for the K1 human sample.                           |
| `r_make_correct_outershells_human3.m`              | Generates the corrected outer shell for the K3 human samples.                          |
| `r_make_correct_outershells_macaque.m`             | Generates the corrected outer shell for the macaque reconstruction.                    |
| `r_make_correct_outershells_mouse.m`               | Generates the corrected outer shell for the mouse reconstruction.                      |
| `r_plot_celldensity.m`                             | Plots cell-density distributions across samples and developmental stages.              |
| `r_verify_celldensity.m`                           | Performs quality control and verification of the calculated cell-density measurements. |

## Usage

1. Generate or obtain the required reconstructed volumes, cell masks, and outer-shell masks.
2. Update the input and output paths in the scripts to match your local file structure.
3. Run the outer-shell correction scripts for the relevant samples.
4. Run `r_celldensity_as_function_of_distance_blastema.m` to calculate cell density by distance from the kidney surface.
5. Use `r_verify_celldensity.m` for quality control and `r_plot_celldensity.m` to generate the final plots.

The main analysis requires sample-specific `.mat` files containing the reconstructed volumes and cell masks, as well as `diameters_nuclei.mat`.

## Output

The analysis generates cell-density measurements across distance bins, intermediate distance maps and projections, and plots used for the cell-density analyses in Figure 2.
