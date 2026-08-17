# Figure 2 – Cell Composition

This folder contains MATLAB scripts used to calculate and visualize cellular composition and cell fractions across the reconstructed kidney samples, as presented in Figure 2.

## Requirements

* MATLAB R2023b
* Statistics and Machine Learning Toolbox

## Scripts

| Script                      | Description                                                            |
| --------------------------- | ---------------------------------------------------------------------- |
| `r_save_total_vols_cells.m` | Calculates total tissue and cellular volumes for each sample.          |
| `r_save_cellfraction.m`     | Calculates the fraction of the reconstructed tissue occupied by cells. |
| `r_save_compositions.m`     | Calculates cellular composition across the reconstructed samples.      |
| `r_plot_cellfrac.m`         | Generates plots of cellular fractions across samples.                  |
| `r_plot_compositions.m`     | Generates plots of cellular composition.                               |
| `r_plot_celldensity_each.m` | Generates cell-density plots for individual samples.                   |
| `create_stacked_plot.m`     | Creates stacked plots summarizing cellular composition.                |

## Input data

The scripts require the processed 3D reconstruction data and cell segmentation/annotation data generated during the reconstruction workflow. The folder also contains `diameters_nuclei.mat`, which provides the manually measured nuclear diameter information used in the analysis.

Before running the scripts, update the input and output paths to match your local file structure.

## Usage

1. Ensure the required reconstructed volumes and cell annotation data are available.
2. Update the file paths in the scripts.
3. Run the `r_save_*` scripts to calculate the cell fractions and compositions.
4. Run the `r_plot_*` scripts to generate the corresponding visualizations.

## Output

The scripts generate cell fractions, cellular composition measurements, and plots used for the composition analyses presented in Figure 2.
