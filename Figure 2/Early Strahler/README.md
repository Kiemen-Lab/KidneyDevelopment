# Figure 2 – Early Strahler

This folder contains MATLAB scripts used to characterize the morphology of the E11 mouse and Macaque G40 kidney (including the branching structure and Strahler order of the developing ureteric tree), and to generate the corresponding visualizations and comparisons presented in Figure 2.

## Requirements

* MATLAB R2023b
* Image Processing Toolbox

## Scripts

| Script                              | Description                                                                 |
| ----------------------------------- | --------------------------------------------------------------------------- |
| `compute_Strahler.m`                | Computes Strahler order from the reconstructed branching tree.              |
| `countBranchGenerations.m`          | Quantifies branch generations within the reconstructed tree.                |
| `r_compute_info_E11.m`              | Calculates morphometric information for the E11 reconstruction.             |
| `r_compute_info_MA40.m`             | Calculates morphometric for the MA40 reconstruction.                        |
| `r_plot_comparison_metanephros.m`   | Generates comparisons of branching characteristics between reconstructions. |
| `r_save_view_renal_vesicles_E11.m`  | Saves views of the E11 renal vesicle/ureteric-tree reconstruction.          |
| `r_save_view_renal_vesicelesMA40.m` | Saves views of the MA40 renal vesicle/ureteric-tree reconstruction.         |
| `r_visualize_ureteric_trees.m`      | Generates 3D visualizations of the reconstructed ureteric trees.            |

## Usage

1. Ensure the required reconstructed data are available.
2. Update the input and output paths in the scripts to match your local file structure.
3. Run the relevant `r_compute_*` scripts to calculate branching information.
4. Use the plotting and visualization scripts to generate the corresponding analyses and figures.

## Output

The scripts generate Strahler orders, branch-generation measurements, quantitative comparisons, and 3D visualizations of the ureteric trees used in Figure 2.
