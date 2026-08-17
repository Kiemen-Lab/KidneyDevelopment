# Figure 1 – Late Kidney Development

This folder contains MATLAB scripts used to process and visualize the 3D kidney reconstructions from late developmental stages, including human, macaque, and mouse samples.

## Requirements

* MATLAB R2023b
* Image Processing Toolbox

## Scripts

| Script                            | Description                                                                         |
| --------------------------------- | ----------------------------------------------------------------------------------- |
| `r_clean_human_vol.m`             | Cleans and processes the reconstructed human kidney volume.                         |
| `r_clean_macaque_vol.m`           | Cleans and processes the reconstructed macaque kidney volume.                       |
| `r_clean_mouse_vol.m`             | Cleans and processes the reconstructed mouse kidney volume.                         |
| `r_save_inside.m`                 | Extracts and saves the internal region of a reconstructed volume.                   |
| `r_save_inside_human.m`           | Performs the corresponding internal-volume extraction for the human reconstruction. |
| `r_visualize_half_volume_human.m` | Generates a 3D visualization of half of the human reconstructed volume.             |
| `r_visualize_half_volume_mouse.m` | Generates a 3D visualization of half of the mouse reconstructed volume.             |

## Usage

1. Open the desired MATLAB script.
2. Update the input and output paths to match your local file structure.
3. Run the scripts in the order required for the analysis.
4. Processed volumes and 3D visualizations are saved to the specified output locations.

The scripts require the corresponding processed 3D reconstruction files as input. These data are not included in the repository and must be obtained separately.

## Output

The scripts generate cleaned 3D reconstruction volumes and visualizations used in the late developmental-stage analyses presented in Figure 1.
