# Figure 1 – Early Kidney Development

This folder contains MATLAB scripts used to visualize and annotate the reconstructed kidney and surrounding tissues during early kidney development (E11), as presented in Figure 1 of the manuscript.

## Requirements

* MATLAB R2023b
* Image Processing Toolbox
* 3D visualization functionality included with MATLAB

## Scripts

| Script                                  | Description                                                                                                    |
| --------------------------------------- | -------------------------------------------------------------------------------------------------------------- |
| `r_make_check_annotation.m`             | Creates an RGB overlay of an H&E image and a corresponding tissue annotation mask for visual quality control.  |
| `r_make_check_annotation_meso.m`        | Creates an overlay for annotation/quality control of the mesonephros.                                          |
| `r_make_check_annotation_meta.m`        | Creates an overlay for annotation/quality control of the metanephros.                                          |
| `r_makecmaps.m`                         | Generates color-map legends used to visualize the different anatomical structures in the reconstructed tissue. |
| `r_visualize_cleaned_entire_E11.m`      | Loads the cleaned E11 3D reconstruction and generates a 3D visualization of the major anatomical structures.   |
| `r_visualize_cleaned_mesonephros_E11.m` | Generates a 3D visualization of the reconstructed E11 mesonephros.                                             |
| `r_visualize_cleaned_metanephros_E11.m` | Generates a 3D visualization of the reconstructed E11 metanephros.                                             |

## Input data

The visualization scripts require the processed 3D reconstruction and/or annotated image data generated during the reconstruction workflow.

The scripts currently contain paths to the input files used during the study. Before running the scripts, replace these paths with the locations of the corresponding input files on your computer.

For example:

```matlab
load('path/to/E11_cleaned.mat')
```

The annotation scripts require the corresponding histological image and annotation mask as TIFF files.

## Usage

1. Install MATLAB R2023b and the required toolboxes.
2. Download or generate the required input data.
3. Open the desired MATLAB script.
4. Update the input file paths to point to your local data.
5. Run the script from MATLAB.

The `r_visualize_cleaned_*` scripts generate interactive 3D visualizations of the reconstructed tissue using MATLAB's 3D visualization tools.

## Expected output

The scripts generate:

* RGB overlays for annotation quality control.
* Color-map legends for anatomical structures.
* Interactive 3D visualizations of the E11 whole kidney, mesonephros, and metanephros.

These scripts were used to generate the visualizations presented in Figure 1 of the manuscript.
