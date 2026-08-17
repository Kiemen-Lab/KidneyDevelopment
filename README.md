
# Overview

This repository contains the code used to generate the analyses and figures presented in the manuscript:

"Cross-species Whole-organ Reconstruction Reveals Conserved Cellular Programs and Divergent Architecture During Kidney Development."

The repository is organized into folders corresponding to the main figures in the manuscript (Figure1-Figure5). Each figure folder contains the scripts used to perform the corresponding analyses and generate the associated figures. 
The DeepLearning folder contains the scripts and models used for tissue segmentation and deep learning-based model inference. 

The repository is intended to facilitate reproduction of the analyses presented in the manuscript. See Installation, Repository Structure, and Reproducing the Analyses sections below for detailed instructions.

# System requirements
Software
MATLAB R2023b (tested)
R Studio 

The code has been tested on:

Windows 10
Windows 11

Most analyses can be run on a standard desktop computer (≥16 GB RAM recommended). An NVIDIA GPU with CUDA support is recommended, but not required, for deep learning inference.

# MATLAB Requirements
The MATLAB-based analyses were developed and tested using MATLAB R2023b and require the following toolboxes: 
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox
- Deep Learning Toolbox 

# Instructions for use
Each figure folder (Figure1-Figure5) contains the scripts used to perform the analyses and generate the corresponding figures in the manuscript.

For each analysis: 
1. Review the documentation provided in the corresponding figure folder.
2. Ensure the required input data and software dependencies are available 
3. Update the input and output directory paths in the relevant scripts to point to the location of your data.
4. Run the scripts in the order specified in the folder documentation.

The scripts can be adapted to analyze user-provided datasets with the appropriate input format. See the documentation within each figure folder for information on required inputs, expected file formats, script functions, and example commands. 
# Data availability

The datasets supporting this study will be made publicly available upon publication of the manuscript.

# License

This repository is distributed under the MIT License. See the LICENSE file for details.

# Citation

If you use this code, please cite the associated manuscript. Citation details will be updated upon publication.
