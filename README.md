# High-precision-3D-grain-boundary-modeling-and-mechanical-simulation-of-S2-columnar-ice
An open-source toolset and compiled 4D-LSM executable for the high-precision 3D grain boundary modeling and micromechanical simulation of S2 columnar ice.
# High-Precision 3D Grain Boundary Modeling and Mechanical Simulation of S2 Columnar Ice

## Overview
This repository provides the open-source Python modeling scripts and executable simulation cases for our research on the micromechanical behavior of S2 columnar ice. It includes the reproducible workflow from 3D Voronoi grain geometry generation to macroscopic mechanical testing using the compiled Four-Dimensional Lattice Spring Model (4D-LSM).

## Repository Structure
This repository contains two main Python modeling scripts and four distinct mechanical simulation cases:

### 1. Python Modeling Scripts
*(Note: You can find the Python scripts in the root directory)*
* **`BooleanOperation.py`** : Generates the initial 3D Voronoi tessellation to simulate the grain structure of S2 columnar ice, allowing precise control over geometric imperfections such as grain growth taper angles and roundness.
* **`GenerateColumnar.py`**: Performs complex geometric Boolean operations on fully closed solid STL structures to obtain the specific specimen shapes required for various mechanical tests.

### 2. Simulation Cases (4D-LSM)
We provide the full configuration and input files (`BALL3D_PSLICE.dat`, `MLSJDat.dat`, etc.) for four classical mechanical tests. Researchers can directly run the 4D-LSM executable within these folders to reproduce our findings:
* **`UniaxialCompressionTest/`**: Simulation of uniaxial compression to study the basic compressive strength and intergranular delamination.
* **`BiaxialCompressionTest/`**: Biaxial compression simulation to investigate the effects of lateral confinement on the ice strength.
* **`TrueTriaxialCompressionTest/`**: True triaxial compression simulation designed to evaluate the material under complex 3D stress states.
* **`BrazilSplitTest/`**: Brazilian splitting test featuring complex arc-shaped boundaries to evaluate the tensile behavior of the columnar ice discs.

## How to Use
1. **Geometry Generation**: Run the Python scripts to generate the initial structural files for the ice grains.
2. **Run Simulation**: Navigate to any of the four simulation case folders. Execute the provided 4D-LSM program alongside the input `.dat` files.
3. **Data Output**: The program computes the inter-particle interaction forces and equivalent stress tensors, outputting the macroscopic mechanical responses (e.g., stress-strain curves) and failure states of the discontinuous joints.

## License
This project is open-sourced under the **GNU General Public License v3.0 (GPL-3.0)**. Anyone who modifies or distributes this code must also release their derived work under the same open-source license.

## Citation
If you find our modeling methodology or simulation tools helpful for your research, please consider citing our paper:
> *[Authors], "High-precision 3D grain boundary modeling and mechanical simulation of S2 columnar ice", [Journal Name], [Year]. DOI: [Your DOI here]*
