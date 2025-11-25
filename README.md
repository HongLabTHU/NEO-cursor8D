# NEO-cursor8D

Source code for the paper "Fine grained two-dimensional cursor control with epidural minimally invasive brain-computer interface". 

## System requirements
- Operating system: Windows 10 / 11 (recommended). MATLAB code can also run on Linux or macOS if toolboxes are available.  
- MATLAB: R2020b or later.  
- Recommended MATLAB Toolboxes: Signal Processing Toolbox, Statistics and Machine Learning Toolbox. Optional: Parallel Computing Toolbox, Optimization Toolbox.  
- Optional: Python 3.8+ for additional analysis/visualization (numpy, scipy, matplotlib, pandas).

## Installation & dependencies
1. Clone or copy the repository to your local machine.
2. MATLAB:
   - Open MATLAB and add the project folder to the MATLAB path:
       addpath(genpath('NEO-cursor8D'))
   - Ensure required toolboxes are installed (use Add-Ons or matlab.addons.installedAddons).
3. Python (optional):
   - Create a virtual environment and install common packages:
     - python -m venv .venv
     - .venv\Scripts\activate
     - pip install numpy scipy matplotlib pandas jupyter

## Quick start — entry points
1. Add project folder to MATLAB path (see above).
2. Typical run commands (replace with the actual entry script found):
     fig1s_single_movement
     fig2ac_single_movement
     ...


