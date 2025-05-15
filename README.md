# CV_Matlab_Comsol

# Overview

This repository contains a collection of tools for electrochemical simulations, focusing on Cyclic Voltammetry (CV) and Electrochemical Impedance Spectroscopy (EIS). The toolkit includes MATLAB scripts for one electron step  CV simulations and a Jupyter notebook for EIS data fitting with custom circuit elements.

## Authors

* Ziwen Zhao [@ziwzh166](https://github.com/ziwzh166)
* Andrew J. Bangnall [@Swagmeister5000](https://github.com/Swagmeister5000)
* Nikolaos Kostopoulos [@ngkostop](https://github.com/ngkostop)
* Alina Sekretareva*[@alina-sekretareva](https://github.com/alina-sekretareva)

## Contents

### Cyclic Voltammetry (CV) Simulations

- [CV_Simulation_BV.m](CV_Simulation_BV.m): CV simulation using Butler-Volmer kinetics
- [CV_Simulation_MH.m](CV_Simulation_MH.m): CV simulation using Marcus-Hush theory
- [CV_SimulationAJBloop.m](CV_SimulationAJBloop.m): Parametric CV simulation with looping functionality
- [CV_comsol_Model1D.mph](CV_comsol_Model1D.mph): COMSOL Multiphysics model for 1D CV simulation based on Butler-Volmer kinetics

### Electrochemical Impedance Spectroscopy (EIS)

- [EISfittingCustomElements_AJB -NewAuto04-24c.ipynb](EISfittingCustomElements_AJB%20-NewAuto04-24c.ipynb): Jupyter notebook for fitting EIS data with custom circuit elements

### Documentation

- [CITATION.cff](CITATION.cff): Citation file format for referencing this work

## Getting Started

### Prerequisites

- MATLAB (for CV simulations)
- COMSOL Multiphysics (for the .mph model)
- Jupyter Notebook with Python (for EIS fitting)
- Required Python packages (listed in the Jupyter notebook)

### Running CV Simulations

1. Open MATLAB and navigate to the repository directory
2. Run one of the CV simulation scripts:

   ```matlab
   CV_Simulation_BV
   ```

   or

   ```matlab
   CV_Simulation_MH
   ```

### Using the EIS Fitting Tool

1. Prior to use please install the following modules:

   ```
   pip install impedance pandas matplotlib
   ```
2. Open the Jupyter notebook:

   ```
   jupyter notebook "EISfittingCustomElements_AJB -NewAuto04-24c.ipynb"
   ```
3. Follow the instructions within the notebook to fit your EIS data

If you use this software in your research, please cite it using the information provided in the [CITATION.cff](CITATION.cff) file.
