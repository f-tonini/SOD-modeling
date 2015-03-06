# SOD-modeling
This repository contains the scripts used to develop a stochastic landscape spread model of forest pathogen *P. ramorum*, causal agent of the emerging infectious disease sudden oak death (SOD). This project is adjusted from the published research article:

Ross K. Meentemeyer, Nik J. Cunniffe, Alex R. Cook, Joao A. N. Filipe, Richard D. Hunter, David M. Rizzo, and Christopher A. Gilligan 2011. Epidemiological modeling of invasion in heterogeneous landscapes: spread of sudden oak death in California (1990–2030). *Ecosphere* 2:art17. [http://dx.doi.org/10.1890/ES10-00192.1] (http://www.esajournals.org/doi/abs/10.1890/ES10-00192.1) 

### *layers*
This folder contains all the GIS layers necessary to test run our code. Subfolders **M** (moisture) and **C** (mean temperature) contain historical weekly raster files and are called by the main R script during execution.

### *scripts*
This folder contains the scripts used in this project. Files with a .r extension are written in R, while a .cpp extension indicates C++ files written using the R package [Rcpp] (http://www.rcpp.org/). The script **SOD\_aniso\_clim.r** implements the spread of *P. ramorum* with the option of accounting for wind direction (anisotropy) and weather suitability.

## Usage

Before running the code make sure to follow these steps:

1. Install R (our code was tested with version 3.0.2) from [here] (http://www.r-project.org/)
2. Open R and install the required packages using the statement:
**install.packages(c("rgdal","raster","lubridate","CircStats","Rcpp"))**
3. Set your working directory (SOD-modeling should be your main folder)
4. Run the script

