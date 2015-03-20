# SOD-modeling
This repository contains the scripts used to develop a stochastic landscape spread model of forest pathogen *P. ramorum*, causal agent of the emerging infectious disease sudden oak death (SOD). This project is adjusted from the published research article:

Ross K. Meentemeyer, Nik J. Cunniffe, Alex R. Cook, Joao A. N. Filipe, Richard D. Hunter, David M. Rizzo, and Christopher A. Gilligan 2011. Epidemiological modeling of invasion in heterogeneous landscapes: spread of sudden oak death in California (1990–2030). *Ecosphere* 2:art17. [http://dx.doi.org/10.1890/ES10-00192.1] (http://www.esajournals.org/doi/abs/10.1890/ES10-00192.1) 

### *layers*
This folder contains all the GIS layers necessary to test run our code. Subfolder **_weather_** contains weekly **m** (=moisture) and **c** (=mean temperature) raster files. These are called by the main R script during execution.

### *scripts*
This folder contains the scripts used in this project. Files with a .r extension are written in R, while a .cpp extension indicates C++ files written using the R package [Rcpp] (http://www.rcpp.org/). The script **SOD\_aniso\_clim.r** implements the spread of *P. ramorum* with the option of accounting for wind direction (anisotropy) and weather suitability.

## Usage

As of now, the main R script is meant to be run from within [**GRASS GIS 7.0.0**] (http://grass.osgeo.org/) using the [rgrass7] (http://cran.r-project.org/web/packages/rgrass7/index.html) package. If you want to use GRASS GIS 6 instead, use the [spgrass6] (http://www.cran.r-project.org/web/packages/spgrass6/index.html) package. However, the latter may require some manual adjustments to the code as functions changed names between versions. 

**NOTE**: **_if you plan on running the code within R only, please checkout the main script and replace the appropriate I/O parameters to read/write raster files and set up initial parameters._**

Before running the code make sure to follow these steps:

1. Install R (our code was tested with version 3.0.2) from [here] (http://www.r-project.org/)

2. Open R and install the required packages using the statement:
**install.packages(c("rgdal","raster","lubridate","CircStats","Rcpp","rgrass7", "optparse"))**

3. Open GRASS GIS 7 and set up a location where to store/import your initial raster file. This step is necessary in order to (a) have a consistent projection across multiple files (all files within a GRASS location should be in the same coordinate system), and (b) make sure all the output raster files will be stores in the same GRASS location
4. Set up the GRASS region to match the imported initial raster file using:
  * **g.region raster = _your\_initial\_raster\_name_**

5. Run the code by invoking the R script from the GRASS terminal prompt using:
    * **Rscript _path\_to\_script_ --_arguments\_list_** 
    
Each parameter in the argument list is specified using either the *-shortflag* or the *--longflag* option.

**-f** or **--file** _your input raster_

**-s** or **--start** _starting year for the simulation_

**-e** or **--end** _last year for the simulation_

**-ss** or **--seasonal** _restrict pathogen spread to the rainy season only YES (default), NO_

**-w** or **--wind** _account for wind direction for the pathogen spread YES, NO (default)_

**-pd** or **--pwdir** _if you account for wind, what should the predominant wind direction be for the area? N (= North), NE (= Northeast), E (= East), etc._
