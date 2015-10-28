#---------------------------------------------------------------------------------------------------------------
# Name:         SOD_aniso_clim.r
# Purpose:      Lattice-based simulation of the climate-driven anisotropic spread of pathogen P. ramorum over a heterogeneous landscape.
# Author:       Francesco Tonini
# Email:        ftonini84@gmail.com
# Created:      01/07/2015
# Copyright:    (c) 2015 by Francesco Tonini
# License:      GNU General Public License (GPL)
# Software:     Tested successfully using R version 3.0.2 (http://www.r-project.org/)
#-----------------------------------------------------------------------------------------------------------------------

#load packages:
suppressPackageStartupMessages(library(raster))    #Raster operation and I/O. Depends R (≥ 2.15.0)
suppressPackageStartupMessages(library(rgdal))     #Geospatial data abstraction library. Depends R (≥ 2.14.0)
suppressPackageStartupMessages(library(lubridate)) #Make dealing with dates a little easier. Depends R (≥ 3.0.0)
suppressPackageStartupMessages(library(CircStats)) #Circular Statistics - Von Mises distribution
suppressPackageStartupMessages(library(Rcpp))      #Seamless R and C++ Integration. Depends R (≥ 3.0.0)
suppressPackageStartupMessages(library(plotrix))   #Add text annotations to plot
suppressPackageStartupMessages(library(ncdf))   #work with NetCDF datasets

##Define the main working directory based on the current script path
setwd("D:\\SOD-modeling")

##Use an external source file w/ all modules (functions) used within this script. 
##Use FULL PATH if source file is not in the same folder w/ this script
source('./scripts/myfunctions_SOD.r')
sourceCpp("./scripts/myCppFunctions.cpp") #for C++ custom functions

##Input rasters: abundance (tree density per hectare)
#----> UMCA
umca_rast <- raster("./layers/UMCA_den_100m.img")
#----> ALL SOD-affected oaks
oaks_rast <- raster("./layers/OAKS_den_100m.img")
#max
mx <- cellStats(oaks_rast, stat='max') 
#----> All live trees
lvtree_rast <- raster("./layers/TPH_den_100m.img")

#raster resolution
res_win <- res(umca_rast)[1]

###################################
### INFECTED AND SUSCEPTIBLES ####

##Initial infection (OAKS):
I_oaks_rast <- raster("./layers/init_2000_cnt.img") 

#define matrices for infected and susceptible species of interest
I_oaks <- as.matrix(I_oaks_rast)
S_oaks <- as.matrix(oaks_rast - I_oaks_rast)
I_umca <- matrix(0, nrow=res_win, ncol=res_win)
S_umca <- as.matrix(umca_rast)

##Initialize infected trees for each species (!!NEEDED UNLESS EMPIRICAL INFO IS AVAILABLE!!)
if(any(S_umca[I_oaks > 0] > 0)) I_umca[I_oaks > 0] <- mapply(function(x,y) ifelse(x > y, min(c(x,y*2)), x), 
                                                             S_umca[I_oaks > 0], I_oaks[I_oaks > 0]) 
##update susceptible matrices by subtracting the initialized infections 
S_umca <- S_umca - I_umca 

##define matrix for immune live trees
N_live <- as.matrix(lvtree_rast)

##background satellite image for plotting
bkr_img <- raster("./layers/ortho_5m_color.tif") 

##Start-End date: 
start <- 2000
end <- 2010

if (start > end) stop('start date must precede end date!!')

#build time series for simulation steps:
dd_start <- as.POSIXlt(as.Date(paste(start,'-01-01',sep='')))
dd_end <- as.POSIXlt(as.Date(paste(end,'-12-31',sep='')))
tstep <- as.character(seq(dd_start, dd_end, 'weeks'))

#create formatting expression for padding zeros depending on total number of steps
formatting_str = paste("%0", floor( log10( length(tstep) ) ) + 1, "d", sep='')

##WEATHER SUITABILITY: read and stack weather suitability raster BEFORE running the simulation

#weather coefficients
mcf.array <- get.var.ncdf(open.ncdf('./layers/weather/weatherCoeff_2000_2014.nc'),  varid = "Mcoef") #M = moisture;
ccf.array <- get.var.ncdf(open.ncdf('./layers/weather/weatherCoeff_2000_2014.nc'),  varid = "Ccoef") #C = temperature;

##Seasonality: Do you want the spread to be limited to certain months?
ss <- 'YES'   #'YES' or 'NO'
if (ss == 'YES') months_msk <- paste('0', 1:9, sep='') #1=January 9=September

##Wind: Do you want the spread to be affected by wind?
wind <- 'YES' #'YES' or 'NO'
pwdir <- 'NE'
spore_rate <- 4.4
nth_output <- 4

#plot background image
plot(bkr_img, xaxs = "i", yaxs = "i")

#plot coordinates for plotting text:
xpos <- (bbox(umca_rast)[1,2] + bbox(umca_rast)[1,1]) / 2
ypos <- bbox(umca_rast)[2,2] - 150

#time counter to access pos index in weather raster stacks
cnt <- 0 

## ----> MAIN SIMULATION LOOP (weekly time steps) <------
for (tt in tstep){
  
  #split date string for raster time stamp
  split_date = unlist(strsplit(tt, '-'))
  
  if (tt == tstep[1]) {
    
    if(!any(S_oaks > 0)) stop('Simulation ended. All oaks are infected!')
    
    ##CALCULATE OUTPUT TO PLOT: 
    # 1) values as % infected
    #I_oaks_rast[] <- ifelse(I_oaks_rast[] == 0, NA, I_oaks_rast[]/oaks_rast[])
    
    # 2) values as number of infected per cell
    I_oaks_rast[] <- ifelse(I_oaks_rast[] == 0, NA, I_oaks_rast[])
    
    # 3) values as 0 (non infected) and 1 (infected) cell
    #I_oaks_rast[] <- ifelse(I_oaks_rast[] > 0, 1, 0) 
    #I_oaks_rast[] <- ifelse(I_oaks_rast[] > 0, 1, NA) 
    
    #PLOT: overlay current plot on background image
    
    #bks <- c(0, 0.25, 0.5, 0.75, 1)
    #my_palette <- colorRampPalette(c("springgreen", "yellow1", "orange", "red1"))(n = 4)
    #image(I_oaks_rast, breaks=bks, col=addalpha(my_palette, 1), add=T, axes=F, box=F, ann=F, legend=F, useRaster=T)
    bks <- seq(0, mx, length = 10)
    image(I_oaks_rast, breaks=bks, col=rev(heat.colors(length(bks)-1, alpha=1)), add=T, axes=F, box=F, ann=F, legend=F, useRaster=T)
    boxed.labels(xpos, ypos, tt, bg="white", border=NA, font=2)
    
    #WRITE TO FILE:
    #writeRaster(I_oaks_rast, filename=paste('./', fOutput, '/', opt$output, '_', sprintf(formatting_str, cnt), sep=''), format='HFA', datatype='FLT4S', overwrite=TRUE) # % infected as output
    #writeRaster(I_oaks_rast, filename=paste('./', fOutput, '/', opt$output, '_', sprintf(formatting_str, cnt), sep=''), format='HFA', datatype='INT1U', overwrite=TRUE) # nbr. infected hosts as output
    #writeRaster(I_oaks_rast, filename=paste('./', fOutput, '/', opt$output, '_', sprintf(formatting_str, cnt), sep=''), format='HFA', datatype='LOG1S', overwrite=TRUE)  # 0=non infected 1=infected output
    
  }else{
    
    #check if there are any susceptible oaks left on the landscape (IF NOT continue LOOP till the end)
    if(!any(S_oaks > 0)) break
    
    #update week counter
    cnt <- cnt + 1
    
    #is current week time step within a spread month (as defined by input parameters)?
    if (ss == 'YES' & !any(substr(tt,6,7) %in% months_msk)) next
    
    #Total weather suitability:
    W <- mcf.array[,,cnt] * ccf.array[,,cnt]
    
    #GENERATE SPORES:  
    #integer matrix
    set.seed(42)
    spores_mat <- SporeGenCpp(I_umca, W, rate = spore_rate) #rate: spores/week for each infected host (4.4 default)
    
    #SPORE DISPERSAL:  
    #'List'
    if (wind == 'YES') {
      
      #Check if predominant wind direction has been specified correctly:
      if (!(pwdir %in% c('N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'))) stop('A predominant wind direction must be specified: N, NE, E, SE, S, SW, W, NW')
      out <- SporeDispCppWind_mh(spores_mat, S_UM=S_umca, S_OK=S_oaks, I_UM=I_umca, I_OK=I_oaks, N_LVE=N_live, 
                                 W, rs=res_win, rtype='Cauchy', scale1=20.57, wdir=pwdir, kappa=2)
    
    }else{
      out <- SporeDispCpp_mh(spores_mat, S_UM=S_umca, S_OK=S_oaks, I_UM=I_umca, I_OK=I_oaks, N_LVE=N_live,
                             W, rs=res_win, rtype='Cauchy', scale1=20.57) ##TO DO
    }  
    
    #update R matrices:
    #UMCA
    S_umca <- out$S_UM 
    I_umca <- out$I_UM 
    #oaks
    S_oaks <- out$S_OK 
    I_oaks <- out$I_OK
    
    ##CALCULATE OUTPUT TO PLOT:
    I_oaks_rast[] <- I_oaks
    
    # 1) values as % infected
    #I_oaks_rast[] <- ifelse(I_oaks_rast[] == 0, NA, I_oaks_rast[]/oaks_rast[])
    
    # 2) values as number of infected per cell
    I_oaks_rast[] <- ifelse(I_oaks_rast[] == 0, NA, I_oaks_rast[])
    
    # 3) values as 0 (non infected) and 1 (infected) cell
    #I_oaks_rast[] <- ifelse(I_oaks_rast[] > 0, 1, 0) 
    #I_oaks_rast[] <- ifelse(I_oaks_rast[] > 0, 1, NA) 
        
    if (cnt %% nth_output == 0){
      
      #PLOT: overlay current plot on background image
      #bks <- c(0, 0.25, 0.5, 0.75, 1)
      #my_palette <- colorRampPalette(c("springgreen", "yellow1", "orange", "red1"))(n = 4)
      #image(I_oaks_rast, breaks=bks, col=addalpha(my_palette, .5), add=T, axes=F, box=F, ann=F, legend=F, useRaster=T)
      bks <- seq(0, mx, length = 10)
      image(I_oaks_rast, breaks=bks, col=rev(heat.colors(length(bks)-1, alpha=1)), add=T, axes=F, box=F, ann=F, legend=F, useRaster=T)
      boxed.labels(xpos, ypos, tt, bg="white", border=NA, font=2)
      
      #WRITE TO FILE:
      #writeRaster(I_oaks_rast, filename=paste('./', fOutput, '/', opt$output, '_', sprintf(formatting_str, cnt), sep=''), format='HFA', datatype='FLT4S', overwrite=TRUE) # % infected as output
      #writeRaster(I_oaks_rast, filename=paste('./', fOutput, '/', opt$output, '_', sprintf(formatting_str, cnt), sep=''), format='HFA', datatype='INT1U', overwrite=TRUE) # nbr. infected hosts as output
      #writeRaster(I_oaks_rast, filename=paste('./', fOutput, '/', opt$output, '_', sprintf(formatting_str, cnt), sep=''), format='HFA', datatype='LOG1S', overwrite=TRUE)  # 0=non infected 1=infected output
      
      
    }
    
  }
  
}

message("Spread model finished")





