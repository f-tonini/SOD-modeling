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

#install packages
#install.packages(c("rgdal","raster","lubridate","CircStats","Rcpp", "rgrass7"))

#load packages
library(raster)  		#Raster operation and I/O. Depends R (≥ 2.15.0)
library(rgdal)	    #Geospatial data abstraction library. Depends R (≥ 2.14.0)
library(lubridate)  #Make dealing with dates a little easier. Depends R (≥ 3.0.0)
library(CircStats)  #Circular Statistics - Von Mises distribution
library(Rcpp)       #Seamless R and C++ Integration. Depends R (≥ 3.0.0)
library(rgrass7)    #Interface Between GRASS 7 GIS and R. Depends R (≥ 2.12)


##Define the main working directory
#setwd("path_to_main_folder")

#Path to folders in which you want to save all your vector & raster files
#fOutput <- 'output'

##Create a physical copy of the subdirectory folder(s) where you save your output
##If the directory already exists it gives a warning, BUT we can suppress it using showWarnings = FALSE
#dir.create(fOutput, showWarnings = FALSE)

##Use an external source file w/ all modules (functions) used within this script. 
##Use FULL PATH if source file is not in the same folder w/ this script
source('./scripts/myfunctions_SOD.r')
sourceCpp("./scripts/myCppFunctions.cpp") #for C++ custom functions

#WEATHER SUITABILITY: read and stack weather suitability raster BEFORE running the simulation
lst <- dir('./layers/weather', pattern='\\.img$', full.names=T)
Mlst <- lst[grep("_m", lst)]
Clst <- lst[grep("_c", lst)]

Mstack <- stack(Mlst) #M = moisture; 
Cstack <- stack(Clst) #C = temperature;


###Input simulation parameters: #####

##Input raster --> HOST INDEX
Nmax_rast <- readRAST(commandArgs()[2]) #in the current version this reads raster in GRASS as a 'sp' R object
Nmax_rast <- raster(Nmax_rast)  #transform 'sp' obj to 'raster' obj
#Nmax_rast <- raster('./layers/HI_100m.img')

#raster resolution
res_win <- res(Nmax_rast)[1]

#clone Smax raster to I (=infected trees) raster and spores (=number of spores)
I_rast <- Nmax_rast 


##Start-End date: 
start <- commandArgs()[3]
end <- commandArgs()[4]
#start <- 2004
#end <- 2008

if (start > end) stop('start date must precede end date!!')

##Seasonality: Do you want the spread to be limited to certain months?
ss <- commandArgs()[5]   #'ON' or 'OFF'
#ss <- 'ON'

if (ss == 'ON') months_msk <- paste('0', 1:9, sep='') #1=January 9=September



set.seed(2000)


########################################################
##INITIAL SOURCES OF INFECTION:
#empty vector with counts of I (=infected trees)
Nmax <- Nmax_rast[]  #integer vector of Nmax
I_lst <- rep(0, length(Nmax))  #integer

#initial sources of infection (integer between 0 and length(Smax))
inf_src <- 2 #integer

#randomly sample the index of cells to be source of infections
#inf.index <- sample(which(HI.raster[] > 0), size = inf.sources)
inf_idx <- sample(which(Nmax > 0), size = inf_src)

#randomize the initial I (=infected trees) counts (this does NOT have to exceed Nmax)
I_lst[inf_idx] <- sapply(Nmax[inf_idx], FUN=function(x) sample(1:x, size=1)) 
###############################################################################


#Susceptibles = Nmax - Infected 
S_lst <- Nmax - I_lst   #integer vector

#integer matrix with susceptible and infected
susceptible <- matrix(S_lst, ncol=ncol(Nmax_rast), nrow=nrow(Nmax_rast), byrow=T)
infected <- matrix(I_lst, ncol=ncol(Nmax_rast), nrow=nrow(Nmax_rast), byrow=T)

cnt <- 0 #time counter to access raster stacks

#build time series:
dd_start <- as.POSIXlt(as.Date(paste(start,'-01-01',sep='')))
dd_end <- as.POSIXlt(as.Date(paste(end,'-12-31',sep='')))
tstep <- as.character(seq(dd_start, dd_end, 'weeks'))


##MAIN SIMULATION LOOP (weekly steps)
for (tt in tstep){
  
  if (tt == tstep[1]) {
    
    if(!any(S_lst > 0)) stop('Simulation ended. There are no more susceptible trees on the landscape!')
    
    ##CALCULATE OUTPUT TO PLOT:
    # 1) values as % infected
    I_rast[] <- ifelse(I_lst == 0, NA, round(I_lst/Nmax, 1))
    
    # 2) values as number of infected per cell
    #I_rast[] <- ifelse(I_lst == 0, NA, I_lst)    
    
    # 3) values as 0 (non infected) and 1 (infected) cell
    #I_rast[] <- ifelse(I_lst > 0, 1, 0) 
    #I_rast[] <- ifelse(I_lst > 0, 1, NA) 
    
    #breakpoints <- c(0, 0.25, 0.5, 0.75, 1)
    #colors <- c("yellow","gold","orange","red")
    #plot(I_rast, breaks=breakpoints, col=colors, main=tt)
    
    #WRITE TO FILE:
    I_rast <- as(I_rast, 'SpatialGridDataFrame')
    writeRAST(I_rast, vname='Infected_0.img', overwrite=TRUE) #write to GRASS raster file
    
    #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='FLT4S', overwrite=TRUE) # % infected as output
    #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='INT1U', overwrite=TRUE) # nbr. infected hosts as output
    #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='LOG1S', overwrite=TRUE)  # 0=non infected 1=infected output
    
  }else{
     
    #check if there are any susceptible trees left on the landscape (IF NOT continue LOOP till the end)
    if(!any(susceptible > 0)) break
    
    #is current week time step within a spread month (as defined by input parameters)?
    if (ss == 'ON' & !any(substr(tt,6,7) %in% months_msk)) next
    
    #update counter
    cnt <- cnt + 1
          
    #Total weather suitability:
    W <- as.matrix(Mstack[[cnt]] * Cstack[[cnt]])
    
    #GENERATE SPORES:  
    #integer matrix
    spores_mat <- SporeGenCpp(infected, W, rate = 4.4) 
    
    #SPORE DISPERSAL:  
    #'List'
    out <- SporeDispCpp(spores_mat, S=susceptible, I=infected, W, rs=res_win, rtype='Cauchy', wtype='Uniform', scale1=20.57)   #wtype='VM', wdir='N', scale1=20.57, kappa=2)
    
    #update R matrices
    susceptible <- out$S 
    infected <- out$I  
    
    #save matrix into raster obj
    I_rast[] <- infected
    
    ##CALCULATE OUTPUT TO PLOT:
    I_lst <- I_rast[]
    
    # 1) values as % infected
    I_rast[] <- ifelse(I_lst == 0, NA, round(I_lst/Nmax, 1))
    
    # 2) values as number of infected per cell
    #I_rast[] <- ifelse(I_lst == 0, NA, I_lst)    
    
    # 3) values as 0 (non infected) and 1 (infected) cell
    #I_rast[] <- ifelse(I_lst > 0, 1, 0) 
    #I_rast[] <- ifelse(I_lst > 0, 1, NA) 
    
    #breakpoints <- c(0, 0.25, 0.5, 0.75, 1)
    #colors <- c("yellow","gold","orange","red")
    #plot(I_rast, breaks=breakpoints, col=colors, main=tt)
    
    #WRITE TO FILE:
    I_rast <- as(I_rast, 'SpatialGridDataFrame')
    writeRAST(I_rast, vname=paste('Infected_', cnt, '.img', sep=''), overwrite=TRUE) #write to GRASS raster file
    
    #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='FLT4S', overwrite=TRUE) # % infected as output
    #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='INT1U', overwrite=TRUE) # nbr. infected hosts as output
    #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='LOG1S', overwrite=TRUE)  # 0=non infected 1=infected output
  }
  

}






