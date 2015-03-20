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
#install.packages(c("rgdal","raster","lubridate","CircStats","Rcpp", "rgrass7", "optparse))

#load packages:
library(raster)  	 #Raster operation and I/O. Depends R (≥ 2.15.0)
library(rgdal)	     #Geospatial data abstraction library. Depends R (≥ 2.14.0)
library(lubridate)  #Make dealing with dates a little easier. Depends R (≥ 3.0.0)
library(CircStats)  #Circular Statistics - Von Mises distribution
library(Rcpp)       #Seamless R and C++ Integration. Depends R (≥ 3.0.0)
library(rgrass7)    #Interface Between GRASS 7 GIS and R. Depends R (≥ 2.12)
library(optparse)   #Parse args from command line

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
option_list = list(
  make_option("--file", action="store", default=NA, type='character', help="input raster file"),
  make_option(c("--start"), action="store", default=NA, type='integer', help="start year"),
  make_option(c("--end"), action="store", default=NA, type='integer', help="end year"),
  make_option(c("--seasonal"), action="store", default="YES", type='character', help="do you want the spread to be seasonal?"),
  make_option(c("--wind"), action="store", default="NO", type='character', help="do you want the spread to be affected by wind?"),
  make_option(c("--pwdir"), action="store", default=NA, type='character', help="predominant wind direction: N,NE,E,SE,S,SW,W,NW")
)

opt = parse_args(OptionParser(option_list=option_list))

#arguments <- commandArgs(TRUE)

##Input raster --> HOST INDEX
Nmax_rast <- readRAST(opt$file)
Nmax_rast <- raster(Nmax_rast)  #transform 'sp' obj to 'raster' obj
#Nmax_rast <- readRAST(arguments[1]) #in the current version this reads raster in GRASS as a 'sp' R object
#Nmax_rast <- raster(Nmax_rast)  #transform 'sp' obj to 'raster' obj
#Nmax_rast <- raster('./layers/HI_100m.img')

#raster resolution
res_win <- res(Nmax_rast)[1]

#clone Smax raster to I (=infected trees) raster and spores (=number of spores)
I_rast <- Nmax_rast 


##Start-End date: 
start <- opt$start
end <- opt$end
#start <- arguments[2]
#end <- arguments[3]

if (start > end) stop('start date must precede end date!!')

##Seasonality: Do you want the spread to be limited to certain months?
ss <- opt$seasonal   #'YES' or 'NO'
#ss <- arguments[4]   #'YES' or 'NO'
if (ss == 'YES') months_msk <- paste('0', 1:9, sep='') #1=January 9=September

##Wind: Do you want the spread to be affected by wind?
#wnd <- arguments[5]  'YES' or 'NO'
if (!(opt$wind %in% c('YES', 'NO'))) stop('You must specify whether you want spread by wind or not: use either YES or NO')


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

# create formatting expression for padding zeros depending on total number of steps
formatting_str = paste("%0", floor( log10( length(tstep) ) ) + 1, "d", sep='')
# grass date formatting
months_names = c('jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec')

##MAIN SIMULATION LOOP (weekly steps)
for (tt in tstep){
  split_date = unlist(strsplit(tt, '-'))
  
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
    I_rast_sp <- as(I_rast, 'SpatialGridDataFrame')
    writeRAST(I_rast_sp, vname=paste('infected_', sprintf(formatting_str, 0), sep=''), overwrite=TRUE) #write to GRASS raster file
	execGRASS('r.timestamp', map=paste('infected_', sprintf(formatting_str, cnt), sep=''), date=paste(split_date[3], months_names[as.numeric(split_date[2])], split_date[1]))
    
    #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='FLT4S', overwrite=TRUE) # % infected as output
    #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='INT1U', overwrite=TRUE) # nbr. infected hosts as output
    #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='LOG1S', overwrite=TRUE)  # 0=non infected 1=infected output
    
  }else{
    
    
    #check if there are any susceptible trees left on the landscape (IF NOT continue LOOP till the end)
    if(!any(susceptible > 0)) break
    
    #update week counter
    cnt <- cnt + 1
    
    #is current week time step within a spread month (as defined by input parameters)?
    if (ss == 'YES' & !any(substr(tt,6,7) %in% months_msk)) next
          
    #Total weather suitability:
    W <- as.matrix(Mstack[[cnt]] * Cstack[[cnt]])
    
    #GENERATE SPORES:  
    #integer matrix
    spores_mat <- SporeGenCpp(infected, W, rate = 4.4) 
    
    #SPORE DISPERSAL:  
    #'List'
    if (opt$wind == 'YES') {
      
      #wdir <- arguments[6] 
      #wdir <- 'N'      
      #Check if predominant wind direction has been specified correctly:
      if (!(opt$pwdir %in% c('N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'))) stop('A predominant wind direction must be specified: N, NE, E, SE, S, SW, W, NW')
      out <- SporeDispCppWind(spores_mat, S=susceptible, I=infected, W, rs=res_win, rtype='Cauchy', scale1=20.57, wdir=opt$pwdir, kappa=2)
    
    }else{
      out <- SporeDispCpp(spores_mat, S=susceptible, I=infected, W, rs=res_win, rtype='Cauchy', scale1=20.57)
    }  
    
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
    I_rast_sp <- as(I_rast, 'SpatialGridDataFrame')
    writeRAST(I_rast_sp, vname=paste('infected_', sprintf(formatting_str, cnt), sep=''), overwrite=TRUE) #write to GRASS raster file
	execGRASS('r.timestamp', map=paste('infected_', sprintf(formatting_str, cnt), sep=''), date=paste(split_date[3], months_names[as.numeric(split_date[2])], split_date[1]))
    
    #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='FLT4S', overwrite=TRUE) # % infected as output
    #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='INT1U', overwrite=TRUE) # nbr. infected hosts as output
    #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='LOG1S', overwrite=TRUE)  # 0=non infected 1=infected output
  }
  

}






