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
#install.packages(c("rgdal","raster","lubridate","CircStats","Rcpp", "rgrass7", "optparse", "plotrix"))

#load packages:
suppressPackageStartupMessages(library(raster))    #Raster operation and I/O. Depends R (≥ 2.15.0)
suppressPackageStartupMessages(library(rgdal))     #Geospatial data abstraction library. Depends R (≥ 2.14.0)
suppressPackageStartupMessages(library(lubridate)) #Make dealing with dates a little easier. Depends R (≥ 3.0.0)
suppressPackageStartupMessages(library(CircStats)) #Circular Statistics - Von Mises distribution
suppressPackageStartupMessages(library(Rcpp))      #Seamless R and C++ Integration. Depends R (≥ 3.0.0)
suppressPackageStartupMessages(library(rgrass7))   #Interface Between GRASS 7 GIS and R. Depends R (≥ 2.12)
suppressPackageStartupMessages(library(optparse))  #Parse args from command line
suppressPackageStartupMessages(library(plotrix))   #Add text annotations to plot

##Define the main working directory based on the current script path
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
base_name <- dirname(sub(file_arg, "", initial_options[grep(file_arg, initial_options)]))
setwd(paste(sep="/", base_name, ".."))

#Path to folders in which you want to save all your vector & raster files
#fOutput <- 'output'

##Create a physical copy of the subdirectory folder(s) where you save your output
##If the directory already exists it gives a warning, BUT we can suppress it using showWarnings = FALSE
#dir.create(fOutput, showWarnings = FALSE)

##Use an external source file w/ all modules (functions) used within this script. 
##Use FULL PATH if source file is not in the same folder w/ this script
source('./scripts/myfunctions_SOD.r')
sourceCpp("./scripts/myCppFunctions.cpp") #for C++ custom functions

###Input simulation parameters: #####
option_list = list(
  make_option(c("-hi","--host_index"), action="store", default=NA, type='character', help="input host index raster map"),
  make_option(c("-src","--sources"), action="store", default=NA, type='character', help="initial sources of infection raster map"),
  make_option(c("-img","--image"), action="store", default=NA, type='character', help="background satellite raster image for plotting"),
  make_option(c("-s","--start"), action="store", default=NA, type='integer', help="start year"),
  make_option(c("-e","--end"), action="store", default=NA, type='integer', help="end year"),
  make_option(c("-ss","--seasonal"), action="store", default="YES", type='character', help="seasonal spread?"),
  make_option(c("-w","--wind"), action="store", default="NO", type='character', help="spread using wind?"),
  make_option(c("-pd","--pwdir"), action="store", default=NA, type='character', help="predominant wind direction: N,NE,E,SE,S,SW,W,NW"),
  make_option(c("-o","--output"), action="store", default=NA, type='character', help="basename for output GRASS raster maps"),
  make_option(c("-n","--nth_output"), action="store", default=1, type='integer', help="output every nth map")
)

opt = parse_args(OptionParser(option_list=option_list))

##Input raster --> HOST INDEX
Nmax_rast <- readRAST(opt$host_index)
Nmax_rast <- raster(Nmax_rast)  #transform 'sp' obj to 'raster' obj
#Nmax_rast <- raster('./layers/HI_100m.img')

#raster resolution
res_win <- res(Nmax_rast)[1]


##Initial sources of infection:
I_rast <- readRAST(opt$sources)
I_rast <- raster(I_rast) #transform 'sp' obj to 'raster' obj
#I_rast <- raster('./layers/init_cnt.img')

#Susceptibles = Host Index (Nmax) - Infected (I_rast)
S_rast <- Nmax_rast - I_rast   

#integer matrix with susceptible and infected
susceptible <- as.matrix(S_rast)
infected <- as.matrix(I_rast)

##background satellite image for plotting
bkr_img <- raster(paste('./layers/', opt$image, sep='')) 

##Start-End date: 
start <- opt$start
end <- opt$end
#start <- arguments[2]
#end <- arguments[3]

if (start > end) stop('start date must precede end date!!')

#build time series for simulation steps:
dd_start <- as.POSIXlt(as.Date(paste(start,'-01-01',sep='')))
dd_end <- as.POSIXlt(as.Date(paste(end,'-12-31',sep='')))
tstep <- as.character(seq(dd_start, dd_end, 'weeks'))

#create formatting expression for padding zeros depending on total number of steps
formatting_str = paste("%0", floor( log10( length(tstep) ) ) + 1, "d", sep='')
#grass date formatting
months_names = c('jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec')


#WEATHER SUITABILITY: read and stack weather suitability raster BEFORE running the simulation
lst <- dir('./layers/weather', pattern='\\.img$', full.names=T)
Mlst <- lst[grep("_m", lst)]
Mlst <- grep(paste(as.character(seq(start,end)), collapse="|"), Mlst, value=TRUE)  #use only the raster files matching the years of interest
Clst <- lst[grep("_c", lst)]
Clst <- grep(paste(as.character(seq(start,end)), collapse="|"), Clst, value=TRUE) #use only the raster files matching the years of interest

Mstack <- stack(Mlst) #M = moisture; 
Cstack <- stack(Clst) #C = temperature;

##Seasonality: Do you want the spread to be limited to certain months?
ss <- opt$seasonal   #'YES' or 'NO'
if (ss == 'YES') months_msk <- paste('0', 1:9, sep='') #1=January 9=September

##Wind: Do you want the spread to be affected by wind?
if (!(opt$wind %in% c('YES', 'NO'))) stop('You must specify whether you want spread by wind or not: use either YES or NO')


#open window screen
windows(width = 10, height = 10, xpos=350, ypos=50, buffered = FALSE)
#quartz()  #use this on Mac OSX
#x11()     #use this on Linux (not tested!)

#plot background image
plot(bkr_img)

#plot coordinates for plotting text:
xpos <- (bbox(I_rast)[1,2] + bbox(I_rast)[1,1]) / 2
ypos <- bbox(I_rast)[2,2] - 150

#time counter to access pos index in weather raster stacks
cnt <- 0 

## ----> MAIN SIMULATION LOOP (weekly time steps) <------
for (tt in tstep){
  
  #split date string for raster time stamp
  split_date = unlist(strsplit(tt, '-'))
  
  if (tt == tstep[1]) {
    
    if(!any(susceptible > 0)) stop('Simulation ended. There are no more susceptible trees on the landscape!')
    
    ##CALCULATE OUTPUT TO PLOT:
    # 1) values as % infected
    I_rast[] <- ifelse(I_rast[] == 0, NA, round(I_rast[]/Nmax_rast[], 1))
    
    # 2) values as number of infected per cell
    #I_rast[] <- ifelse(I_rast[] == 0, NA, I_rast[])
    
    # 3) values as 0 (non infected) and 1 (infected) cell
    #I_rast[] <- ifelse(I_rast[] > 0, 1, 0) 
    #I_rast[] <- ifelse(I_rast[] > 0, 1, NA) 
    
    #PLOT: overlay current plot on background image
    bks <- c(0, 0.25, 0.5, 0.75, 1)
    #colors <- c("yellow","gold","orange","red")
    image(I_rast, breaks=bks, col=rev(heat.colors(length(bks)-1, alpha=1)), add=T, axes=F, box=F, ann=F, legend=F, useRaster=T)
    boxed.labels(xpos, ypos, tt, bg="white", border=NA, font=2)
    
    #WRITE TO FILE:
    I_rast_sp <- as(I_rast, 'SpatialGridDataFrame')
    writeRAST(I_rast_sp, vname=paste(opt$output, '_', sprintf(formatting_str, 0), sep=''), overwrite=TRUE) #write to GRASS raster file
	  execGRASS('r.timestamp', map=paste(opt$output, '_', sprintf(formatting_str, cnt), sep=''), date=paste(split_date[3], months_names[as.numeric(split_date[2])], split_date[1]))
    
    #writeRaster(I_rast, filename=paste('./', fOutput, '/', opt$output, '_', sprintf(formatting_str, 0), sep=''), format='HFA', datatype='FLT4S', overwrite=TRUE) # % infected as output
    #writeRaster(I_rast, filename=paste('./', fOutput, '/', opt$output, '_', sprintf(formatting_str, 0), sep=''), format='HFA', datatype='INT1U', overwrite=TRUE) # nbr. infected hosts as output
    #writeRaster(I_rast, filename=paste('./', fOutput, '/', opt$output, '_', sprintf(formatting_str, 0), sep=''), format='HFA', datatype='LOG1S', overwrite=TRUE)  # 0=non infected 1=infected output
    
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
    # 1) values as % infected
    I_rast[] <- ifelse(I_rast[] == 0, NA, round(I_rast[]/Nmax_rast[], 1))
    
    # 2) values as number of infected per cell
    #I_rast[] <- ifelse(I_rast[] == 0, NA, I_rast[])
    
    # 3) values as 0 (non infected) and 1 (infected) cell
    #I_rast[] <- ifelse(I_rast[] > 0, 1, 0) 
    #I_rast[] <- ifelse(I_rast[] > 0, 1, NA) 
        
    if (cnt %% opt$nth_output == 0){
      
      #PLOT: overlay current plot on background image
      bks <- c(0, 0.25, 0.5, 0.75, 1)
      #colors <- c("yellow","gold","orange","red")
      image(I_rast, breaks=bks, col=rev(heat.colors(length(bks)-1, alpha=.5)), add=T, axes=F, box=F, ann=F, legend=F, useRaster=T)
      boxed.labels(xpos, ypos, tt, bg="white", border=NA, font=2)
      
      #WRITE TO FILE:
      I_rast_sp <- as(I_rast, 'SpatialGridDataFrame')
      writeRAST(I_rast_sp, vname=paste(opt$output, '_', sprintf(formatting_str, cnt), sep=''), overwrite=TRUE) #write to GRASS raster file
      execGRASS('r.timestamp', map=paste(opt$output, '_', sprintf(formatting_str, cnt), sep=''), date=paste(split_date[3], months_names[as.numeric(split_date[2])], split_date[1]))  
      
      #writeRaster(I_rast, filename=paste('./', fOutput, '/', opt$output, '_', sprintf(formatting_str, cnt), sep=''), format='HFA', datatype='FLT4S', overwrite=TRUE) # % infected as output
      #writeRaster(I_rast, filename=paste('./', fOutput, '/', opt$output, '_', sprintf(formatting_str, cnt), sep=''), format='HFA', datatype='INT1U', overwrite=TRUE) # nbr. infected hosts as output
      #writeRaster(I_rast, filename=paste('./', fOutput, '/', opt$output, '_', sprintf(formatting_str, cnt), sep=''), format='HFA', datatype='LOG1S', overwrite=TRUE)  # 0=non infected 1=infected output
    }
    
  }
  
}

message("Spread model finished")





