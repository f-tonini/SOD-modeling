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
#source('./scripts/myfunctions_SOD.r')
sourceCpp("./scripts/myCppFunctions.cpp") #for C++ custom functions

###Input simulation parameters: #####
option_list = list(
  make_option(c("-u","--umca"), action="store", default=NA, type='character', help="input bay laurel (UMCA) raster map"),
  make_option(c("-qg","--quag"), action="store", default=NA, type='character', help="input coast live oak (QUAG) raster map"),
  make_option(c("-qk","--quke"), action="store", default=NA, type='character', help="input black oak (QUKE) raster map"),
  make_option(c("-lvt","--livetree"), action="store", default=NA, type='character', help="input live tree (all) raster map"),
  make_option(c("-src","--sources"), action="store", default=NA, type='character', help="initial sources of infection raster map"),
  make_option(c("-img","--image"), action="store", default=NA, type='character', help="background satellite raster image for plotting"),
  make_option(c("-s","--start"), action="store", default=NA, type='integer', help="start year"),
  make_option(c("-e","--end"), action="store", default=NA, type='integer', help="end year"),
  make_option(c("-ss","--seasonal"), action="store", default="YES", type='character', help="seasonal spread?"),
  make_option(c("-w","--wind"), action="store", default="NO", type='character', help="spread using wind?"),
  make_option(c("-pd","--pwdir"), action="store", default=NA, type='character', help="predominant wind direction: N,NE,E,SE,S,SW,W,NW"),
  make_option(c("-spr","--spore_rate"), action="store", default=4.4, type='numeric', help="spore production rate per week for each infected tree"),
  make_option(c("-o","--output"), action="store", default=NA, type='character', help="basename for output GRASS raster maps"),
  make_option(c("-n","--nth_output"), action="store", default=1, type='integer', help="output every nth map"),
  make_option(c("-scn","--scenario"), action="store", default=NA, type='character', help="future weather scenario")
)

opt = parse_args(OptionParser(option_list=option_list))

##Input rasters: abundance (tree density per hectare)
#----> UMCA
umca_rast <- readRAST(opt$umca)
umca_rast <- round(raster(umca_rast))  #transform 'sp' obj to 'raster' obj
#----> QUAG
quag_rast <- readRAST(opt$quag)
quag_rast <- round(raster(quag_rast))  
#----> QUAG
quke_rast <- readRAST(opt$quke)
quke_rast <- round(raster(quke_rast))
#----> All live trees
livetree_rast <- readRAST(opt$livetree)
livetree_rast <- round(raster(livetree_rast))
livetree_rast <- livetree_rast - (umca_rast + quag_rast + quke_rast)

#overall oak host abundance raster
oaks_rast <- quag_rast + quke_rast

#raster resolution
res_win <- res(umca_rast)[1]

##Initial sources of infection (UMCA):
I_umca_rast <- readRAST(opt$sources)
I_umca_rast <- raster(I_umca_rast) 

#initialize zero-infection rasters for oaks:
I_quag_rast <- I_umca_rast  ##clone raster
I_quag_rast[] <- 0          ##init to 0 
I_quke_rast <- I_umca_rast  ##clone raster
I_quke_rast[] <- 0          ##init to 0

#overall oak infection raster:
I_oaks_rast <- I_quag_rast + I_quke_rast

#Susceptibles = Current Abundance - Infected (I_rast)
S_umca_rast <- umca_rast - I_umca_rast

#integer matrix with susceptible and infected
susceptible_umca <- as.matrix(S_umca_rast)
infected_umca <- as.matrix(I_umca_rast)
susceptible_quag <- as.matrix(quag_rast)
infected_quag <- as.matrix(I_quag_rast)
susceptible_quke <- as.matrix(quke_rast)
infected_quke <- as.matrix(I_quke_rast)
#live tree (all) matrix
livetree_matr <- as.matrix(livetree_rast)

##background satellite image for plotting
bkr_img <- raster(paste('./layers/', opt$image, sep='')) 

##Start-End date: 
start <- opt$start
end <- opt$end

if (start > end) stop('start date must precede end date!!')

#build time series for simulation steps:
dd_start <- as.POSIXlt(as.Date(paste(start,'-01-01',sep='')))
dd_end <- as.POSIXlt(as.Date(paste(end,'-12-31',sep='')))
tstep <- as.character(seq(dd_start, dd_end, 'weeks'))

#create formatting expression for padding zeros depending on total number of steps
formatting_str = paste("%0", floor( log10( length(tstep) ) ) + 1, "d", sep='')
#grass date formatting
months_names = c('jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec')

##WEATHER SUITABILITY: read and stack weather suitability raster BEFORE running the simulation
#list of ALL weather layers
lst <- dir('./layers/weather', pattern='\\.img$', full.names=T)

#strip and read the last available historical year
last_yr <- as.numeric(unlist(strsplit(basename(tail(lst, n=1)),'_'))[1])
if (start > last_yr) stop('start simulation date needs to be within the range of available historical data')

#sublist of weather coefficients
Mlst <- lst[grep("_m", lst)] #M = moisture; 
Clst <- lst[grep("_c", lst)] #C = temperature;

#read csv table with ranked historical years
weather_rank <- read.table('./layers/weather/WeatherRanked.csv', header = T, stringsAsFactors = F, sep=',')

##future weather scenarios: if current simulation year is past the available data, use future climate 
##(1) reading GCM projections from another raster stack  ##TODO!!!
##(2) using available historical COMPLETE years (1990-2007), ranked by suitability

#CASE 1: no future weather needed
if(is.na(opt$scenario) & end <= last_yr){
  Mlst <- grep(paste(as.character(seq(start,end)), collapse="|"), Mlst, value=TRUE)  #use only the raster files matching the years of interest
  Clst <- grep(paste(as.character(seq(start,end)), collapse="|"), Clst, value=TRUE) #use only the raster files matching the years of interest  
  Mstack <- stack(Mlst)  
  Cstack <- stack(Clst) 
#CASE 2: RANDOM future weather scenario  
}else if(opt$scenario == 'random'){
  if(end <= last_yr) stop('you specified a future weather scenario BUT the end year is not a future year!')
  yrs <- sample(weather_rank[,1], size = end - last_yr, replace = T)
  Mlst_1 <- grep(paste(as.character(seq(start,end)), collapse="|"), Mlst, value=TRUE)
  Mlst_2 <- unlist(lapply(yrs, FUN=function(x){grep(as.character(x), Mlst, value=TRUE)}))
  Mlst <- c(Mlst_1, Mlst_2)
  Mstack <- stack(Mlst) 
  Clst_1 <- grep(paste(as.character(seq(start,end)), collapse="|"), Clst, value=TRUE)
  Clst_2 <- unlist(lapply(yrs, FUN=function(x){grep(as.character(x), Clst, value=TRUE)}))
  Clst <- c(Clst_1, Clst_2)
  Cstack <- stack(Clst) 
#CASE 3: UPPER 50% (favorable) future weather scenario
}else if(opt$scenario == 'favorable'){
  yrs <- sample(weather_rank[1:9,1], size = end - last_yr, replace = T)
  if(end <= last_yr) stop('you specified a future weather scenario BUT the end year is not a future year!')
  yrs <- sample(weather_rank[,1], size = end - last_yr, replace = T)
  Mlst_1 <- grep(paste(as.character(seq(start,end)), collapse="|"), Mlst, value=TRUE)
  Mlst_2 <- unlist(lapply(yrs, FUN=function(x){grep(as.character(x), Mlst, value=TRUE)}))
  Mlst <- c(Mlst_1, Mlst_2)
  Mstack <- stack(Mlst) 
  Clst_1 <- grep(paste(as.character(seq(start,end)), collapse="|"), Clst, value=TRUE)
  Clst_2 <- unlist(lapply(yrs, FUN=function(x){grep(as.character(x), Clst, value=TRUE)}))
  Clst <- c(Clst_1, Clst_2)
  Cstack <- stack(Clst) 
#CASE 4: LOWER 50% (unfavorable) future weather scenario
}else{
  yrs <- sample(weather_rank[10:18, 1], size = end - last_yr, replace = T)
  if(end <= last_yr) stop('you specified a future weather scenario BUT the end year is not a future year!')
  yrs <- sample(weather_rank[,1], size = end - last_yr, replace = T)
  Mlst_1 <- grep(paste(as.character(seq(start,end)), collapse="|"), Mlst, value=TRUE)
  Mlst_2 <- unlist(lapply(yrs, FUN=function(x){grep(as.character(x), Mlst, value=TRUE)}))
  Mlst <- c(Mlst_1, Mlst_2)
  Mstack <- stack(Mlst) 
  Clst_1 <- grep(paste(as.character(seq(start,end)), collapse="|"), Clst, value=TRUE)
  Clst_2 <- unlist(lapply(yrs, FUN=function(x){grep(as.character(x), Clst, value=TRUE)}))
  Clst <- c(Clst_1, Clst_2)
  Cstack <- stack(Clst) 
}

##Seasonality: Do you want the spread to be limited to certain months?
ss <- opt$seasonal   #'YES' or 'NO'
if (ss == 'YES') months_msk <- paste('0', 1:9, sep='') #1=January 9=September

##Wind: Do you want the spread to be affected by wind?
if (!(opt$wind %in% c('YES', 'NO'))) stop('You must specify whether you want spread by wind or not: use either YES or NO')

#open window screen
windows(width = 10, height = 10, xpos = 350, ypos = 50, buffered = FALSE)
#quartz()  #use this on Mac OSX
#x11()     #use this on Linux (not tested!)

#plot background image
plot(bkr_img, xaxs = "i", yaxs = "i")

#plot coordinates for plotting text:
xpos <- (bbox(I_umca_rast)[1,2] + bbox(I_umca_rast)[1,1]) / 2
ypos <- bbox(I_umca_rast)[2,2] - 150

#time counter to access pos index in weather raster stacks
cnt <- 0 

## ----> MAIN SIMULATION LOOP (weekly time steps) <------
for (tt in tstep){
  
  #split date string for raster time stamp
  split_date = unlist(strsplit(tt, '-'))
  
  if (tt == tstep[1]) {
    
    if(!any(susceptible_quag > 0) & !any(susceptible_quke > 0)) stop('Simulation ended. There are no more susceptible oaks on the landscape!')
    
    ##CALCULATE OUTPUT TO PLOT: 
    # 1) values as % infected
    I_oaks_rast[] <- ifelse(I_oaks_rast[] == 0, NA, round(I_oaks_rast[]/oaks_rast[], 1))
    
    # 2) values as number of infected per cell
    #I_oaks_rast[] <- ifelse(I_oaks_rast[] == 0, NA, I_oaks_rast[])
    
    # 3) values as 0 (non infected) and 1 (infected) cell
    #I_oaks_rast[] <- ifelse(I_oaks_rast[] > 0, 1, 0) 
    #I_oaks_rast[] <- ifelse(I_oaks_rast[] > 0, 1, NA) 
    
    #PLOT: overlay current plot on background image
    bks <- c(0, 0.25, 0.5, 0.75, 1)
    #colors <- c("yellow","gold","orange","red")
    image(I_oaks_rast, breaks=bks, col=rev(heat.colors(length(bks)-1, alpha=1)), add=T, axes=F, box=F, ann=F, legend=F, useRaster=T)
    boxed.labels(xpos, ypos, tt, bg="white", border=NA, font=2)
    
    #WRITE TO FILE:
    I_oaks_rast_sp <- as(I_oaks_rast, 'SpatialGridDataFrame')
    writeRAST(I_oaks_rast_sp, vname=paste(opt$output, '_', sprintf(formatting_str, 0), sep=''), overwrite=TRUE) #write to GRASS raster file
	  execGRASS('r.timestamp', map=paste(opt$output, '_', sprintf(formatting_str, cnt), sep=''), date=paste(split_date[3], months_names[as.numeric(split_date[2])], split_date[1]))
    
    #writeRaster(I_oaks_rast, filename=paste('./', fOutput, '/', opt$output, '_', sprintf(formatting_str, 0), sep=''), format='HFA', datatype='FLT4S', overwrite=TRUE) # % infected as output
    #writeRaster(I_oaks_rast, filename=paste('./', fOutput, '/', opt$output, '_', sprintf(formatting_str, 0), sep=''), format='HFA', datatype='INT1U', overwrite=TRUE) # nbr. infected hosts as output
    #writeRaster(I_oaks_rast, filename=paste('./', fOutput, '/', opt$output, '_', sprintf(formatting_str, 0), sep=''), format='HFA', datatype='LOG1S', overwrite=TRUE)  # 0=non infected 1=infected output
    
  }else{
    
    
    #check if there are any susceptible trees left on the landscape (IF NOT continue LOOP till the end)
    if(!any(susceptible_quag > 0) & !any(susceptible_quke > 0)) break
    
    #update week counter
    cnt <- cnt + 1
    
    #is current week time step within a spread month (as defined by input parameters)?
    if (ss == 'YES' & !any(substr(tt,6,7) %in% months_msk)) next
          
    #Total weather suitability:
    W <- as.matrix(Mstack[[cnt]] * Cstack[[cnt]])
    
    #GENERATE SPORES:  
    #integer matrix
    spores_mat <- SporeGenCpp(infected_umca, W, rate = opt$spore_rate) #rate: spores/week for each infected host (4.4 default)
    
    #SPORE DISPERSAL:  
    #'List'
    if (opt$wind == 'YES') {
      
      #Check if predominant wind direction has been specified correctly:
      if (!(opt$pwdir %in% c('N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'))) stop('A predominant wind direction must be specified: N, NE, E, SE, S, SW, W, NW')
      out <- SporeDispCppWind_mh(spores_mat, S_UM=susceptible_umca, S_QA=susceptible_quag, S_QK=susceptible_quke, 
                                I_UM=infected_umca, I_QA=infected_quag, I_QK=infected_quke, LVT=livetree_matr, W, rs=res_win, rtype='Cauchy', scale1=20.57, wdir=opt$pwdir, kappa=2)
    
    }else{
      out <- SporeDispCpp_mh(spores_mat, S_UM=susceptible_umca, S_QA=susceptible_quag, S_QK=susceptible_quke, 
                             I_UM=infected_umca, I_QA=infected_quag, I_QK=infected_quke, LVT=livetree_matr, W, rs=res_win, rtype='Cauchy', scale1=20.57) ##TO DO
    }  
    
    #update R matrices:
    #UMCA
    susceptible_umca <- out$S_UM 
    infected_umca <- out$I_UM 
    #QUAG
    susceptible_quag <- out$S_QA 
    infected_quag <- out$I_QA
    #QUKE
    susceptible_quke <- out$S_QK 
    infected_quke <- out$I_QK 
    
    ##CALCULATE OUTPUT TO PLOT:
    I_oaks_rast[] <- infected_quag + infected_quke
    
    # 1) values as % infected
    I_oaks_rast[] <- ifelse(I_oaks_rast[] == 0, NA, round(I_oaks_rast[]/oaks_rast[], 1))
    
    # 2) values as number of infected per cell
    #I_oaks_rast[] <- ifelse(I_oaks_rast[] == 0, NA, I_oaks_rast[])
    
    # 3) values as 0 (non infected) and 1 (infected) cell
    #I_oaks_rast[] <- ifelse(I_oaks_rast[] > 0, 1, 0) 
    #I_oaks_rast[] <- ifelse(I_oaks_rast[] > 0, 1, NA) 
        
    if (cnt %% opt$nth_output == 0){
      
      #PLOT: overlay current plot on background image
      bks <- c(0, 0.25, 0.5, 0.75, 1)
      #colors <- c("yellow","gold","orange","red")
      image(I_oaks_rast, breaks=bks, col=rev(heat.colors(length(bks)-1, alpha=.5)), add=T, axes=F, box=F, ann=F, legend=F, useRaster=T)
      boxed.labels(xpos, ypos, tt, bg="white", border=NA, font=2)
      
      #WRITE TO FILE:
      I_oaks_rast_sp <- as(I_oaks_rast, 'SpatialGridDataFrame')
      writeRAST(I_oaks_rast_sp, vname=paste(opt$output, '_', sprintf(formatting_str, 0), sep=''), overwrite=TRUE) #write to GRASS raster file
      execGRASS('r.timestamp', map=paste(opt$output, '_', sprintf(formatting_str, cnt), sep=''), date=paste(split_date[3], months_names[as.numeric(split_date[2])], split_date[1]))
      
      #writeRaster(I_oaks_rast, filename=paste('./', fOutput, '/', opt$output, '_', sprintf(formatting_str, 0), sep=''), format='HFA', datatype='FLT4S', overwrite=TRUE) # % infected as output
      #writeRaster(I_oaks_rast, filename=paste('./', fOutput, '/', opt$output, '_', sprintf(formatting_str, 0), sep=''), format='HFA', datatype='INT1U', overwrite=TRUE) # nbr. infected hosts as output
      #writeRaster(I_oaks_rast, filename=paste('./', fOutput, '/', opt$output, '_', sprintf(formatting_str, 0), sep=''), format='HFA', datatype='LOG1S', overwrite=TRUE)  # 0=non infected 1=infected output
      
    }
    
  }
  
}

message("Spread model finished")





