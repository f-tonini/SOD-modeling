#---------------------------------------------------------------------------------------------------------------
# Name:         SOD_aniso_clim_RGRASS.r
# Purpose:      Lattice-based simulation of the climate-driven anisotropic spread of pathogen P. ramorum over a heterogeneous landscape.
# Author:       Francesco Tonini
# Email:        ftonini84@gmail.com
# Created:      01/07/2015
# Copyright:    (c) 2015 by Francesco Tonini
# License:      GNU General Public License (GPL)
# Software:     Tested successfully using R version 3.0.2 (http://www.r-project.org/)
#-----------------------------------------------------------------------------------------------------------------------

#install packages
#install.packages(c("rgdal","raster","lubridate","CircStats","Rcpp", "rgrass7", "optparse", "plotrix", "ncdf4"))

#load packages:
suppressPackageStartupMessages(library(raster))    #Raster operation and I/O. Depends R (≥ 2.15.0)
suppressPackageStartupMessages(library(rgdal))     #Geospatial data abstraction library. Depends R (≥ 2.14.0)
suppressPackageStartupMessages(library(lubridate)) #Make dealing with dates a little easier. Depends R (≥ 3.0.0)
suppressPackageStartupMessages(library(CircStats)) #Circular Statistics - Von Mises distribution
suppressPackageStartupMessages(library(Rcpp))      #Seamless R and C++ Integration. Depends R (≥ 3.0.0)
suppressPackageStartupMessages(library(rgrass7))   #Interface Between GRASS 7 GIS and R. Depends R (≥ 2.12)
suppressPackageStartupMessages(library(optparse))  #Parse args from command line
suppressPackageStartupMessages(library(plotrix))   #Add text annotations to plot
suppressPackageStartupMessages(library(ncdf4))   #work with NetCDF datasets

##Define the main working directory based on the current script path
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
base_name <- dirname(sub(file_arg, "", initial_options[grep(file_arg, initial_options)]))
setwd(paste(sep="/", base_name, ".."))

##Use an external source file w/ all modules (functions) used within this script. 
##Use FULL PATH if source file is not in the same folder w/ this script
source('./scripts/myfunctions_SOD.r')
sourceCpp("./scripts/myCppFunctions.cpp") #for C++ custom functions

###Input simulation parameters: #####
option_list = list(
  make_option(c("-u","--umca"), action="store", default=NA, type='character', help="input bay laurel (UMCA) raster map"),
  make_option(c("-ld","--lide"), action="store", default=NA, type='character', help="input Tanoak (LIDE) raster map"),
  make_option(c("-ac","--acma"), action="store", default=NA, type='character', help="input Bigleaf Maple (ACMA) raster map"),
  make_option(c("-ar","--arme"), action="store", default=NA, type='character', help="input Madrone (ARME) raster map"),
  make_option(c("-ae","--aeca"), action="store", default=NA, type='character', help="input California Buckeye (AECA) raster map"),
  make_option(c("-psm","--psme"), action="store", default=NA, type='character', help="input Douglas fir (PSME) raster map"),
  make_option(c("-se","--sese"), action="store", default=NA, type='character', help="input Redwood (SESE) raster map"),
  make_option(c("-ok","--oaks"), action="store", default=NA, type='character', help="input SOD-oaks raster map"),  
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
  make_option(c("-scn","--scenario"), action="store", default=NA, type='character', help="future weather scenario: random, favorable, unfavorable")
)

opt = parse_args(OptionParser(option_list=option_list))

##################################################################
##INPUT RASTERS: abundance (tree density per hectare)

#--- All live trees ---
lvtree_rast <- readRAST(opt$livetree)
lvtree_rast <- raster(lvtree_rast) #transform 'sp' obj to 'raster' obj

#--- P. ramorum ONLY affected species: ---
#UMCA
umca_rast <- readRAST(opt$umca)
umca_rast <- raster(umca_rast)  
#SESE3 (Redwood)
sese_rast <- readRAST(opt$sese)
sese_rast <- raster(sese_rast)
#ACMA3 (Bigleaf Maple)
acma_rast <- readRAST(opt$acma)
acma_rast <- raster(acma_rast)
#AECA (California Buckeye)
aeca_rast <- readRAST(opt$aeca)
aeca_rast <- raster(aeca_rast)
#ARME (Madrone)
arme_rast <- readRAST(opt$arme)
arme_rast <- raster(arme_rast)
#PSME (Douglas fir)
psme_rast <- readRAST(opt$psme)
psme_rast <- raster(psme_rast)

#--- P. ramorum host and SOD-affected species: ---
#LIDE3 (Tanoak)
lide_rast <- readRAST(opt$lide)
lide_rast <- raster(lide_rast)

#--- SOD-affected species: ---
#oaks (QUAG, QUCH2, QUKE)
oaks_rast <- readRAST(opt$oaks)
oaks_rast <- raster(oaks_rast)

##########################################################
#### SPREAD SCORE ####
#Read table with P.ramorum host spread potential scores:
Hscore <- read.table("./HostScore.txt", header = T, stringsAsFactors = F, sep=',')
#subset only species with score > 0 and that are present in the study area
Hscore <- Hscore[Hscore$Host_score > 0 & !is.na(Hscore$LEMMA_ABBR), ]
#calculate and add relative score weights:
Hscore$weight <- round(Hscore$Host_score / Hscore$Host_score[Hscore$LEMMA_ABBR == "UMCA"], 2)
#create a list with species and associated weight:
wLst <- split(Hscore$weight, Hscore$LEMMA_ABBR)
######################################################

#max value for plotting
mx <- cellStats(oaks_rast, stat='max') 
#raster resolution (all rasters are the same resolution so take any as reference)
res_win <- res(umca_rast)[1]

###################################
### INFECTED AND SUSCEPTIBLES ####

##Initial infection (OAKS):
I_oaks_rast <- readRAST(opt$sources)  #measured
I_oaks_rast <- raster(I_oaks_rast) 

#define matrices for infected and susceptible species of interest
I_oaks <- as.matrix(I_oaks_rast)
S_oaks <- as.matrix(oaks_rast - I_oaks_rast)
I_umca <- matrix(0, nrow=res_win, ncol=res_win)
S_umca <- as.matrix(umca_rast)
I_lide <- matrix(0, nrow=res_win, ncol=res_win)
S_lide <- as.matrix(lide_rast)
I_acma <- matrix(0, nrow=res_win, ncol=res_win)
S_acma <- as.matrix(acma_rast)
I_arme <- matrix(0, nrow=res_win, ncol=res_win)
S_arme <- as.matrix(arme_rast)
I_aeca <- matrix(0, nrow=res_win, ncol=res_win)
S_aeca <- as.matrix(aeca_rast)
I_psme <- matrix(0, nrow=res_win, ncol=res_win)
S_psme <- as.matrix(psme_rast)
I_sese <- matrix(0, nrow=res_win, ncol=res_win)
S_sese <- as.matrix(sese_rast)

##Initialize infected trees for each species (!!NEEDED UNLESS EMPIRICAL INFO IS AVAILABLE!!)
if(any(S_umca[I_oaks > 0] > 0)) I_umca[I_oaks > 0] <- mapply(function(x,y) ifelse(x > y, min(c(x,y*2)), x), 
                                                             S_umca[I_oaks > 0], I_oaks[I_oaks > 0]) 
if(any(S_lide[I_oaks > 0] > 0)) I_lide[I_oaks > 0] <- mapply(function(x,y) ifelse(x > y, min(c(x,y*2)), x), 
                                                             S_lide[I_oaks > 0], I_oaks[I_oaks > 0])
if(any(S_acma[I_oaks > 0] > 0)) I_acma[I_oaks > 0] <- mapply(function(x,y) ifelse(x > y, min(c(x,y)), x), 
                                                             S_acma[I_oaks > 0], I_oaks[I_oaks > 0])
if(any(S_arme[I_oaks > 0] > 0)) I_arme[I_oaks > 0] <- mapply(function(x,y) ifelse(x > y, min(c(x,y)), x), 
                                                             S_arme[I_oaks > 0], I_oaks[I_oaks > 0])
if(any(S_aeca[I_oaks > 0] > 0)) I_aeca[I_oaks > 0] <- mapply(function(x,y) ifelse(x > y, min(c(x,y)), x), 
                                                             S_aeca[I_oaks > 0], I_oaks[I_oaks > 0])
if(any(S_psme[I_oaks > 0] > 0)) I_psme[I_oaks > 0] <- mapply(function(x,y) ifelse(x > y, min(c(x,y)), x), 
                                                             S_psme[I_oaks > 0], I_oaks[I_oaks > 0])
if(any(S_sese[I_oaks > 0] > 0)) I_sese[I_oaks > 0] <- mapply(function(x,y) ifelse(x > y, min(c(x,y)), x), 
                                                             S_sese[I_oaks > 0], I_oaks[I_oaks > 0])

##update susceptible matrices by subtracting the initialized infections 
S_umca <- S_umca - I_umca 
S_lide <- S_lide - I_lide
S_acma <- S_acma - I_acma
S_arme <- S_arme - I_arme
S_aeca <- S_aeca - I_aeca
S_psme <- S_psme - I_psme
S_sese <- S_sese - I_sese

##define matrix for immune live trees
N_live <- as.matrix(lvtree_rast)

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
#weather coefficients
mcf.array <- ncvar_get(nc_open('./layers/weather/weatherCoeff_2000_2014.nc'),  varid = "Mcoef") #M = moisture;
ccf.array <- ncvar_get(nc_open('./layers/weather/weatherCoeff_2000_2014.nc'),  varid = "Ccoef") #C = temperature;

##Seasonality: Do you want the spread to be limited to certain months?
ss <- 'YES'   #'YES' or 'NO'
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
    I_oaks_rast_sp <- as(I_oaks_rast, 'SpatialGridDataFrame')
    writeRAST(I_oaks_rast_sp, vname=paste(opt$output, '_', sprintf(formatting_str, cnt), sep=''), overwrite=TRUE) #write to GRASS raster file
	  execGRASS('r.timestamp', map=paste(opt$output, '_', sprintf(formatting_str, cnt), sep=''), date=paste(split_date[3], months_names[as.numeric(split_date[2])], split_date[1]))

    #writeRaster(I_oaks_rast, filename=paste('./', fOutput, '/', opt$output, '_', sprintf(formatting_str, cnt), sep=''), format='HFA', datatype='FLT4S', overwrite=TRUE) # % infected as output
    #writeRaster(I_oaks_rast, filename=paste('./', fOutput, '/', opt$output, '_', sprintf(formatting_str, cnt), sep=''), format='HFA', datatype='INT1U', overwrite=TRUE) # nbr. infected hosts as output
    #writeRaster(I_oaks_rast, filename=paste('./', fOutput, '/', opt$output, '_', sprintf(formatting_str, cnt), sep=''), format='HFA', datatype='LOG1S', overwrite=TRUE)  # 0=non infected 1=infected output
    
  }else{
    
    #check if there are any susceptible oaks left on the landscape (IF NOT continue LOOP till the end)
    if(!any(S_oaks > 0)) break
    
    #update week counter
    cnt <- cnt + 1
    
    #is current week time step within a spread month (as defined by input parameters)?
    if (ss == 'YES' & !any(substr(tt,6,7) %in% months_msk)) {
      if (cnt == length(tstep) - 1) {
        #WRITE TO FILE:
        I_oaks_rast_sp <- as(I_oaks_rast, 'SpatialGridDataFrame')
        writeRAST(I_oaks_rast_sp, vname=paste(opt$output, '_', sprintf(formatting_str, cnt), sep=''), overwrite=TRUE) #write to GRASS raster file
        execGRASS('r.timestamp', map=paste(opt$output, '_', sprintf(formatting_str, cnt), sep=''), date=paste(split_date[3], months_names[as.numeric(split_date[2])], split_date[1]))
        I_umca_rast <- I_oaks_rast
        I_umca_rast[] <- I_umca
        I_umca_rast[] <- ifelse(I_umca_rast[] == 0, NA, I_umca_rast[])
        I_umca_rast_sp <- as(I_umca_rast, 'SpatialGridDataFrame')
        writeRAST(I_umca_rast_sp, vname=paste(opt$output, 'umca_', sprintf(formatting_str, cnt), sep=''), overwrite=TRUE) #write to GRASS raster file
        execGRASS('r.timestamp', map=paste(opt$output, 'umca_', sprintf(formatting_str, cnt), sep=''), date=paste(split_date[3], months_names[as.numeric(split_date[2])], split_date[1]))        
      }
      next
    }
          
    #Total weather suitability:
    W <- mcf.array[,,cnt] * ccf.array[,,cnt]
    
    #GENERATE SPORES:  
    #integer matrix
    set.seed(42)
    spores_mat <- SporeGenCpp_MH(I_umca, I_lide, I_acma, I_arme, I_aeca, 
                                 I_psme, I_sese, wLst, W, rate = opt$spore_rate) #rate: spores/week for each infected host (4.4 default)
    
    #SPORE DISPERSAL:  
    #'List'
    if (opt$wind == 'YES') {
      
      #Check if predominant wind direction has been specified correctly:
      if (!(opt$pwdir %in% c('N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'))) stop('A predominant wind direction must be specified: N, NE, E, SE, S, SW, W, NW')
      out <- SporeDispCppWind_MH(spores_mat, S_UM=S_umca, S_LD=S_lide, S_AC=S_acma, S_AR=S_arme, S_AE=S_aeca, 
                                 S_PS=S_psme, S_SE=S_sese, S_OK=S_oaks, I_UM=I_umca, I_LD=I_lide, I_AC=I_acma, 
                                 I_AR=I_arme, I_AE=I_aeca, I_PS=I_psme, I_SE=I_sese, I_OK=I_oaks, N_LVE=N_live, 
                                 W, rs=res_win, rtype='Cauchy', scale1=20.57, wdir=opt$pwdir, kappa=2)
      
    }else{
      out <- SporeDispCpp_MH(spores_mat, S_UM=S_umca, S_LD=S_lide, S_AC=S_acma, S_AR=S_arme, S_AE=S_aeca, 
                             S_PS=S_psme, S_SE=S_sese, S_OK=S_oaks, I_UM=I_umca, I_LD=I_lide, I_AC=I_acma, 
                             I_AR=I_arme, I_AE=I_aeca, I_PS=I_psme, I_SE=I_sese, I_OK=I_oaks, N_LVE=N_live,
                             W, rs=res_win, rtype='Cauchy', scale1=20.57) ##TO DO
    }  
    
    #update R matrices:
    #UMCA
    S_umca <- out$S_UM 
    I_umca <- out$I_UM 
    #LIDE
    S_lide <- out$S_LD 
    I_lide <- out$I_LD     
    #ACMA
    S_acma <- out$S_AC 
    I_acma <- out$I_AC     
    #ARME
    S_arme <- out$S_AR 
    I_arme <- out$I_AR     
    #AECA
    S_aeca <- out$S_AE 
    I_aeca <- out$I_AE     
    #PSME
    S_psme <- out$S_PS 
    I_psme <- out$I_PS     
    #SESE
    S_sese <- out$S_SE 
    I_sese <- out$I_SE    
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
        
    if (cnt %% opt$nth_output == 0){
      
      #PLOT: overlay current plot on background image
      #bks <- c(0, 0.25, 0.5, 0.75, 1)
      #my_palette <- colorRampPalette(c("springgreen", "yellow1", "orange", "red1"))(n = 4)
      #image(I_oaks_rast, breaks=bks, col=addalpha(my_palette, .5), add=T, axes=F, box=F, ann=F, legend=F, useRaster=T)
      bks <- seq(0, mx, length = 10)
      image(I_oaks_rast, breaks=bks, col=rev(heat.colors(length(bks)-1, alpha=1)), add=T, axes=F, box=F, ann=F, legend=F, useRaster=T)
      boxed.labels(xpos, ypos, tt, bg="white", border=NA, font=2)
      
      #WRITE TO FILE:
      I_oaks_rast_sp <- as(I_oaks_rast, 'SpatialGridDataFrame')
      writeRAST(I_oaks_rast_sp, vname=paste(opt$output, '_', sprintf(formatting_str, cnt), sep=''), overwrite=TRUE) #write to GRASS raster file
      execGRASS('r.timestamp', map=paste(opt$output, '_', sprintf(formatting_str, cnt), sep=''), date=paste(split_date[3], months_names[as.numeric(split_date[2])], split_date[1]))
      
      #writeRaster(I_oaks_rast, filename=paste('./', fOutput, '/', opt$output, '_', sprintf(formatting_str, cnt), sep=''), format='HFA', datatype='FLT4S', overwrite=TRUE) # % infected as output
      #writeRaster(I_oaks_rast, filename=paste('./', fOutput, '/', opt$output, '_', sprintf(formatting_str, cnt), sep=''), format='HFA', datatype='INT1U', overwrite=TRUE) # nbr. infected hosts as output
      #writeRaster(I_oaks_rast, filename=paste('./', fOutput, '/', opt$output, '_', sprintf(formatting_str, cnt), sep=''), format='HFA', datatype='LOG1S', overwrite=TRUE)  # 0=non infected 1=infected output
      
    }
    
  }
  
}

message("Spread model finished")





