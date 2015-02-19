#--------------------------------------------------------------------------------
# Name:         SOD_original.r
# Purpose:      Lattice-based simulation of the spread of pathogen P. ramorum over a heterogeneous landscape, using focal ('moving window') approach.
# Author:       Francesco Tonini
# Email:        ftonini84@gmail.com
# Created:      01/07/2015
# Copyright:    (c) 2015 by Francesco Tonini
# License:      GNU General Public License (GPL)
# Software:     Tested successfully using R version 3.0.2 (http://www.r-project.org/)
#-----------------------------------------------------------------------------------------------------------------------

#set.seed(2000)

##Define the main working directory
##Make sure to specify the appropriate path using either / or \\ to specify the path 
#setwd("path_to_desired_folder")
setwd("D:/SOD Simulation")

#Path to folders in which you want to save all your vector & raster files
fOutput <- 'output'

##Create a physical copy of the subdirectory folder(s) where you save your output
##If the directory already exists it gives a warning, BUT we can suppress it using showWarnings = FALSE
dir.create(fOutput, showWarnings = FALSE)

##Use an external source file w/ all modules (functions) used within this script. 
##Use FULL PATH if source file is not in the same folder w/ this script
source('./scripts/myfunctions_SOD.R')

##Load all required libraries
print('Loading required libraries...')
load.packages()

###Let's set all simulation parameters:
start_time <- 1
end_time <- 12

#ELEVATION HAS TO BE IN METERS!! IF NOT, PLEASE CONVERT UNITS USING LINE BELOW
#DEM <- raster('./layers/DEM_Sonoma90.img')
DEM <- raster('./layers/DEM_test.img')
#DEM <- DEM*0.3048  #feet to meters

WindStrgth <- HI.raster
WindStrgth[] <- 4 #constant wind strength always

##read study area raster with count of S (=susceptible trees)
#HI.raster <- raster('./layers/HI_Sonoma.img')
HI.raster <- raster('./layers/test.img')

#max number of trees that can be infected in each cell
Nmax <- HI.raster[]

#initialize list with counts of I (=infected trees)
I.lst <- rep(0, length(Nmax))

#initialize raster object for I (=infected trees) to plot at each time step
I.raster <- HI.raster  

#initialize list with counts of pathogen spores
spores.lst <- rep(0, length(Nmax))

#I.new.raster <- HI.raster   

#USE ONLY AS PART OF STRATEGY (1)************see within FOR LOOP***********************
#extract xy coords of cell centroids from study area raster
#xyAOI <- xyFromCell(HI.raster, 1:ncell(HI.raster))
#x <- xyAOI[,1]
#y <- xyAOI[,2]
#create matrix with all combinations of cell index (row-column) in the study area
#xyAOIPos <- expand.grid(1:ncol(HI.raster),1:nrow(HI.raster))
#create 'owin' object defining the spatial window of the study area (required by 'spatstat' pkg)
#AOI_owin <- owin(c(extent(HI.raster)[1],extent(HI.raster)[2]), c(extent(HI.raster)[3],extent(HI.raster)[4]), mask=matrix(TRUE, nrow(HI.raster),ncol(HI.raster))) 
#*****************************************************************************************************************************************************************


#define how many sources of infections you want to initialize on the landscape
inf.sources <- 2 #number of pixels

#randomly sample the index of cells to be source of infections
#inf.index <- sample(which(HI.raster[] > 0), size = inf.sources)
inf.index <- sample(which(Nmax > 0), size = inf.sources)

#randomize the initial I (=infected trees) counts (this does NOT have to exceed Nmax)
I.lst[inf.index] <- sapply(Nmax[inf.index], FUN=function(x)sample(1:x, size=1)) 
#I.raster[inf.index] <- sapply(Nmax[inf.index], FUN=function(x) sample(1:x, size=1)) 

#update list w/ counts of S (=susceptible trees)
S.lst <- Nmax - I.lst

### DISPERSAL KERNEL:

#parameters for Cauchy mixture estimated via MCMC in Meentemeyer et al. 2011 (20.57 meters, 9.5 km, gamma=0.9947)
#probKernel <- Kernel.gen(x = HI.raster, d, type = "Cauchy Mixture", scale=c(20.57,9500), gamma=0.9947)   #use "exp" for exponential kernel and "gauss" for gaussian kernel  
#probKernel <- Kernel.gen(x = HI.raster, d, type = "Cauchy", scale=20.57)   #use "Exponential" for exponential kernel and "Gauss" for gaussian kernel  

d <- 500
probKernel <- kernel2D(HI.raster, d) #d = meters or units of analysis
dist <- kernel2DToDistance(probKernel, HI.raster) ##USE d for CIRCULAR KERNELS (not implemented yet)! 

DEM.buff <- buffer.areas(DEM[], HI.raster, width.row = nrow(probKernel), width.col = ncol(probKernel))
susceptible <- buffer.areas(S.lst, HI.raster, width.row = nrow(probKernel), width.col = ncol(probKernel))
infected <- buffer.areas(I.lst, HI.raster, width.row = nrow(probKernel), width.col = ncol(probKernel))
spores <- buffer.areas(spores.lst, HI.raster, width.row = nrow(probKernel), width.col = ncol(probKernel))

##LOOP for each month (or whatever chosen time unit)
for (tt in start_time:end_time){
  
  print(paste('Time', tt))
  
  if (tt == start_time) {
    
    I.raster[] <- round(I.lst/Nmax,1)
    #writeRaster(I.raster, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='INT1U', overwrite=TRUE) #nbr. infected hosts as output
    #writeRaster(I.raster > 0, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='LOG1S', overwrite=TRUE)  #0=non infected 1=infected output
    plot(I.raster, main=paste("Time ", tt, sep=''))
    
  }else{
    
    #IF THERE AREN'T ANY SUSCEPTIBLES LEFT ON THE LANDSCAPE, SKIP TO NEXT TIME STEP...
    #************************************************************************************
    
    #if(!any(S.lst >0)){
    #  plot(I.raster, main=paste("Time ", tt, sep=''))
    #  next 
    #} 
    
    if(!any(susceptible >0, na.rm=T)){
      plot(I.raster, main=paste("Time ", tt, sep=''))
      next 
    }  
    #********************************************************************************
    
    
    #Within each infected cell (I > 0) draw random number of infections ~Poisson(lambda=rate of spore production) for each infected host. 
    #Take SUM for total infections produced by each cell.  
    spores.lst[I.lst > 0] <- sapply(I.lst[I.lst > 0], FUN=new.infections.gen, rate=4.4*4)  #4.4 * 4 convert approximately to spores/month    
    
    
    #STRATEGY (1)******************************************************************** 
    #GENERATE N DISP DISTANCES AND ANGLE FOR EACH SPORE UNIT (slower!!)
    #new.spores.lst <- radialDisp(xyAOIPos, AOI_owin)
    #update list of I (=infected trees) with dispersed spores
    #I.lst[S.lst > 0] <- I.lst[S.lst > 0] + new.spores.lst[S.lst > 0] 
    #********************************************************************************
    
    #STRATEGY (2)******************************************************************** 
    #USE GETFOCALVALUES() INSTEAD OF MOVING WINDOW (slower!!)
    #I.new.raster[] <- spores.lst
    #I.new.raster[] <- myFocal.fun(I.new.raster, probKernel, ngb=dim(probKernel))
    #only the cells containing SUSCEPTIBLES can be challenged by these new infections
    #I.lst[S.lst > 0] <- I.lst[S.lst > 0] + I.new.raster[S.lst > 0] 
    #********************************************************************************
    
    #STRATEGY (3)********************************************************************************** 
    #USE FOCAL OPERATION AND MOVING WINDOW w/ KERNEL WEIGHTS 
    #focal moving window:
    #I.new.raster[] <- spores.lst
    #I.new.raster <- round(focal(I.new.raster, w=probKernel, fun=sum, pad=TRUE, padValue=0))    
    
    #out <- mapply(FUN=randomU, I.new.raster[], I.lst, S.lst)
    #I.lst <- unlist(out["I",])
    #S.lst <- unlist(out["S",])
    
    
    #STRATEGY (4)*************************************************************************************
    #2D KERNELS ARE CENTERED ON TOP OF SPORE CELLS
    #THIS IS OPPOSITE TO MOVING WINDOW IN THAT WE DO NOT LOOK AT SPORES DISPERSING INTO A GIVEN CELL
    #BUT RATHER AT SPORES DISPERSING FROM EACH INFECTED CELL. THIS IS ONLY FASTER IF LOOPS ARE CONVERTED TO C++
    #(using with 'Rcpp' package!)
    
    #save spores.lst into right position of spore buffered landscape matrix
    spores[(1 + floor(nrow(probKernel)/2)):(nrow(spores)-floor(nrow(probKernel)/2)),  (1 + floor(ncol(probKernel)/2)):(ncol(spores)-floor(ncol(probKernel)/2))] <- matrix(spores.lst, ncol=ncol(HI.raster), byrow=T)
    
    #c <- 0.8 (c <= 1 fat-tailed kernel)
    out <- dispFunction(spores, susceptible, infected, probKernel, type='Cauchy', c=0.5, scale=20.57) #change this function to C++
    susceptible <- out$S
    infected <- out$I
    
    
    I.raster[] <- infected[(1 + floor(nrow(probKernel)/2)):(nrow(infected)-floor(nrow(probKernel)/2)), (1 + floor(ncol(probKernel)/2)):(ncol(infected)-floor(ncol(probKernel)/2))]
    I.lst <- I.raster[]
    
    I.raster[] <- round(I.raster[]/Nmax,1)
    
    
    #USE THE FOLLOWING ONLY FOR STRATEGY (1) OR (2)!!!!
    #only the cells containing SUSCEPTIBLES can be challenged by these new infections
    #I.lst[S.lst > 0] <- I.lst[S.lst > 0] + I.new.raster[S.lst > 0]     
    #if I (=infected trees) exceeds Nmax then override I with Nmax
    #if(any(I.lst > Nmax)) I.lst[I.lst > Nmax & Nmax > 0] <- Nmax[I.lst > Nmax & Nmax > 0]  
    #update counts of S (=susceptible trees) after new infections
    #S.lst <- Nmax - I.lst
    
    
    #I.raster[] <- round(I.lst/Nmax,1)
    
    #writeRaster(I.raster, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='INT1U', overwrite=TRUE) #nbr. infected hosts as output
    #writeRaster(I.raster > 0, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='LOG1S', overwrite=TRUE)  #0=non infected 1=infected output
    plot(I.raster, main=paste("Time ", tt, sep=''))
    #Sys.sleep(1)
    
  }  
  
  
}  