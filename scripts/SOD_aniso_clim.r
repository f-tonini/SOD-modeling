#--------------------------------------------------------------------------------
# Name:         SOD_aniso.r
# Purpose:      Lattice-based simulation of the spread of pathogen P. ramorum over a heterogeneous landscape.
# Author:       Francesco Tonini
# Email:        ftonini84@gmail.com
# Created:      01/07/2015
# Copyright:    (c) 2015 by Francesco Tonini
# License:      GNU General Public License (GPL)
# Software:     Tested successfully using R version 3.0.2 (http://www.r-project.org/)
#-----------------------------------------------------------------------------------------------------------------------

#install packages
install.packages(c("rgdal","raster","lubridate","CircStats","spatstat","Rcpp"))

#load packages
library(raster)  		#Raster operation and I/O
library(rgdal)	
library(lubridate)
library(CircStats)  #Von Mises distribution
library(spatstat)
library(Rcpp)


set.seed(2000)

##Define the main working directory
##Make sure to specify the appropriate path using either / or \\ to specify the path 
#setwd("path_to_desired_folder")
setwd("D:\\TangibleLandscape\\SOD-modeling")

#Path to folders in which you want to save all your vector & raster files
fOutput <- 'output'

##Create a physical copy of the subdirectory folder(s) where you save your output
##If the directory already exists it gives a warning, BUT we can suppress it using showWarnings = FALSE
dir.create(fOutput, showWarnings = FALSE)

##Use an external source file w/ all modules (functions) used within this script. 
##Use FULL PATH if source file is not in the same folder w/ this script
source('./scripts/myfunctions_SOD.r')
sourceCpp("./scripts/myCppFunctions.cpp") #for C++ custom functions


###Input simulation parameters:
start <- 2004
end <- 2006

dd_start <- as.POSIXlt(as.Date(paste(start,'-01-01',sep='')))
dd_end <- as.POSIXlt(as.Date(paste(end,'-12-31',sep='')))

tstep <- seq(dd_start, dd_end, 'weeks')

#months <- c("January","February","March","April","May","June","July","August","September","October","November","December")
#months_msk <- c("January","February","March","April","May","June","July","August","September")
#tstep <- sapply(months, FUN=function(x) paste(x,start:end,sep=''))
#tstep <- c(t(s))

#MOISTURE
Mstack <- stack()

for (yr in seq(start,end)){
  
  #moisture coeff
  Mlst <- dir(paste('./layers/M/', yr, sep=''), pattern='\\.img$', full.names=T)
  #ordered sequence
  Mlst <- Mlst[match( paste('m',seq(1,52),sep=''),  sub("^([^.]*).*", "\\1", basename(Mlst)) )]
  Mlst_rstck <- stack(Mlst)
  
  Mstack <- stack(Mstack, Mlst_rstck)
  
}

#TEMPERATURE
Cstack <- stack()

for (yr in seq(start,end)){
  
  #moisture coeff
  Clst <- dir(paste('./layers/C/', yr, sep=''), pattern='\\.img$', full.names=T)
  Clst_rstck <- stack(Clst)
  
  Cstack <- stack(Cstack, Clst_rstck)
  
}



##read initial raster (host index) with counts of max available "Susceptible" trees per cell
##counts are integers [0, 100]
Nmax_rast <- raster('./layers/HI_100m.img')
Nmax <- Nmax_rast[]  #integer vector of Smax

#raster resolution
res_win <- res(Nmax_rast)[1]

#empty vector with counts of I (=infected trees)
I_lst <- rep(0, length(Nmax))  #integer

#clone Smax raster to I (=infected trees) raster and spores (=number of spores)
I_rast <- Nmax_rast 

#empty vector with counts of pathogen spores
spores_lst <- rep(0, length(Nmax))  #integer

########################################################
##SOURCES OF INFECTION:
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

#initialize counter to index number of layers in the RasterStack
cnt <- 0

##LOOP for each month (or whatever chosen time unit)
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
    
    breakpoints <- c(0, 0.25, 0.5, 0.75, 1)
    colors <- c("yellow","gold","orange","red")
    plot(I_rast, breaks=breakpoints, col=colors, main=tt)
    
    #WRITE TO FILE:
    #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='FLT4S', overwrite=TRUE) # % infected as output
    #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='INT1U', overwrite=TRUE) # nbr. infected hosts as output
    #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='LOG1S', overwrite=TRUE)  # 0=non infected 1=infected output
    
  }else{
      
    #is current week time step within a spread month? spread month defined as 1-9 (Jan-Sep)
    if (!any(month(tt) %in% seq(1,9))) next
    
    #check if there are any susceptible trees left on the landscape (IF NOT continue LOOP till the end)
    if(!any(susceptible > 0)){
      breakpoints <- c(0, 0.25, 0.5, 0.75, 1)
      colors <- c("yellow","gold","orange","red")
      plot(I_rast, breaks=breakpoints, col=colors, main=tt)
      #WRITE TO FILE:
      #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='FLT4S', overwrite=TRUE) # % infected as output
      #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='INT1U', overwrite=TRUE) # nbr. infected hosts as output
      #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='LOG1S', overwrite=TRUE)  # 0=non infected 1=infected output
      next 
    }
    
    #update counter
    cnt <- cnt + 1 
    
    #Within each infected cell (I > 0) draw random number of infections ~Poisson(lambda=rate of spore production) for each infected host. 
    #Take SUM for total infections produced by each cell.
    
    #compute weather suitability coefficients for current time step
    W <- as.matrix(Mstack[[cnt]] * Cstack[[cnt]])
    
    ##LOOP TO GENERATE SPORES
    #integer vector
    spores_mat <- SporeGen(infected, W, rate = 4.4)
    
    
    #SPORE DISPERSAL:
    out <- SporeDispCpp(spores_mat, S=susceptible, I=infected, rs=res_win, rtype='Cauchy', wtype='VM', wdir='N', scale1=20.57, kappa=2)
    
    
    #(SporeDisp2 seems faster in R alone!!)
    #out <- SporeDisp(spores_mat, S=susceptible, I=infected, rs=res_win, rtype='Cauchy', scale=20.57, wtype='Uniform')  #C++ functions to CONVERT
    #out <- SporeDisp(spores_mat, S=susceptible, I=infected, rs=res_win, rtype='Cauchy', scale=20.57, wtype='VM', wdir='N', kappa=2)  #C++ functions to CONVERT
    #out <- SporeDisp2(spores_mat, S=susceptible, I=infected, rs=res_win, rtype='Cauchy', scale=20.57, wtype='Uniform')  #C++ functions to CONVERT
    #out <- SporeDisp2(spores_mat, S=susceptible, I=infected, rs=res_win, rtype='Cauchy', scale=20.57, wtype='VM', wdir='N', kappa=2)  #C++ functions to CONVERT
    
    susceptible <- out$S 
    infected <- out$I  
    
    I_rast[] <- infected
    I_lst <- I_rast[]
    
    #wipe out spores vector to use in following time steps
    spores_lst <- rep(0, length(Nmax))  #integer
    
    ##CALCULATE OUTPUT TO PLOT:
    # 1) values as % infected
    I_rast[] <- ifelse(I_lst == 0, NA, round(I_lst/Nmax, 1))
    
    # 2) values as number of infected per cell
    #I_rast[] <- ifelse(I_lst == 0, NA, I_lst)    
    
    # 3) values as 0 (non infected) and 1 (infected) cell
    #I_rast[] <- ifelse(I_lst > 0, 1, 0) 
    #I_rast[] <- ifelse(I_lst > 0, 1, NA) 
    
    breakpoints <- c(0, 0.25, 0.5, 0.75, 1)
    colors <- c("yellow","gold","orange","red")
    plot(I_rast, breaks=breakpoints, col=colors, main=tt)
    
    #WRITE TO FILE:
    #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='FLT4S', overwrite=TRUE) # % infected as output
    #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='INT1U', overwrite=TRUE) # nbr. infected hosts as output
    #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='LOG1S', overwrite=TRUE)  # 0=non infected 1=infected output
  }
  
  
}





