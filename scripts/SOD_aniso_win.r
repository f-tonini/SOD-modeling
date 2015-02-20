#--------------------------------------------------------------------------------
# Name:         SOD_aniso_win.r
# Purpose:      Lattice-based simulation of the spread of pathogen P. ramorum over a heterogeneous landscape, using moving window approach.
# Author:       Francesco Tonini
# Email:        ftonini84@gmail.com
# Created:      01/07/2015
# Copyright:    (c) 2015 by Francesco Tonini
# License:      GNU General Public License (GPL)
# Software:     Tested successfully using R version 3.0.2 (http://www.r-project.org/)
#-----------------------------------------------------------------------------------------------------------------------

set.seed(2000)

##Define the main working directory
##Make sure to specify the appropriate path using either / or \\ to specify the path 
#setwd("path_to_desired_folder")
setwd("D:/TangibleLandscape")

#Path to folders in which you want to save all your vector & raster files
fOutput <- 'output'

##Create a physical copy of the subdirectory folder(s) where you save your output
##If the directory already exists it gives a warning, BUT we can suppress it using showWarnings = FALSE
dir.create(fOutput, showWarnings = FALSE)

##Use an external source file w/ all modules (functions) used within this script. 
##Use FULL PATH if source file is not in the same folder w/ this script
source('./scripts/myfunctions_SOD.r')

##Load all required libraries
print('Loading required libraries...')
load.packages()


###Input simulation parameters:
start_time <- 1
end_time <- 12

##read initial raster (host index) with counts of max available "Susceptible" trees per cell
##counts are integers [0, 100]
Nmax_rast <- raster('./layers/HI_100m.img')
Nmax <- Nmax_rast[]  #integer vector of Smax

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

#max cut-off distance for the 'moving window'
d <- 5000 #meters or units of analysis

kernelWin <- kernel2D(Nmax_rast, d) #d = meters or units of analysis
dstMat <- kernel2DToDistance(kernelWin, Nmax_rast)

#buffer matrix by half kernel size on all sides
susceptible <- buffer.areas(S_lst, Nmax_rast, width.row = nrow(kernelWin), width.col = ncol(kernelWin))
infected <- buffer.areas(I_lst, Nmax_rast, width.row = nrow(kernelWin), width.col = ncol(kernelWin))
spores <- buffer.areas(spores_lst, Nmax_rast, width.row = nrow(kernelWin), width.col = ncol(kernelWin))


##LOOP for each month (or whatever chosen time unit)
for (tt in start_time:end_time){
  
  print(paste('Time', tt))
  
  if (tt == start_time) {
    
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
    plot(I_rast, breaks=breakpoints, col=colors, main=paste("Time ", tt, sep=''))
    
    #WRITE TO FILE:
    #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='FLT4S', overwrite=TRUE) # % infected as output
    #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='INT1U', overwrite=TRUE) # nbr. infected hosts as output
    #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='LOG1S', overwrite=TRUE)  # 0=non infected 1=infected output
    
  }else{
    
    #check if there are any susceptible trees left on the landscape (IF NOT continue LOOP till the end)
    if(!any(susceptible > 0)){
      breakpoints <- c(0, 0.25, 0.5, 0.75, 1)
      colors <- c("yellow","gold","orange","red")
      plot(I_rast, breaks=breakpoints, col=colors, main=paste("Time ", tt, sep=''))
      #WRITE TO FILE:
      #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='FLT4S', overwrite=TRUE) # % infected as output
      #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='INT1U', overwrite=TRUE) # nbr. infected hosts as output
      #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='LOG1S', overwrite=TRUE)  # 0=non infected 1=infected output
      next 
    }
    
    #Within each infected cell (I > 0) draw random number of infections ~Poisson(lambda=rate of spore production) for each infected host. 
    #Take SUM for total infections produced by each cell.
    
    #integer vector
    spores_lst[I_lst > 0] <- sapply(I_lst[I_lst > 0], FUN=new.infections.gen, rate=4.4*4)  #4.4 * 4 convert approximately to spores/month  
        
    #save spores_lst into right position of spore buffered landscape matrix
    spores[(1 + floor(nrow(kernelWin)/2)):(nrow(spores)-floor(nrow(kernelWin)/2)),  
           (1 + floor(ncol(kernelWin)/2)):(ncol(spores)-floor(ncol(kernelWin)/2))] <- matrix(spores_lst, ncol=ncol(Nmax_rast), byrow=T)
    
    #disperse spores over landscape
    #c <- 0.8 (c <= 1 fat-tailed kernel)
    spores <- dispFunction(bkgr = spores, dst = dstMat, type='Cauchy', c=0.5, scale=20.57) #change this function to C++
    
    out <- infectFun(bkgr = spores, S = susceptible, I = infected)
    
    susceptible <- out$S
    infected <- out$I
    
    I_rast[] <- infected[(1 + floor(nrow(kernelWin)/2)):(nrow(infected)-floor(nrow(kernelWin)/2)), (1 + floor(ncol(kernelWin)/2)):(ncol(infected)-floor(ncol(kernelWin)/2))]
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
    plot(I_rast, breaks=breakpoints, col=colors, main=paste("Time ", tt, sep=''))
    
    #WRITE TO FILE:
    #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='FLT4S', overwrite=TRUE) # % infected as output
    #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='INT1U', overwrite=TRUE) # nbr. infected hosts as output
    #writeRaster(I_rast, filename=paste('./',fOutput,'/Infected_', tt, '.img',sep=''), format='HFA', datatype='LOG1S', overwrite=TRUE)  # 0=non infected 1=infected output
    
  }  
  
  
}  