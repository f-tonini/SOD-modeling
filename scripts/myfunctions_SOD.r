#--------------------------------------------------------------------------------
# Name:         myfunctions_SOD.r
# Purpose:      Modules (functions) called by the main script
# Author:       Francesco Tonini
# Email:        ftonini84@gmail.com
# Created:      01/07/2015
# Copyright:    (c) 2015 by Francesco Tonini
# License:      GNU General Public License (GPL)
# Software:     Tested successfully using R version 3.0.2 (http://www.r-project.org/)
#-----------------------------------------------------------------------------------------------------------------------

addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}


cauchy.weight <- function(rs, d, scale, gamma = NULL) {
  
  #univariate Cauchy
  #f(x) = 1 / (pi * scale * (1 + (d/scale)^2)) 
  
  #bivariate Cauchy from http://en.wikipedia.org/wiki/Cauchy_distribution
  #f(x) = (1/2*pi) * (scale / (p[,1] + p[,2] + scale^2)^1.5) 
  
  if(is.null(scale)) stop('Parameter scale needs to be specified')
  
  nx <- 1 + 2 * floor(d/rs[1])
  ny <- 1 + 2 * floor(d/rs[2])
  
  m <- matrix(ncol=nx, nrow=ny)
  m[ceiling(ny/2), ceiling(nx/2)] <- 1
  
  xr <- (nx * rs[1]) / 2
  yr <- (ny * rs[2]) / 2 
  r <- raster(m, xmn=-xr[1], xmx=xr[1], ymn=-yr[1], ymx=yr[1], crs='+proj=utm +zone=1')
  dist <- as.matrix(distance(r))
  #p <- xyFromCell(r, 1:ncell(r))^2   #use this and comment 'dist' if using bivar Cauchy
    
  if (length(scale) == 1) {
    if (!is.null(gamma)) stop('The parameter gamma has to be NULL when using a single Cauchy pdf')
    
    #using PDF
    #half-cauchy: f(x) = 2/pi * (scale/(r^2 + scale^2))  #where scale is the median disp distance here...
    m <- 2/pi * (scale/(dist^2 + scale^2))
        
    #using bivariate cauchy from http://en.wikipedia.org/wiki/Cauchy_distribution
    #m <- (1/2*pi) * (scale / (p[,1] + p[,2] + scale^2)^1.5) 
    #m <- matrix(m, ncol=nx, nrow=ny, byrow=TRUE)
    
    #sum of weights should add up to 1
    #m / sum(m)
    m
    
  }else{   #Cauchy mixture
    
    if (is.null(gamma)) stop('The parameter gamma has to be specified when using a mixture Cauchy pdf')
    if (gamma >= 1 | gamma <= 0) stop('The parameter gamma must range between (0-1)')
    scale1 <- scale[1]
    scale2 <- scale[2]
    
    m <- (gamma * 2/pi * (scale[1]/(dist^2 + scale[1]^2)) ) + ( (1-gamma) * 2/pi * (scale[2]/(dist^2 + scale[2]^2)) )
    
    #using bivariate cauchy from http://en.wikipedia.org/wiki/Cauchy_distribution
    #m <- gamma * ((1/2*pi) * (scale1 / (p[,1] + p[,2] + scale1^2)^1.5)) + (1-gamma) * ((1/2*pi) * (scale2 / (p[,1] + p[,2] + scale2^2)^1.5))
    #m <- matrix(m, ncol=nx, nrow=ny, byrow=TRUE)
    
    #sum of weights should add up to 1
    m / sum(m)    
    
  }
  
}  


gauss.weight <- function(rs, d, sigma=NULL, alpha=NULL) {
  
  #bivariate Gaussian filter, according to http://en.wikipedia.org/wiki/Gaussian_filter
  #f(x) = 1/(2*pi*sigma^2) * exp(-(p[,1]+p[,2])/(2*sigma^2))
  
  #two-parameter dispersal kernel as defined in Clark, J.S., 1998. Why trees migrate so fast [...]
  #f(x) = ( c / 2 * alpha * gamma(1/c) ) * exp( -abs(d/alpha)^c  ) 
  
  nx <- 1 + 2 * floor(d/rs[1])
  ny <- 1 + 2 * floor(d/rs[2])
  
  m <- matrix(ncol=nx, nrow=ny)
  
  xr <- (nx * rs[1]) / 2
  yr <- (ny * rs[2]) / 2 
  r <- raster(m, xmn=-xr[1], xmx=xr[1], ymn=-yr[1], ymx=yr[1], crs='+proj=utm +zone=1')
  dist <- as.matrix(distance(r))  
  #p <- xyFromCell(r, 1:ncell(r))^2   #use this and comment 'dist' if using bivar Gauss
      
  #using Clark, J.S., 1998
  if (length(alpha) != 1) stop('The parameter alpha must be a single number (average distance)')
  if (is.null(alpha)) stop('The parameter alpha must be specified (average distance)')
    
  c <- 2
  m <- ( c / 2 * alpha * gamma(1/c) ) * exp( -abs(dist/alpha)^c  ) 
  m <- matrix(m, ncol=nx, nrow=ny, byrow=TRUE)
  
  #using bivariate filter
  #if (length(sigma) != 1) stop('The parameter sigma must be a single number')
  #if (is.null(sigma)) stop('The parameter sigma must be specified')
  #m <- 1/(2*pi*sigma^2) * exp(-(p[,1]+p[,2])/(2*sigma^2))
  #m <- matrix(m, ncol=nx, nrow=ny, byrow=TRUE)
  
  #sum of weights should add up to 1
  m / sum(m)
  
}  
  

exp.weight <- function(rs, d, alpha) {
  
  #two-parameter dispersal kernel as defined in Clark, J.S., 1998. Why trees migrate so fast [...]
  #f(x) = ( c / 2 * alpha * gamma(1/c) ) * exp( -abs(d/alpha)^c  ) 
  
  nx <- 1 + 2 * floor(d/rs[1])
  ny <- 1 + 2 * floor(d/rs[2])
  
  m <- matrix(ncol=nx, nrow=ny)

  xr <- (nx * rs[1]) / 2
  yr <- (ny * rs[2]) / 2 
  r <- raster(m, xmn=-xr[1], xmx=xr[1], ymn=-yr[1], ymx=yr[1], crs='+proj=utm +zone=1')
  dist <- as.matrix(distance(r))  
 
  #using Clark, J.S., 1998
  if (length(alpha) != 1) stop('The parameter alpha must be a single number (average distance)')
  if (is.null(alpha)) stop('The parameter alpha must be specified (average distance)')
  
  c <- 1
  m <- ( c / 2 * alpha * gamma(1/c) ) * exp( -abs(dist/alpha)^c  ) 
  m <- matrix(m, ncol=nx, nrow=ny, byrow=TRUE)
  
  #sum of weights should add up to 1
  m / sum(m)
  
}  
  
 
Kernel.gen <- function(x, d, type = c('Cauchy', 'Cauchy Mixture', 'Exponential', 'Gauss'), scale = NULL, gamma = NULL, sigma = NULL, alpha = NULL){
  
  type <- match.arg(type)
  x <- res(x)
  
  if (type == 'Cauchy'){
    cauchy.weight(x, d, scale)
  }else if (type =='Cauchy Mixture'){
    cauchy.weight(x, d, scale, gamma)   
  }else if (type =='Exponential'){
    exp.weight(x, d, alpha = alpha)    
  }else if (type == 'Gauss'){
    gauss.weight(x, d, alpha = alpha)  
  }
  
}


WindKernel.gen <- function(x, d, mu, kappa){
  
  rs <- res(x)

  nx <- 1 + 2 * floor(d/rs[1])
  ny <- 1 + 2 * floor(d/rs[2])
  
  m <- matrix(ncol=nx, nrow=ny)
  
  xr <- (nx * rs[1]) / 2
  yr <- (ny * rs[2]) / 2 
  r <- raster(m, xmn=-xr[1], xmx=xr[1], ymn=-yr[1], ymx=yr[1], crs='+proj=utm +zone=1')
  r[ceiling(ny/2), ceiling(nx/2)] <- 1
  
  dir <- direction(r, from=TRUE)
  dir[ceiling(ny/2), ceiling(nx/2)] <- 0
  
  m <- as.matrix(dvm(dir, mu, kappa))
  
  #sum of weights should add up to 1
  #m / sum(m)
  m
  
}


distGen <- function(n, type = c('Cauchy', 'Cauchy Mixture', 'Exponential', 'Gauss'), mean = NULL, sd = NULL, scale = NULL, gamma=NULL) {
  
  type <- match.arg(type)
  
  if (type == 'Cauchy'){
    
    dist <- abs(rcauchy(n, scale = scale))
    
  }else if (type =='Cauchy Mixture'){
    
    if(length(scale) != 2) stop('The parameter scale must have two values, one for each component of the mixture!')
    if (is.null(gamma)) stop('The parameter gamma has to be specified when using a mixture Cauchy pdf!')
    if (gamma >= 1 | gamma <= 0) stop('The parameter gamma must range between (0-1)')
    f <- sample(1:2, n, replace=TRUE, prob=c(gamma, 1-gamma))
    dist <- abs(sapply(f, FUN=function(x) if(x==1) rcauchy(1, scale = scale[1]) else if (x==2) rcauchy(1, scale = scale[2])))
    
  }else if (type =='Exponential'){
    
    if (is.null(mean)) stop('The parameter mean must be specified when using an Exponential pdf!')
    dist <- rexp(n, rate = 1/mean)
    
  }else if (type == 'Gauss'){
    
    if (is.null(mean) | is.null(sd)) stop('The parameters mean and standard deviation (sd) must be specified when using an Gaussian pdf!')
    dist <- abs(rnorm(n, mean = mean, sd = sd))
    
  }
  
  return(dist)
  
}

#use Von Mises pdf to generate angular data
angleGen <- function(n, mean=0, k=0){
  
  #mean: mean direction in radians of the von Mises distribution
  #k: concentration parameter of the von Mises distribution. Small k --> large variance; large k --> small variance
  #if(k>0) theta <- rvm(n, mean, k) else 
  theta <- runif(n, 0, 360)
  return(theta)
  
}


findRowPos <- function(x, w){
  
  nrow <- nrow(HI.raster) - nearest.raster.point(x[1], x[2], w)[[1]] + 1
  return(nrow)
  
}


findColPos <- function(x, w){
  
  ncol <- nearest.raster.point(x[1], x[2], w)[[2]]
  return(ncol)
  
}


spores.count <- function(x, y){
  
  sum(x[1] == y[,1] & x[2] == y[,2])
  
}


#Within each infected cell (I > 0) draw random number of infections ~Poisson(lambda=rate of spore production) for each infected host. 
#Take SUM for total infections produced by each cell. 

new.infections.gen <- function(x, rate, clim=NULL){
  
    if(is.null(clim)) x <- rpois(x,lambda=rate) else x <- rpois(x,lambda=rate * clim)

    x <- sum(x)
    return(x)
  
}


myFocal.fun <- function(u, v, ngb){
  
  spores <- getValuesFocal(u, ngb=ngb)
  spores[is.na(spores)] <- 0
  prob <- c(t(v))
  
  new.I <- apply(spores, MARGIN = 1, FUN = function(x) ifelse(any(x>0), x %*% prob, 0))

  return(new.I)
  
}


radialDisp <- function(xyAOIPos, AOI_owin){
  
  #generate dispersal distances and angles from user-defined pdf
  distance <- distGen(n = sum(spores.lst), type = 'Cauchy', scale=20.57, gamma=0.9947)
  theta <- angleGen(n = sum(spores.lst))
  
  #set of new xy coords of dispersed spores
  x.new <- rep(x[spores.lst>0], spores.lst[spores.lst>0]) + sin(theta) * distance
  y.new <- rep(y[spores.lst>0], spores.lst[spores.lst>0]) + cos(theta) * distance
  
  #remove/ignore xy coords falling outside study area    
  condx <- which(x.new < extent(HI.raster)[1] | x.new > extent(HI.raster)[2])
  condy <- which(y.new < extent(HI.raster)[3] | y.new > extent(HI.raster)[4])
  
  x.new <- x.new[-unique(c(condx,condy))]
  y.new <- y.new[-unique(c(condx,condy))]
  
  #get index (position) of cells in study area containing new xy coords of dispersed spores
  yPos <- apply(cbind(x.new,y.new), 1, FUN=findRowPos, AOI_owin)
  xPos <- apply(cbind(x.new,y.new), 1, FUN=findColPos, AOI_owin)
  
  #create list of spore counts (spores in each cell after dispersal)
  new.spores.lst <- apply(xyAOIPos, 1, FUN=spores.count, cbind(xPos, yPos))
  
  return(new.spores.lst)  
    
}


SporeDisp <- function(x, S, I, W=NULL, rs, rtype=c('Cauchy', 'Cauchy Mixture', 'Exponential', 'Gauss'), mean = NULL, sd = NULL, scale = NULL, gamma=NULL, 
                      wtype=c('Uniform', 'VM'), wdir=c('N','NE','E','SE','S','SW','W','NW'), kappa=NULL){
  
  out <- list()
  
  rtype <- match.arg(rtype) #radial distribution type
  wtype <- match.arg(wtype) #wind distribution type
  
  #LOOP THROUGH EACH CELL of the input matrix 'x' (this should be the study area)
  for(row in 1:nrow(x)){
    
    for(col in 1:ncol(x)){
      
      if(x[row,col] > 0){  #if [row,col] > 0 disperse spores from it
        
        for(sp in 1:x[row,col]){
          
          ##GENERATE DISTANCES
          if (rtype=='Cauchy'){
            if(is.null(scale)) stop('Kernel of type "Cauchy" requires a scale parameter!')
            dist <- abs(rcauchy(1, scale = scale))
          }else if (rtype=='Cauchy Mixture'){
            if(length(scale) != 2) stop('The parameter scale must have two values, one for each component of the mixture!')
            if (is.null(gamma)) stop('The parameter gamma has to be specified when using a mixture Cauchy pdf!')
            if (gamma >= 1 | gamma <= 0) stop('The parameter gamma must range between (0-1)')
            f <- sample(1:2, 1, prob=c(gamma, 1-gamma))
            if(f == 1) dist <- abs(rcauchy(1, scale = scale[1])) else dist <- abs(rcauchy(1, scale = scale[2])) 
          }else if (rtype =='Exponential'){
            if (is.null(mean)) stop('The parameter mean must be specified when using an Exponential pdf!')
            dist <- rexp(1, rate = 1/mean)
          }else if (rtype == 'Gauss'){
            if (is.null(mean) | is.null(sd)) stop('The parameters mean and standard deviation (sd) must be specified when using an Gaussian pdf!')
            dist <- abs(rnorm(1, mean = mean, sd = sd))
          }
          
          ##GENERATE ANGLES
          if (wtype=='Uniform'){
            theta <- runif(1, -pi, pi)
          }else if(wtype=='VM'){
            if(is.null(kappa) | kappa <= 0) stop('kappa must be greater than zero!')
            #Check predominant windDir
            if(wdir == 'N') {
              theta <- rvm(1, mean = 0 * (pi/180), kappa)  #kappa=concentration
            }else if(wdir == 'NE'){
              theta <- rvm(1, mean = 45 * (pi/180), kappa)  #kappa=concentration
            }else if(wdir == 'E'){
              theta <- rvm(1, mean = 90 * (pi/180), kappa)  #kappa=concentration
            }else if(wdir == 'SE'){
              theta <- rvm(1, mean = 135 * (pi/180), kappa)  #kappa=concentration
            }else if(wdir == 'S'){
              theta <- rvm(1, mean = 180 * (pi/180), kappa)  #kappa=concentration
            }else if(wdir == 'SW'){
              theta <- rvm(1, mean = 225 * (pi/180), kappa)  #kappa=concentration
            }else if(wdir == 'W'){
              theta <- rvm(1, mean = 270 * (pi/180), kappa)  #kappa=concentration
            }else if(wdir == 'NW'){
              theta <- rvm(1, mean = 315 * (pi/180), kappa)  #kappa=concentration
            }            
          }
          
          #distance in cell counts
          #distc <- round(dist / res.win)
          
          row0 <- row - round((dist * cos(theta)) / rs)
          col0 <- col + round((dist * sin(theta)) / rs)
          
          if (row0 < 1 | row0 > nrow(x)) next     ## outside the region ##
          if (col0 < 1 | col0 > ncol(x)) next     ## outside the region ##
          
          if(S[row0, col0] > 0){  
            PropS <- S[row0, col0] / (S[row0, col0] + I[row0, col0])
            U <- runif(1)
            if(!is.null(W)) Prob <- PropS * W[row0, col0] else Prob <- PropS    #added weather

            if (U < Prob){
              I[row0, col0] <- I[row0, col0] + 1 
              S[row0, col0] <- S[row0, col0] - 1
            } 
          }#ENF IF
          
        }#END FOR LOOP OVER ALL SPORES WITHIN CURRENT CELL
        
      }#END IF
      
    }    
  }#END LOOP OVER ALL CELLS
  
  out$I <- I 
  out$S <- S
  
  return(out) 
  
}


SporeDisp2 <- function(x, S, I, W=NULL, rs, rtype=c('Cauchy', 'Cauchy Mixture', 'Exponential', 'Gauss'), mean = NULL, sd = NULL, scale = NULL, gamma=NULL, 
                       wtype=c('Uniform', 'VM'), wdir=c('N','NE','E','SE','S','SW','W','NW'), kappa=NULL){
  
  out <- list()
  
  rtype <- match.arg(rtype) #radial distribution type
  wtype <- match.arg(wtype) #wind distribution type
    
  #LOOP THROUGH EACH CELL of the input matrix 'x' (this should be the study area)
  for(row in 1:nrow(x)){
    
    for(col in 1:ncol(x)){
      
      if(x[row,col] > 0){  #if [row,col] > 0 disperse spores from it
          
        ##GENERATE DISTANCES ALL AT ONCE
        if (rtype=='Cauchy'){
          if(is.null(scale)) stop('Kernel of type "Cauchy" requires a scale parameter!')
          dist <- abs(rcauchy(x[row,col], scale = scale))
        }else if (rtype=='Cauchy Mixture'){
          if(length(scale) != 2) stop('The parameter scale must have two values, one for each component of the mixture!')
          if (is.null(gamma)) stop('The parameter gamma has to be specified when using a mixture Cauchy pdf!')
          if (gamma >= 1 | gamma <= 0) stop('The parameter gamma must range between (0-1)')
          f <- sample(1:2, x[row,col], prob=c(gamma, 1-gamma))
          dist <- ifelse(f == 1, abs(rcauchy(1, scale = scale[1])), abs(rcauchy(1, scale = scale[2])))
        }else if (rtype =='Exponential'){
          if (is.null(mean)) stop('The parameter mean must be specified when using an Exponential pdf!')
          dist <- rexp(x[row,col], rate = 1/mean)
        }else if (rtype == 'Gauss'){
          if (is.null(mean) | is.null(sd)) stop('The parameters mean and standard deviation (sd) must be specified when using an Gaussian pdf!')
          dist <- abs(rnorm(x[row,col], mean = mean, sd = sd))
        }
          
        ##GENERATE ANGLES ALL AT ONCE
        if (wtype=='Uniform'){
          theta <- runif(x[row,col], -pi, pi)
        }else if(wtype=='VM'){
          if(is.null(kappa) | kappa <= 0) stop('kappa must be greater than zero!')
          #Check predominant windDir
          if(wdir == 'N') {
            theta <- rvm(x[row,col], mean = 0 * (pi/180), kappa)  #kappa=concentration
          }else if(wdir == 'NE'){
            theta <- rvm(x[row,col], mean = 45 * (pi/180), kappa)  #kappa=concentration
          }else if(wdir == 'E'){
            theta <- rvm(x[row,col], mean = 90 * (pi/180), kappa)  #kappa=concentration
          }else if(wdir == 'SE'){
            theta <- rvm(x[row,col], mean = 135 * (pi/180), kappa)  #kappa=concentration
          }else if(wdir == 'S'){
            theta <- rvm(x[row,col], mean = 180 * (pi/180), kappa)  #kappa=concentration
          }else if(wdir == 'SW'){
            theta <- rvm(x[row,col], mean = 225 * (pi/180), kappa)  #kappa=concentration
          }else if(wdir == 'W'){
            theta <- rvm(x[row,col], mean = 270 * (pi/180), kappa)  #kappa=concentration
          }else if(wdir == 'NW'){
            theta <- rvm(x[row,col], mean = 315 * (pi/180), kappa)  #kappa=concentration
          }            
        }
          
        #distance in cell counts
        #distc <- round(dist / res.win)
          
        row0 <- row - round((dist * cos(theta)) / rs)
        col0 <- col + round((dist * sin(theta)) / rs)
          
        bool <- ifelse(row0 < 1 | row0 > nrow(x) | col0 < 1 | col0 > ncol(x), 1, 0)
        
        ## outside the region ##
        if(any(bool==1)){
          row0 <- row0[-which(bool==1)]
          col0 <- col0[-which(bool==1)]
        }
        
        #if all spores fall outside region continue loop to next cell
        if (length(row0) == 0 | length(col0) == 0) next 
        
        ##LOOP THROUGH ROW0 and COL0 TO ADD DISPERSED SPORES
        for(i in 1:length(row0)){
          if(S[row0[i], col0[i]] > 0){
            PropS <- round(S[row0[i], col0[i]] / (S[row0[i], col0[i]] + I[row0[i], col0[i]]), 2)
            U <- runif(1)
            if(!is.null(W)) Prob <- PropS * W[row0[i], col0[i]] else Prob <- PropS    #added weather
            if (U < Prob){    
              I[row0[i], col0[i]] <- I[row0[i], col0[i]] + 1 
              S[row0[i], col0[i]] <- S[row0[i], col0[i]] - 1
            }#ENF IF 
          }#ENF IF          
        }#END LOOP OVER DISP SPORES FROM CURRENT SOURCE CELL
          
      }#END IF
      
    }    
  }#END LOOP OVER ALL CELLS
  
  out$I <- I
  out$S <- S
  
  return(out) 
  
}


randomU <- function(x, y, z){
  
  out <- list()
  
  if(z > 0 & x > 0){
    
    currentPropS <- round( z / (z + y), 2)
    
    for (i in 1:x){    
      r <- runif(1)
      if (r < currentPropS){
        y <- y + 1 
        z <- z - 1
        currentPropS <- round( z / (z + y), 2)
      }
    }
    
  } 
  
  out$I <- y
  out$S <- z
  
  return(out)  
  
}


buffer.areas <- function(x, rs, width.row, width.col){
  
  z <- matrix(x, ncol=ncol(rs), byrow=T)
  buffer.col <- matrix(NA, nrow=nrow(rs), ncol=floor((width.col - 1)/2))
  z <- cbind(buffer.col, z, buffer.col)
  buffer.row <- matrix(NA, ncol=ncol(z), nrow=floor((width.row - 1)/2))
  z <- rbind(buffer.row, z, buffer.row)
  
  return(z)
  
}

kernel2D <- function(x, d){
  
  rs <- res(x)
  if(d <= res(x)[1] | d <= res(x)[2]) stop('d must be >= grain of study area')
  
  nx <- 1 + 2 * floor(d/rs[1])
  ny <- 1 + 2 * floor(d/rs[2])
  
  probKernel <- matrix(NA, ncol=nx, nrow=ny)
  probKernel[ceiling(ny/2), ceiling(nx/2)] <- 1
  
  return(probKernel)
  
}


kernel2DToDistance <- function(x, r){
  
  rs <- res(r)
  
  xr <- (ncol(x) * rs[1]) / 2
  yr <- (nrow(x) * rs[2]) / 2 
  rast <- raster(x, xmn=-xr[1], xmx=xr[1], ymn=-yr[1], ymx=yr[1], crs='+proj=utm +zone=1')
  
  z <- as.matrix(distance(rast))
  
  return(z)  
  
}


dispFunction <- function(bkgr, dst, type=c('Cauchy','Cauchy Mixture','Gauss',"Exponential"), c=NULL, scale=NULL, gamma=NULL, sd=NULL, rate=NULL){
  
  type <- match.arg(type)
  
  hw_col <- floor(ncol(dst)/2)
  hw_row <- floor(nrow(dst)/2)
  
  if(is.null(c)) stop('c must be specified in the dispersal function!')
  if(c != 2 & c > 1) stop('c must be either <= 1 or = 2')  
  
  #Start moving the kernel over all cells of study area. BUT only run dispersal on top of source cells (cells with bkgr)
  for(row in 1:nrow(bkgr)){
    
    for(col in 1:ncol(bkgr)){
      
      if(!is.na(bkgr[row,col]) & bkgr[row,col] > 0) {
        
        #************************
        #RADIAL DISPERSAL
        #************************
        
        #sample a value for alpha using a radial dispersal pdf (this step is necessary to add stochasticity)
        if(type == 'Cauchy'){
          if(c > 1) stop('for type == Cauchy the value of c must be < 1')
          #sample alpha from chosen distribution of interest
          alpha <- abs(rcauchy(1, location=0, scale=scale))          
        }else if(type == 'Cauchy Mixture'){
          if(c > 1) stop('for type == "Cauchy" the value of c must be < 1')
          if(length(scale) != 2) stop('The parameter scale must have two values, one for each component of the mixture!')
          if (gamma >= 1 | gamma <= 0) stop('The parameter gamma must range between (0-1)')
          f <- sample(1:2, 1, prob=c(gamma, 1-gamma))
          if(f==1) alpha <- abs(rcauchy(1, location=0, scale = scale[1])) else if (f==2) alpha <- abs(rcauchy(1, location=0, scale = scale[2]))
        }else if(type == 'Gauss'){
          if(c != 2) stop('for type = "Gauss" the value of c must be = 2')
          alpha <- abs(rnorm(1, location=0, sd=sd)) 
        }else if(type == 'Exponential'){
          if(c != 1) stop('for type = "Exponential" the value of c must be = 1')
          alpha <- rexp(1, rate=1/rate) 
        }
        
        #now we can calculate the kernel weights using the value of alpha and set of distances from central source cell
        kernelWt <- ( c / 2 * alpha * gamma(1/c) ) * exp( -abs(dst/alpha)^c  )  ##CHECK IF THIS NEEDS TO BE CONVERTED INTO ANOTHER LOOP IN C++
        
        #sum of weights should add up to 1
        Wt <- kernelWt/sum(kernelWt)
         
        #MULTIPLY THE SPORES IN THE SOURCE CELL BY EACH VALUE IN DISP KERNEL
        #spores_disp <- round(bkgr[row,col] * Wt)
        spores_disp <- bkgr[row,col] * Wt
        
        bkgr[(row-hw_row):(row+hw_row), (col-hw_col):(col+hw_col)] <- bkgr[(row-hw_row):(row+hw_row), (col-hw_col):(col+hw_col)] + spores_disp
        
      } #END IF   
      
    }  
  } #END LOOP OVER EACH CELL INSIDE STUDY BUFFERED STUDY AREA
  
  return(round(bkgr))
}


infectFun <- function(bkgr, S, I){
  
  out <- list()
  
  #LOOP OVER EACH CELL OF STUDY AREA AND FOR EACH NON-ZERO SPORE CELL (dispersed spores) RUN THE PROBABILISTIC 
  #LOOP TO SEE HOW MANY SPORES TURN INTO AN INFECTION
  for(row in 1:nrow(bkgr)){
    for(col in 1:ncol(bkgr)){
      
      if(!is.na(S[row,col]) & S[row,col] > 0 & bkgr[row,col] > 0){  
        currentPropS <- round(S[row,col] / (S[row,col] + I[row,col]), 2)
        for(spNbr in 1:bkgr[row,col]){
          U <- runif(1)
          if (U < currentPropS){
            I[row,col] <- I[row,col] + 1 
            S[row,col] <- S[row,col] - 1
            currentPropS <- round( S[row,col] / (S[row,col] + I[row,col]), 2)
          } 
        } #END LOOP OVER EACH SPORE IN CELL [row,col]            
      } #ENF IF
      
    }
  }
  
  out$S <- S
  out$I <- I
  return(out)
  
}


climGen <- function(ls, start, end, Wrank, scn=NA){
  
  #CASE 1: no future weather needed
  if(is.na(scn) & end <= last_yr){
    ls <- grep(paste(as.character(seq(start,end)), collapse="|"), ls, value=TRUE)  #use only the raster files matching the years of interest
    ls.stack <- stack(ls)  
    #CASE 2: RANDOM future weather scenario  
  }else if(scn == 'random'){
    if(end <= last_yr) stop('you specified a future weather scenario BUT the end year is not a future year!')
    yrs <- sample(weather_rank[,1], size = end - last_yr, replace = T)
    ls_1 <- grep(paste(as.character(seq(start,end)), collapse="|"), ls, value=TRUE)
    ls_2 <- unlist(lapply(yrs, FUN=function(x){grep(as.character(x), ls, value=TRUE)}))
    ls <- c(ls_1, ls_2)
    ls.stack <- stack(ls) 
    #CASE 3: UPPER 50% (favorable) future weather scenario
  }else if(scn == 'favorable'){
    if(end <= last_yr) stop('you specified a future weather scenario BUT the end year is not a future year!')
    yrs <- sample(weather_rank[1:9,1], size = end - last_yr, replace = T)
    ls_1 <- grep(paste(as.character(seq(start,end)), collapse="|"), ls, value=TRUE)
    ls_2 <- unlist(lapply(yrs, FUN=function(x){grep(as.character(x), ls, value=TRUE)}))
    ls <- c(ls_1, ls_2)
    ls.stack <- stack(ls) 
    #CASE 4: LOWER 50% (unfavorable) future weather scenario
  }else{
    if(end <= last_yr) stop('you specified a future weather scenario BUT the end year is not a future year!')
    yrs <- sample(weather_rank[10:18, 1], size = end - last_yr, replace = T)
    ls_1 <- grep(paste(as.character(seq(start,end)), collapse="|"), ls, value=TRUE)
    ls_2 <- unlist(lapply(yrs, FUN=function(x){grep(as.character(x), ls, value=TRUE)}))
    ls <- c(ls_1, ls_2)
    ls.stack <- stack(ls) 
  }
  
  return(ls.stack)
  
}


