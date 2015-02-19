// [[Rcpp::depends(RcppArmadillo)]]

// the following header also include Rcpp and Armadillo headers
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
SEXP SporeDisp(SEXP SPORES, SEXP SUSC, SEXP INF, String type, float res_win, float mean=NA_REAL, float sd=NA_REAL, float scale1=NA_REAL, float scale2=NA_REAL, float gamma=NA_REAL){

  // input parameters //
  IntegerMatrix spores (SPORES);
  IntegerMatrix S (SUSC);
  IntegerMatrix I (INF);
  
  // internal variables //
  IntegerMatrix sporedisp = clone<IntegerMatrix>(spores);
  
  int nrow = sporedisp.nrow(); 
  int ncol = sporedisp.ncol();
  int dist;
  int row0;
  int col0;
  
  float theta;
  float PropS;
  
  RNGScope scope;
  
  for (int row = 0; row < nrow; row++) {
    for (int col = 0; col < ncol; col++){
      
      if(spores(row,col) > 0){
        
        for(int sp = 1; (sp = spores(row,col)); sp++){
          
          if (type == "Cauchy"){
            dist = abs(R::rcauchy(0, scale1));
          }else if (type =='Cauchy Mixture'){
            if( is_na(scale1) || is_na(scale2) ) stop('The parameter scale must have two values, one for each component of the mixture!');
            if (is_na(gamma)) stop('The parameter gamma has to be specified when using a mixture Cauchy pdf!');
            if (gamma >= 1 | gamma <= 0) stop('The parameter gamma must range between (0-1)');
            NumericVector x = NumericVector::create(1,2);
            NumericVector prob = NumericVector::create(gamma,1-gamma);
            int size = 1;           
            int f = RcppArmadillo::sample(x, size=size, prob=prob);
            if(f == 1) dist <- abs(rcauchy(1, scale = scale[1])) else dist <- abs(rcauchy(1, scale = scale[2])) 
          }else if (type =='Exponential'){
            if (is.null(mean)) stop('The parameter mean must be specified when using an Exponential pdf!')
            dist = abs(R::rexp(0, 1/mean));
          }else if (type == 'Gauss'){
            if (is.null(mean) | is.null(sd)) stop('The parameters mean and standard deviation (sd) must be specified when using an Gaussian pdf!')
            dist <- abs(rnorm(1, mean = mean, sd = sd))
          }
                    
          //angle
          theta = R::runif(-PI, PI);
                    
          row0 = row - floor(((dist * cos(theta)) + (res_win/2)) / res_win);
          col0 = col + floor(((dist * sin(theta)) + (res_win/2)) / res_win);
          
          if (row0 < 1 || row0 >= nrow) continue;     //outside the region
          if (col0 < 1 || col0 >= ncol) continue;     //outside the region
          
          if(S(row0, col0) > 0){  
            PropS = S(row0, col0) / (S(row0, col0) + I(row0, col0));
            float U = rand();
            if (U < PropS){
              I(row0, col0) =+ 1; 
              S(row0, col0) =- 1;
            } 
          }//ENF IF
          
        
        }//END LOOP OVER ALL SPORES IN CURRENT CELL GRID
       
       
      }//END IF S > 0  
      
    }   
  }//END LOOP OVER ALL GRID CELLS


/*
  out$I <- infected 
  out$S <- susceptible
  
  return(out) 
*/

  return List::create(Named("I")=I, Named("S")=S);
  
  
}



