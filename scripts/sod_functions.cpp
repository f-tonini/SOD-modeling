#include <Rcpp.h>
using namespace Rcpp;

/**
type=c('Cauchy','Cauchy Mixture','Gauss',"Exponential")
# bkgr, S, I matrix, integer
# kernel matrix, integer
# c float
# scale float or list of two floats
# gamma, sd, rate fl
 */
// [[Rcpp::export]]
List dispFunction(NumericMatrix bkgr, NumericMatrix S, NumericMatrix I,
    NumericMatrix kernel, String type,
    float c, NumericVector scale,
    float gamma, float sd, float rate) {
  
  int hw_col = floor(kernel.ncol() / 2);
  int hw_row = floor(kernel.nrow() / 2);
  
  //if(is.null(c)) stop("c must be specified in the dispersal function!");
  if(c != 2 & c > 1) stop("c must be either <= 1 or = 2");
  
  // Start moving the kernel over all cells of study area. BUT only run dispersal on top of source cells (cells with bkgr)
  for(int row = 1; row < bkgr.nrow(); row++){
    
    for(int col = 1; col < bkgr.ncol(); col++){
      
      if(!NumericMatrix::is_na(bkgr[row, col]) & bkgr[row, col] > 0) {
        
        float alpha;

		Function rcauchy("rcauchy");

        if(type == "Cauchy"){
          if(c > 1)
			stop("for type == Cauchy the value of c must be < 1");
          // sample alpha from chosen distribution of interest
          // TODO: is [] and the assignment safe here?
          alpha = abs(rcauchy(1, 0, scale))[0];
        }else if(type == "Cauchy Mixture"){
          if(c > 1)
			stop("for type == 'Cauchy' the value of c must be < 1");
          if(scale.size() != 2)
			stop("The parameter scale must have two values, one for each component of the mixture!");
          if (gamma >= 1 | gamma <= 0)
			stop("The parameter gamma must range between (0-1)");
		  Function sample("sample");
		  // TODO: do this in some more effective way
          NumericVector fv = sample(Range(1, 2), 1, false, (gamma, 1-gamma));
          int f = fv[0];
          if(f==1)
			alpha = abs(rcauchy(1, 0, scale[1]))[0];
		  else if (f==2)
		    alpha = abs(rcauchy(1, 0, scale[2]))[0];
        }else if(type == "Gauss"){
          if(c != 2)
			stop("for type = 'Gauss' the value of c must be = 2");
          Function rnorm("rnorm");
          alpha = abs(rnorm(1, 0, sd))[0];
        }else if(type == "Exponential"){
          if(c != 1)
			stop("for type = 'Exponential' the value of c must be = 1");
          alpha = rexp(1, rate=1/rate)[0];
        }
		Function gamma_function("gamma");
		// CHECK IF THIS NEEDS TO BE CONVERTED INTO ANOTHER LOOP IN C++
		// what is dist?
        kernel = ( c / 2 * alpha * as<float>(gamma_function(1/c)) );// * exp( -pow(abs(dist/alpha), c) );
        // sum of weights should add up to 1
        // TODO: cppFunction("NumericVector fun(NumericVector a){ a = a / sum(a); return a;}")
        kernel = kernel;//sum(kernel);
        
        // MULTIPLY THE SPORES IN THE SOURCE CELL BY EACH VALUE IN DISP KERNEL
        Function round("round");
        NumericVector spores_disp = round(bkgr[row,col] * kernel);
        
        // loop thru values of infected and susceptibles falling within 2D kernel
        for(int i = row-hw_row; i <= row+hw_row; i++){
          for(int j = col-hw_col; j <= col+hw_col; j++){
            
            // LOOP OVER EACH CELL WITHIN 2D KERNEL (SKIP CELL IF THERE ARE NO Susceptibles & spores PRESENT IN IT)
            if(!NumericMatrix::is_na(S[i,j]) && S[i,j] > 0 && spores_disp[i,j] > 0){
              float currentPropS = as<float>(round(S[i,j] / (S[i,j] + I[i,j]), 2));
              for (int spNbr = 1;
                   spNbr < spores_disp[i-(row-hw_row) + 1,j-(col-hw_col) + 1];
                   spNbr++){
                Function runif("runif");
                float U = as<NumericVector>(runif(1))[0];
                if (U < currentPropS){
                  I[i,j] = I[i,j] + 1;
                  S[i,j] = S[i,j] - 1;
                  currentPropS = as<float>(round( S[i,j] / (S[i,j] + I[i,j]), 2));
                } 
              } // END LOOP OVER EACH SPORE DISPERSED IN CELL i,j wITHIN 2D KERNEL             
            } // ENF IF
          } // END LOOP OVER EACH CELL FALLING WITHIN 2D KERNEL        
        }
      } // END IF
    }  
  } // END LOOP OVER EACH CELL INSIDE STUDY BUFFERED STUDY AREA
  return List::create(
    _["S"] = S, 
    _["I"] = I
  );
}
