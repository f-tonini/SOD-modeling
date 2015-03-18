#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar



//Within each infected cell (I > 0) draw random number of infections ~Poisson(lambda=rate of spore production) for each infected host. 
//Take SUM for total infections produced by each cell. 

// [[Rcpp::export]]
IntegerMatrix SporeGenCpp(IntegerMatrix I, NumericMatrix W, double rate){
  
  
  // internal variables
  int nrow = I.nrow(); 
  int ncol = I.ncol();
  
  IntegerMatrix SP = clone<IntegerMatrix>(I);
  //Function rpois("rpois");
  
  RNGScope scope;
  
  // LOOP THROUGH EACH INFECTED CELL AND GENERATE AMOUNT OF SPORES
  for (int row = 0; row < nrow; row++) {
    for (int col = 0; col < ncol; col++){
      
      if (I(row, col) > 0){  //if infected > 0, generate spores proportional to production rate * weather suitability

        double lambda = rate * W(row, col);
        NumericVector inf = rpois(I(row, col), lambda);
        int s = sum(inf);
        SP(row, col) = s;

      }
   
    }
  }
  
  
  return SP;


}

// [[Rcpp::export]]
List SporeDispCpp(IntegerMatrix x, IntegerMatrix S, IntegerMatrix I, NumericMatrix W,   //use different name than the functions in myfunctions_SOD.r
                double rs, String rtype, String wtype, 
                String wdir=NA_STRING,
                double scale1=NA_REAL, int kappa=NA_REAL, double scale2=NA_REAL,
                double mean=NA_REAL, double sd=NA_REAL,
                double gamma=NA_REAL){

  // internal variables //
  int nrow = x.nrow(); 
  int ncol = x.ncol();
  int row0;
  int col0;

  double dist = 0;
  double theta = 0;
  double PropS;

  //for Rcpp random numbers
  RNGScope scope;
  
  //Function rcauchy("rcauchy");
  //Function rexp("rexp");
  //Function rnorm("rnorm");
  //Function runif("runif");
  Function rvm("rvm");
  Function sample("sample");

  //LOOP THROUGH EACH CELL of the input matrix 'x' (this should be the study area)
  for (int row = 0; row < nrow; row++) {
    for (int col = 0; col < ncol; col++){
      
      if(x(row,col) > 0){  //if spores in cell (row,col) > 0, disperse
        
        for(int sp = 1; (sp <= x(row,col)); sp++){
          
          //GENERATE DISTANCES
          if (rtype == "Cauchy"){

            dist = abs(R::rcauchy(0, scale1));
          
          }else if (rtype == "Cauchy Mixture"){
            
            if (gamma >= 1 || gamma <= 0) stop("The parameter gamma must range between (0-1)");
            
            //TO DO: does it have to be NumericVector even for size=1 ? ******************
            NumericVector fv = sample(Range(1, 2), 1, false, (gamma, 1-gamma));
            int f = fv[0];
            if(f == 1) 
              dist = abs(R::rcauchy(0, scale1));
            else if (f==2) 
              dist = abs(R::rcauchy(0, scale2));

          }else if (rtype == "Exponential"){
          
            if (mean <= 0) stop("The parameter mean must be greater than zero!");
            dist = abs(R::rexp(1/mean));

          }else if (rtype == "Gauss"){
            
            if (mean <= 0 || sd <= 0) stop("The parameters mean and standard deviation (sd) must be greater than zero!");
            dist = abs(R::rnorm(mean, sd));

          }
                    
          //GENERATE ANGLES
          if (wtype=="Uniform"){

            theta = R::runif(-PI, PI);
            
          }else if(wtype=="VM"){
            
            if(kappa <= 0) 
              stop("kappa must be greater than zero!");
            
            if (wdir == "N") 
              theta = as<double>(rvm(1, 0 * (PI/180), kappa));  // kappa=concentration
            else if (wdir == "NE")
              theta = as<double>(rvm(1, 45 * (PI/180), kappa));  // kappa=concentration
            else if(wdir == "E")
              theta = as<double>(rvm(1, 90 * (PI/180), kappa));  // kappa=concentration
            else if(wdir == "SE")
              theta = as<double>(rvm(1, 135 * (PI/180), kappa));  // kappa=concentration
            else if(wdir == "S")
              theta = as<double>(rvm(1, 180 * (PI/180), kappa));  // kappa=concentration
            else if(wdir == "SW")
              theta = as<double>(rvm(1, 225 * (PI/180), kappa));  // kappa=concentration
            else if(wdir == "W")
              theta = as<double>(rvm(1, 270 * (PI/180), kappa));  // kappa=concentration
            else if(wdir == "NW")
              theta = as<double>(rvm(1, 315 * (PI/180), kappa));  // kappa=concentration

          }
          
          
          row0 = row - round((dist * cos(theta)) / rs);
          col0 = col + round((dist * sin(theta)) / rs);
          
          
          if (row0 < 1 || row0 >= nrow) continue;     //outside the region
          if (col0 < 1 || col0 >= ncol) continue;     //outside the region
          
          
          if(S(row0, col0) > 0){  
            PropS = double(S(row0, col0)) / (S(row0, col0) + I(row0, col0));
            double U = R::runif(0,1);
            double Prob = PropS * W(row0, col0);

            if (U < Prob){   //weather suitability affects prob success!  
              I(row0, col0) = I(row0, col0) + 1; 
              S(row0, col0) = S(row0, col0) - 1;
            } 
          }//ENF IF
          
        
        }//END LOOP OVER ALL SPORES IN CURRENT CELL GRID
       
       
      }//END IF S > 0  
      
    }   
  }//END LOOP OVER ALL GRID CELLS

  //return List::create(Named("I")=I, Named("S")=S);
  return List::create(
    _["S"] = S, 
    _["I"] = I
  );
  
}
