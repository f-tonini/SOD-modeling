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

//spore generation function for multi-host patho-system with HOST competency weights
// [[Rcpp::export]]
IntegerMatrix SporeGenCpp_MH(IntegerMatrix I_UM, IntegerMatrix I_LD, IntegerMatrix I_AC,
                             IntegerMatrix I_AR, IntegerMatrix I_AE, IntegerMatrix I_PS,
                             IntegerMatrix I_SE, List HW, NumericMatrix W, 
                             double rate){
  
  
  // internal variables
  int nrow = I_UM.nrow(); 
  int ncol = I_UM.ncol();
  NumericVector inf;
  
  NumericVector hw_um_v = as<NumericVector>(HW["UMCA"]);
  NumericVector hw_ld_v = as<NumericVector>(HW["LIDE3"]);
  NumericVector hw_ac_v = as<NumericVector>(HW["ACMA3"]);
  NumericVector hw_ar_v = as<NumericVector>(HW["ARME"]);
  NumericVector hw_ae_v = as<NumericVector>(HW["AECA"]);
  NumericVector hw_ps_v = as<NumericVector>(HW["PSME"]);
  NumericVector hw_se_v = as<NumericVector>(HW["SESE3"]);
  
  // initialize a new matrix of zeroes by cloning one of the infected matrices
  IntegerMatrix SP = clone<IntegerMatrix>(I_UM);
  
  RNGScope scope;
  
  //LOOP THROUGH EACH CELL OF LANDSCAPE MATRIX
  for (int row = 0; row < nrow; row++) {
    for (int col = 0; col < ncol; col++){
      
      int s = 0;
      
      // IF IN THE CURRENT CELLS THERE ARE SOME INFECTED TREES, THEN GENERATE SPORES 
      // PROPORTIONALLY TO THE HOST WEIGHT AND LOCAL WEATHER SUITABILITY
      if (I_UM(row, col) > 0){ 
        double lambda = rate * W(row, col) * hw_um_v[0];
        inf = rpois(I_UM(row, col), lambda);
        s = s + sum(inf);
      }
      
      if (I_LD(row, col) > 0){
        double lambda = rate * W(row, col) * hw_ld_v[0];
        inf = rpois(I_LD(row, col), lambda);
        s = s + sum(inf);
      }
      
      if (I_AC(row, col) > 0){
        double lambda = rate * W(row, col) * hw_ac_v[0];
        inf = rpois(I_AC(row, col), lambda);
        s = s + sum(inf);
      }
      
      if (I_AR(row, col) > 0){
        double lambda = rate * W(row, col) * hw_ar_v[0];
        inf = rpois(I_AR(row, col), lambda);
        s = s + sum(inf);
      }
      
      if (I_AE(row, col) > 0){
        double lambda = rate * W(row, col) * hw_ae_v[0];
        inf = rpois(I_AE(row, col), lambda);
        s = s + sum(inf);
      }
      
      if (I_PS(row, col) > 0){
        double lambda = rate * W(row, col) * hw_ps_v[0];
        inf = rpois(I_PS(row, col), lambda);
        s = s + sum(inf);
      }
      
      if (I_SE(row, col) > 0){
        double lambda = rate * W(row, col) * hw_se_v[0];
        inf = rpois(I_SE(row, col), lambda);
        s = s + sum(inf);
      }
      
      SP(row, col) = s;
      
    }
  }//END OF LOOP
  
  
  return SP;
  
  
}

// [[Rcpp::export]]
List SporeDispCpp(IntegerMatrix x, IntegerMatrix S, IntegerMatrix I, NumericMatrix W,   //use different name than the functions in myfunctions_SOD.r
                  double rs, String rtype, double scale1,
                  double scale2=NA_REAL,  //default values
                  double gamma=NA_REAL){  //default values
  
  // internal variables //
  int nrow = x.nrow(); 
  int ncol = x.ncol();
  int row0;
  int col0;
  
  double dist;
  double theta;
  double PropS;
  
  //for Rcpp random numbers
  RNGScope scope;
  
  Function sample("sample");
  //Function rcauchy("rcauchy");  
  
  //LOOP THROUGH EACH CELL of the input matrix 'x' (this should be the study area)
  for (int row = 0; row < nrow; row++) {
    for (int col = 0; col < ncol; col++){
      
      if(x(row,col) > 0){  //if spores in cell (row,col) > 0, disperse
        
        for(int sp = 1; (sp <= x(row,col)); sp++){
          
          //GENERATE DISTANCES:
          if (rtype == "Cauchy") 
            dist = abs(R::rcauchy(0, scale1));
          else if (rtype == "Cauchy Mixture"){
            if (gamma >= 1 || gamma <= 0) stop("The parameter gamma must range between (0-1)");
            NumericVector fv = sample(Range(1, 2), 1, false, NumericVector::create(gamma, 1-gamma));
            int f = fv[0];
            if(f == 1) 
              dist = abs(R::rcauchy(0, scale1));
            else
              dist = abs(R::rcauchy(0, scale2));
          }
          else 
            stop("The parameter rtype must be set to either 'Cauchy' or 'Cauchy Mixture'");
          
          //GENERATE ANGLES (using Uniform distribution):
          theta = R::runif(-PI, PI);
          
          //calculate new row and col position for the dispersed spore unit (using dist and theta)
          row0 = row - round((dist * cos(theta)) / rs);
          col0 = col + round((dist * sin(theta)) / rs);
          
          
          if (row0 < 0 || row0 >= nrow) continue;     //outside the region
          if (col0 < 0 || col0 >= ncol) continue;     //outside the region
          
          //if susceptibles are present in current cell, calculate prob of infection
          if(S(row0, col0) > 0){  
            PropS = double(S(row0, col0)) / (S(row0, col0) + I(row0, col0));
            double U = R::runif(0,1);
            double Prob = PropS * W(row0, col0); //weather suitability affects prob success! 
            
            //if U < Prob then one unit becomes infected
            if (U < Prob){    
              I(row0, col0) = I(row0, col0) + 1; //update infected
              S(row0, col0) = S(row0, col0) - 1; //update susceptible
            } 
          }//ENF IF
          
          
        }//END LOOP OVER ALL SPORES IN CURRENT CELL GRID
        
        
      }//END IF  
      
    }   
  }//END LOOP OVER ALL GRID CELLS
  
  //return List::create(Named("I")=I, Named("S")=S);
  return List::create(
    _["S"] = S, 
    _["I"] = I
  );
  
}

// [[Rcpp::export]]
List SporeDispCppWind(IntegerMatrix x, IntegerMatrix S, IntegerMatrix I, NumericMatrix W,   //use different name than the functions in myfunctions_SOD.r
                      double rs, String rtype, double scale1, 
                      String wdir, int kappa,
                      double scale2=NA_REAL,  //default values
                      double gamma=NA_REAL){  //default values
  
  // internal variables //
  int nrow = x.nrow(); 
  int ncol = x.ncol();
  int row0;
  int col0;
  
  double dist;
  double theta;
  double PropS;
  
  //for Rcpp random numbers
  RNGScope scope;
  
  //Function rcauchy("rcauchy");
  Function rvm("rvm");
  Function sample("sample");
  
  //LOOP THROUGH EACH CELL of the input matrix 'x' (this should be the study area)
  for (int row = 0; row < nrow; row++) {
    for (int col = 0; col < ncol; col++){
      
      if(x(row,col) > 0){  //if spores in cell (row,col) > 0, disperse
        
        for(int sp = 1; (sp <= x(row,col)); sp++){
          
          //GENERATE DISTANCES:
          if (rtype == "Cauchy") 
            dist = abs(R::rcauchy(0, scale1));
          else if (rtype == "Cauchy Mixture"){
            if (gamma >= 1 || gamma <= 0) stop("The parameter gamma must range between (0-1)");
            NumericVector fv = sample(Range(1, 2), 1, false, NumericVector::create(gamma, 1-gamma));
            int f = fv[0];
            if(f == 1) 
              dist = abs(R::rcauchy(0, scale1));
            else 
              dist = abs(R::rcauchy(0, scale2));
          }
          else 
            stop("The parameter rtype must be set to either 'Cauchy' or 'Cauchy Mixture'");
          
          //GENERATE ANGLES (using Von Mises distribution):
          if(kappa <= 0)  // kappa=concentration
            stop("kappa must be greater than zero!");
          
          //predominant wind dir
          if (wdir == "N") 
            theta = as<double>(rvm(1, 0 * (PI/180), kappa));  
          else if (wdir == "NE")
            theta = as<double>(rvm(1, 45 * (PI/180), kappa));  
          else if(wdir == "E")
            theta = as<double>(rvm(1, 90 * (PI/180), kappa));  
          else if(wdir == "SE")
            theta = as<double>(rvm(1, 135 * (PI/180), kappa));  
          else if(wdir == "S")
            theta = as<double>(rvm(1, 180 * (PI/180), kappa));  
          else if(wdir == "SW")
            theta = as<double>(rvm(1, 225 * (PI/180), kappa));  
          else if(wdir == "W")
            theta = as<double>(rvm(1, 270 * (PI/180), kappa));  
          else
            theta = as<double>(rvm(1, 315 * (PI/180), kappa));
          
          
          //calculate new row and col position for the dispersed spore unit (using dist and theta)
          row0 = row - round((dist * cos(theta)) / rs);
          col0 = col + round((dist * sin(theta)) / rs);
          
          
          if (row0 < 0 || row0 >= nrow) continue;     //outside the region
          if (col0 < 0 || col0 >= ncol) continue;     //outside the region
          
          //if susceptibles are present in current cell, calculate prob of infection
          if(S(row0, col0) > 0){  
            PropS = double(S(row0, col0)) / (S(row0, col0) + I(row0, col0));
            double U = R::runif(0,1);
            double Prob = PropS * W(row0, col0); //weather suitability affects prob success! 
            
            //if U < Prob then one unit becomes infected
            if (U < Prob){    
              I(row0, col0) = I(row0, col0) + 1; //update infected
              S(row0, col0) = S(row0, col0) - 1; //update susceptible
            } 
          }//ENF IF
          
          
        }//END LOOP OVER ALL SPORES IN CURRENT CELL GRID
        
        
      }//END IF 
      
    }   
  }//END LOOP OVER ALL GRID CELLS
  
  //return List::create(Named("I")=I, Named("S")=S);
  return List::create(
    _["S"] = S, 
    _["I"] = I
  );
  
}

// [[Rcpp::export]]
List SporeDispCpp_MH(IntegerMatrix x, 
                     IntegerMatrix S_UM, IntegerMatrix S_LD, //SUSCEPTIBLE
                     IntegerMatrix S_AC, IntegerMatrix S_AR, 
                     IntegerMatrix S_AE, IntegerMatrix S_PS, 
                     IntegerMatrix S_SE, IntegerMatrix S_OK, 
                     IntegerMatrix I_UM, IntegerMatrix I_LD, //INFECTED
                     IntegerMatrix I_AC, IntegerMatrix I_AR, 
                     IntegerMatrix I_AE, IntegerMatrix I_PS, 
                     IntegerMatrix I_SE, IntegerMatrix I_OK, 
                     IntegerMatrix N_LVE,
                     NumericMatrix W,   //use different name than the functions in myfunctions_SOD.r
                     double rs, String rtype, double scale1, 
                     double scale2=NA_REAL,  //default values
                     double gamma=NA_REAL)	//default values
{  
  
  // internal variables //
  int nrow = x.nrow(); 
  int ncol = x.ncol();
  int row0;
  int col0;
  
  double dist;
  double theta;
  double PropS;
  
  //for Rcpp random numbers
  RNGScope scope;
  
  Function sample("sample");
  //Function rcauchy("rcauchy");  
  
  //LOOP THROUGH EACH CELL of the input matrix 'x' (this should be the study area)
  for (int row = 0; row < nrow; row++) {
    for (int col = 0; col < ncol; col++){
      
      if(x(row,col) > 0){  //if spores in cell (row,col) > 0, disperse
        
        for(int sp = 1; (sp <= x(row,col)); sp++){
          
          //GENERATE DISTANCES:
          if (rtype == "Cauchy") 
            dist = abs(R::rcauchy(0, scale1));
          else if (rtype == "Cauchy Mixture"){
            if (gamma >= 1 || gamma <= 0) stop("The parameter gamma must range between (0-1)");
            NumericVector fv = sample(Range(1, 2), 1, false, NumericVector::create(gamma, 1-gamma));
            int f = fv[0];
            if(f == 1) 
              dist = abs(R::rcauchy(0, scale1));
            else 
              dist = abs(R::rcauchy(0, scale2));
          }
          else stop("The parameter rtype must be set to either 'Cauchy' or 'Cauchy Mixture'");
          
          //GENERATE ANGLES (using Uniform distribution):
          theta = R::runif(-PI, PI);
          
          //calculate new row and col position for the dispersed spore unit (using dist and theta)
          row0 = row - round((dist * cos(theta)) / rs);
          col0 = col + round((dist * sin(theta)) / rs);
          
          if (row0 < 0 || row0 >= nrow) continue;     //outside the region
          if (col0 < 0 || col0 >= ncol) continue;     //outside the region
          
          //if distance is within same pixel challenge all SOD hosts, otherwise challenge UMCA only
          if (row0 == row && col0 == col){
            
            //if susceptible hosts are present in cell, calculate prob of infection
            if(S_UM(row0, col0) > 0 || S_OK(row0, col0) > 0 || S_LD(row0, col0) > 0 || S_AC(row0, col0) > 0 ||
               S_AR(row0, col0) > 0 || S_AE(row0, col0) > 0 || S_PS(row0, col0) > 0 || S_SE(row0, col0) > 0){
              
              //WHAT IS THE PROBABILITY THAT A SPORE HITS ANY SUSCEPTIBLE HOST?
              PropS = double(S_UM(row0, col0) + S_OK(row0, col0) + S_LD(row0, col0) + S_AC(row0, col0) + 
                             S_AR(row0, col0) + S_AE(row0, col0) + S_PS(row0, col0) + S_SE(row0, col0) ) / N_LVE(row0, col0);
                
              double U = R::runif(0,1);			  
                
              //if U < PropS then one of the susceptible trees will get hit
              if (U < PropS){
                
                //WHICH SUSCEPTIBLE TREE GETS HIT? CALCULATE PROBABILITY WEIGHTS
                double Prob_UM = double (S_UM(row0, col0)) / N_LVE(row0, col0); 
                double Prob_LD = double (S_LD(row0, col0)) / N_LVE(row0, col0);
                double Prob_AC = double (S_AC(row0, col0)) / N_LVE(row0, col0);
                double Prob_AR = double (S_AR(row0, col0)) / N_LVE(row0, col0);
                double Prob_AE = double (S_AE(row0, col0)) / N_LVE(row0, col0);
                double Prob_PS = double (S_PS(row0, col0)) / N_LVE(row0, col0);
                double Prob_SE = double (S_SE(row0, col0)) / N_LVE(row0, col0);
                double Prob_OK = double (S_OK(row0, col0)) / N_LVE(row0, col0);
                
                //sample which of the three hosts will be hit given probability weights
                IntegerVector sv = sample(seq_len(8), 1, false, 
                                          NumericVector::create(Prob_UM, Prob_LD, Prob_AC, Prob_AR, 
                                                                Prob_AE, Prob_PS, Prob_SE, Prob_OK));
                
                //WHAT IS THE PROBABILITY THAT A SPORE TURNS INTO AN INFECTION, GIVEN THAT A SUSCEPTIBLE TREE HAS BEEN HIT?
                //This depends on both the local weather AND the host competency score (relative to UMCA and OAKS)
                int s = sv[0];
                if (s == 1){
                  double ProbINF = Prob_UM * W(row0, col0);
                  if (U < ProbINF){
                    I_UM(row0, col0) = I_UM(row0, col0) + 1; //update infected UMCA
                    S_UM(row0, col0) = S_UM(row0, col0) - 1; //update susceptible UMCA
                  }
                }else if (s == 2){
                  double ProbINF = Prob_LD * W(row0, col0);
                  if (U < ProbINF){
                    I_LD(row0, col0) = I_LD(row0, col0) + 1; //update infected LIDE
                    S_LD(row0, col0) = S_LD(row0, col0) - 1; //update susceptible LIDE                    
                  }
                }else if (s == 3){
                  double ProbINF = Prob_AC * W(row0, col0) * 0.1;  //hardcoded coeff to decrease prob.infection (transmission)
                  if (U < ProbINF){
                    I_AC(row0, col0) = I_AC(row0, col0) + 1; //update infected ACMA
                    S_AC(row0, col0) = S_AC(row0, col0) - 1; //update susceptible ACMA     
                  }
                }else if (s == 4){
                  double ProbINF = Prob_AR * W(row0, col0) * 0.1;  //hardcoded coeff to decrease prob.infection (transmission)
                  if (U < ProbINF){
                    I_AR(row0, col0) = I_AR(row0, col0) + 1; //update infected ARME
                    S_AR(row0, col0) = S_AR(row0, col0) - 1; //update susceptible ARME     				
                  }
                }else if (s == 5){
                  double ProbINF = Prob_AE * W(row0, col0) * 0.1;  //hardcoded coeff to decrease prob.infection (transmission)
                  if (U < ProbINF){
                    I_AE(row0, col0) = I_AE(row0, col0) + 1; //update infected AECA
                    S_AE(row0, col0) = S_AE(row0, col0) - 1; //update susceptible AECA  				
                  }
                }else if (s == 6){
                  double ProbINF = Prob_PS * W(row0, col0) * 0.1;  //hardcoded coeff to decrease prob.infection (transmission)
                  if (U < ProbINF){
                    I_PS(row0, col0) = I_PS(row0, col0) + 1; //update infected PSME
                    S_PS(row0, col0) = S_PS(row0, col0) - 1; //update susceptible PSME  				
                  }
                }else if (s == 7){
                  double ProbINF = Prob_SE * W(row0, col0) * 0.1;  //hardcoded coeff to decrease prob.infection (transmission)
                  if (U < ProbINF){
                    I_SE(row0, col0) = I_SE(row0, col0) + 1; //update infected SESE
                    S_SE(row0, col0) = S_SE(row0, col0) - 1; //update susceptible SESE 				
                  }
                }else{
                  double ProbINF = Prob_OK * W(row0, col0) * 0.75; //hardcoded coeff to decrease prob.infection (transmission)
                  if (U < ProbINF){
                    I_OK(row0, col0) = I_OK(row0, col0) + 1; //update infected OAKS
                    S_OK(row0, col0) = S_OK(row0, col0) - 1; //update susceptible OAKS  				
                  }
                }
              }//END IF CHECK vs UNIFORM NUMBER	                
              
            }//END IF NO SUSCEPTIBLE PRESENT IN CELL
          
          }else{  //IF DISTANCE IS OUTSIDE THE SAME CELL
              
            if(S_UM(row0, col0) > 0 || S_LD(row0, col0) > 0){  //IF SUSCEPTIBLE HOST IS AVAILABLE (UMCA OR LIDE)
              
              //WHAT IS THE PROBABILITY THAT A SPORE HITS ANY SUSCEPTIBLE HOST?
              PropS = double(S_UM(row0, col0) + S_LD(row0, col0)) / N_LVE(row0, col0);
                 
              double U = R::runif(0,1);		

              //if U < ProbS then one of the susceptible trees will get hit
              if (U < PropS){
                
                //WHICH SUSCEPTIBLE TREE GETS HIT? CALCULATE PROBABILITY WEIGHTS
                double Prob_UM = double (S_UM(row0, col0)) / N_LVE(row0, col0); 
                double Prob_LD = double (S_LD(row0, col0)) / N_LVE(row0, col0);

                //sample which of the three hosts will be hit given probability weights
                IntegerVector sv = sample(seq_len(2), 1, false, 
                                          NumericVector::create(Prob_UM, Prob_LD));
                
                //WHAT IS THE PROBABILITY THAT A SPORE TURNS INTO AN INFECTION, GIVEN THAT A SUSCEPTIBLE TREE HAS BEEN HIT?
                //This depends on both the local weather AND the host competency score (relative to UMCA and OAKS)
                int s = sv[0];
                if (s == 1){
                  double ProbINF = Prob_UM * W(row0, col0);
                  if (U < ProbINF){
                    I_UM(row0, col0) = I_UM(row0, col0) + 1; //update infected UMCA
                    S_UM(row0, col0) = S_UM(row0, col0) - 1; //update susceptible UMCA
                  }
                }else{
                  double ProbINF = Prob_LD * W(row0, col0);
                  if (U < ProbINF){
                    I_LD(row0, col0) = I_LD(row0, col0) + 1; //update infected LIDE
                    S_LD(row0, col0) = S_LD(row0, col0) - 1; //update susceptible LIDE                    
                  }                
                }  
              }//END IF U < PROB
            }//ENF IF S > 0  
          } //END IF DISTANCE CHECK
          
        }//END LOOP OVER ALL SPORES IN CURRENT CELL GRID     
      }//END IF SPORES IN CURRENT CELL
    }   
  }//END LOOP OVER ALL GRID CELLS
    
  //return List::create(Named("I")=I, Named("S")=S);
  return List::create(
    _["S_UM"] = S_UM, 
    _["I_UM"] = I_UM,
    _["S_OK"] = S_OK, 
    _["I_OK"] = I_OK,
    _["S_LD"] = S_LD, 
    _["I_LD"] = I_LD,
    _["S_AC"] = S_AC, 
    _["I_AC"] = I_AC,
    _["S_AR"] = S_AR, 
    _["I_AR"] = I_AR,
    _["S_AE"] = S_AE, 
    _["I_AE"] = I_AE,	
    _["S_PS"] = S_PS, 
    _["I_PS"] = I_PS,
    _["S_SE"] = S_SE, 
    _["I_SE"] = I_SE	
  );
    
} //END OF FUNCTION				  
  
// [[Rcpp::export]]
List SporeDispCppWind_MH(IntegerMatrix x, 
                         IntegerMatrix S_UM, IntegerMatrix S_LD, //SUSCEPTIBLE
                         IntegerMatrix S_AC, IntegerMatrix S_AR, 
                         IntegerMatrix S_AE, IntegerMatrix S_PS, 
                         IntegerMatrix S_SE, IntegerMatrix S_OK, 
                         IntegerMatrix I_UM, IntegerMatrix I_LD, //INFECTED
                         IntegerMatrix I_AC, IntegerMatrix I_AR, 
                         IntegerMatrix I_AE, IntegerMatrix I_PS, 
                         IntegerMatrix I_SE, IntegerMatrix I_OK, 
                         IntegerMatrix N_LVE,
                         NumericMatrix W,   //use different name than the functions in myfunctions_SOD.r
                         double rs, String rtype, double scale1, 
                         String wdir, int kappa,
                         double scale2=NA_REAL,  //default values
                         double gamma=NA_REAL)   //default values
{ 
    
  // internal variables //
  int nrow = x.nrow(); 
  int ncol = x.ncol();
  int row0;
  int col0;
    
  double dist;
  double theta;
  double PropS;
    
  //for Rcpp random numbers
  RNGScope scope;
    
  //Function rcauchy("rcauchy");
  Function rvm("rvm");
  Function sample("sample");
    
  //LOOP THROUGH EACH CELL of the input matrix 'x' (this should be the study area)
  for (int row = 0; row < nrow; row++) {
    for (int col = 0; col < ncol; col++){
        
      if(x(row,col) > 0){  //if spores in cell (row,col) > 0, disperse
          
        for(int sp = 1; (sp <= x(row,col)); sp++){
            
          //GENERATE DISTANCES:
          if (rtype == "Cauchy") 
            dist = abs(R::rcauchy(0, scale1));
          else if (rtype == "Cauchy Mixture"){
            if (gamma >= 1 || gamma <= 0) stop("The parameter gamma must range between (0-1)");
            NumericVector fv = sample(Range(1, 2), 1, false, NumericVector::create(gamma, 1-gamma));
            int f = fv[0];
            if(f == 1) 
              dist = abs(R::rcauchy(0, scale1));
            else 
              dist = abs(R::rcauchy(0, scale2));
          }
          else stop("The parameter rtype must be set to either 'Cauchy' or 'Cauchy Mixture'");
            
          //GENERATE ANGLES (using Von Mises distribution):
          if(kappa <= 0)  // kappa=concentration
            stop("kappa must be greater than zero!");
            
          //predominant wind dir
          if (wdir == "N") 
            theta = as<double>(rvm(1, 0 * (PI/180), kappa));  
          else if (wdir == "NE")
            theta = as<double>(rvm(1, 45 * (PI/180), kappa));  
          else if(wdir == "E")
            theta = as<double>(rvm(1, 90 * (PI/180), kappa));  
          else if(wdir == "SE")
            theta = as<double>(rvm(1, 135 * (PI/180), kappa));  
          else if(wdir == "S")
            theta = as<double>(rvm(1, 180 * (PI/180), kappa));  
          else if(wdir == "SW")
            theta = as<double>(rvm(1, 225 * (PI/180), kappa));  
          else if(wdir == "W")
            theta = as<double>(rvm(1, 270 * (PI/180), kappa));  
          else
            theta = as<double>(rvm(1, 315 * (PI/180), kappa));
            
          //calculate new row and col position for the dispersed spore unit (using dist and theta)
          row0 = row - round((dist * cos(theta)) / rs);
          col0 = col + round((dist * sin(theta)) / rs);
            
          if (row0 < 0 || row0 >= nrow) continue;     //outside the region
          if (col0 < 0 || col0 >= ncol) continue;     //outside the region
            
          //if distance is within same pixel challenge all SOD hosts, otherwise challenge UMCA only
          if (row0 == row && col0 == col){
              
            //if susceptible hosts are present in cell, calculate prob of infection
            if(S_UM(row0, col0) > 0 || S_OK(row0, col0) > 0 || S_LD(row0, col0) > 0 || S_AC(row0, col0) > 0 ||
               S_AR(row0, col0) > 0 || S_AE(row0, col0) > 0 || S_PS(row0, col0) > 0 || S_SE(row0, col0) > 0){
                
              //WHAT IS THE PROBABILITY THAT A SPORE HITS ANY SUSCEPTIBLE HOST?
              PropS = double(S_UM(row0, col0) + S_OK(row0, col0) + S_LD(row0, col0) + S_AC(row0, col0) + 
                             S_AR(row0, col0) + S_AE(row0, col0) + S_PS(row0, col0) + S_SE(row0, col0) ) / N_LVE(row0, col0);
                
              double U = R::runif(0,1);			  
                
              //if U < PropS then one of the susceptible trees will get hit
              if (U < PropS){
                  
                //WHICH SUSCEPTIBLE TREE GETS HIT? CALCULATE PROBABILITY WEIGHTS
                double Prob_UM = double (S_UM(row0, col0)) / N_LVE(row0, col0); 
                double Prob_LD = double (S_LD(row0, col0)) / N_LVE(row0, col0);
                double Prob_AC = double (S_AC(row0, col0)) / N_LVE(row0, col0);
                double Prob_AR = double (S_AR(row0, col0)) / N_LVE(row0, col0);
                double Prob_AE = double (S_AE(row0, col0)) / N_LVE(row0, col0);
                double Prob_PS = double (S_PS(row0, col0)) / N_LVE(row0, col0);
                double Prob_SE = double (S_SE(row0, col0)) / N_LVE(row0, col0);
                double Prob_OK = double (S_OK(row0, col0)) / N_LVE(row0, col0);
                  
                //sample which of the three hosts will be hit given probability weights
                IntegerVector sv = sample(seq_len(8), 1, false, 
                                          NumericVector::create(Prob_UM, Prob_LD, Prob_AC, Prob_AR, 
                                                                Prob_AE, Prob_PS, Prob_SE, Prob_OK));
                  
                int s = sv[0];
                if (s == 1){
                  double ProbINF = Prob_UM * W(row0, col0);
                  if (U < ProbINF){
                    I_UM(row0, col0) = I_UM(row0, col0) + 1; //update infected UMCA
                    S_UM(row0, col0) = S_UM(row0, col0) - 1; //update susceptible UMCA
                  }
                }else if (s == 2){
                  double ProbINF = Prob_LD * W(row0, col0);
                  if (U < ProbINF){
                    I_LD(row0, col0) = I_LD(row0, col0) + 1; //update infected LIDE
                    S_LD(row0, col0) = S_LD(row0, col0) - 1; //update susceptible LIDE                    
                  }
                }else if (s == 3){
                  double ProbINF = Prob_AC * W(row0, col0) * 0.1;  //hardcoded coeff to decrease prob.infection (transmission)
                  if (U < ProbINF){
                    I_AC(row0, col0) = I_AC(row0, col0) + 1; //update infected ACMA
                    S_AC(row0, col0) = S_AC(row0, col0) - 1; //update susceptible ACMA     
                  }
                }else if (s == 4){
                  double ProbINF = Prob_AR * W(row0, col0) * 0.1;  //hardcoded coeff to decrease prob.infection (transmission)
                  if (U < ProbINF){
                    I_AR(row0, col0) = I_AR(row0, col0) + 1; //update infected ARME
                    S_AR(row0, col0) = S_AR(row0, col0) - 1; //update susceptible ARME     				
                  }
                }else if (s == 5){
                  double ProbINF = Prob_AE * W(row0, col0) * 0.1;  //hardcoded coeff to decrease prob.infection (transmission)
                  if (U < ProbINF){
                    I_AE(row0, col0) = I_AE(row0, col0) + 1; //update infected AECA
                    S_AE(row0, col0) = S_AE(row0, col0) - 1; //update susceptible AECA  				
                  }
                }else if (s == 6){
                  double ProbINF = Prob_PS * W(row0, col0) * 0.1;  //hardcoded coeff to decrease prob.infection (transmission)
                  if (U < ProbINF){
                    I_PS(row0, col0) = I_PS(row0, col0) + 1; //update infected PSME
                    S_PS(row0, col0) = S_PS(row0, col0) - 1; //update susceptible PSME  				
                  }
                }else if (s == 7){
                  double ProbINF = Prob_SE * W(row0, col0) * 0.1;  //hardcoded coeff to decrease prob.infection (transmission)
                  if (U < ProbINF){
                    I_SE(row0, col0) = I_SE(row0, col0) + 1; //update infected SESE
                    S_SE(row0, col0) = S_SE(row0, col0) - 1; //update susceptible SESE 				
                  }
                }else{
                  double ProbINF = Prob_OK * W(row0, col0) * 0.75; //hardcoded coeff to decrease prob.infection (transmission)
                  if (U < ProbINF){
                    I_OK(row0, col0) = I_OK(row0, col0) + 1; //update infected OAKS
                    S_OK(row0, col0) = S_OK(row0, col0) - 1; //update susceptible OAKS  				
                  }
                }
              }//END IF CHECK vs UNIFORM NUMBER	                
              
            }//END IF NO SUSCEPTIBLE PRESENT IN CELL
            
          }else{  //IF DISTANCE IS OUTSIDE THE SAME CELL
            
            if(S_UM(row0, col0) > 0 || S_LD(row0, col0) > 0){  //IF SUSCEPTIBLE HOST IS AVAILABLE (UMCA OR LIDE)
              
              //WHAT IS THE PROBABILITY THAT A SPORE HITS ANY SUSCEPTIBLE HOST?
              PropS = double(S_UM(row0, col0) + S_LD(row0, col0)) / N_LVE(row0, col0);
              
              double U = R::runif(0,1);		
              
              //if U < ProbS then one of the susceptible trees will get hit
              if (U < PropS){
                
                //WHICH SUSCEPTIBLE TREE GETS HIT? CALCULATE PROBABILITY WEIGHTS
                double Prob_UM = double (S_UM(row0, col0)) / N_LVE(row0, col0); 
                double Prob_LD = double (S_LD(row0, col0)) / N_LVE(row0, col0);
                
                //sample which of the three hosts will be hit given probability weights
                IntegerVector sv = sample(seq_len(2), 1, false, 
                                          NumericVector::create(Prob_UM, Prob_LD));
                
                //WHAT IS THE PROBABILITY THAT A SPORE TURNS INTO AN INFECTION, GIVEN THAT A SUSCEPTIBLE TREE HAS BEEN HIT?
                //This depends on both the local weather AND the host competency score (relative to UMCA and OAKS)
                int s = sv[0];
                if (s == 1){
                  double ProbINF = Prob_UM * W(row0, col0);
                  if (U < ProbINF){
                    I_UM(row0, col0) = I_UM(row0, col0) + 1; //update infected UMCA
                    S_UM(row0, col0) = S_UM(row0, col0) - 1; //update susceptible UMCA
                  }
                }else{
                  double ProbINF = Prob_LD * W(row0, col0);
                  if (U < ProbINF){
                    I_LD(row0, col0) = I_LD(row0, col0) + 1; //update infected LIDE
                    S_LD(row0, col0) = S_LD(row0, col0) - 1; //update susceptible LIDE                    
                  }                
                }  
              }//END IF U < PROB
            }//ENF IF S > 0  
          } //END IF DISTANCE CHECK
          
        }//END LOOP OVER ALL SPORES IN CURRENT CELL GRID     
      }//END IF SPORES IN CURRENT CELL
    }   
  }//END LOOP OVER ALL GRID CELLS
  
  //return List::create(Named("I")=I, Named("S")=S);
  return List::create(
    _["S_UM"] = S_UM, 
    _["I_UM"] = I_UM,
    _["S_OK"] = S_OK, 
    _["I_OK"] = I_OK,
    _["S_LD"] = S_LD, 
    _["I_LD"] = I_LD,
    _["S_AC"] = S_AC, 
    _["I_AC"] = I_AC,
    _["S_AR"] = S_AR, 
    _["I_AR"] = I_AR,
    _["S_AE"] = S_AE, 
    _["I_AE"] = I_AE,	
    _["S_PS"] = S_PS, 
    _["I_PS"] = I_PS,
    _["S_SE"] = S_SE, 
    _["I_SE"] = I_SE	
  );
  
} //END OF FUNCTION		
    