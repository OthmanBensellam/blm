#include <iostream>
#include <string>
#include <cmath>
#include "jlparser/parser.hpp"
#include "MonteCarlo.hpp"
#include "Option.hpp"
#include "Geometric.hpp"
#include "Basket.hpp"
#include "Performance.hpp"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_random.h"
#include "PricingResults.hpp"
#include "tools.hpp"

using namespace std;

int main(int argc, char **argv){

  if (argc != 2) {
    std::cerr << "Wrong number of arguments. Exactly one arguments is required" << std::endl;
    std::exit(0);
  }

  PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
  pnl_rng_sseed(rng, time(NULL));

  double T, r, strike;
  PnlVect *spot, *sigma, *divid, *coeffs;

  string type;
  int size;
  int dates;
  int n_samples;
  int degree;
  double rho;
  double prixi;
  double prix;
  bool c = false;
  Option* opt;

  //Lécture du fichier d'entrée
  char *infile = argv[1];
  Param *P = new Parser(infile);

  P->extract("correlation", rho);
  P->extract("degree for polynomial regression", degree);
  P->extract("option type", type);
  P->extract("maturity", T);
  P->extract("model size", size);
  P->extract("dates", dates);
  P->extract("spot", spot, size);
  P->extract("volatility", sigma, size);
  P->extract("interest rate", r);
  if (P->extract("dividend rate", divid, size, true) == false)
  {
    divid = pnl_vect_create_from_zero(size);
  }
  P->extract("strike", strike);
  P->extract("MC iterations", n_samples);

  /*Determination du type de l'option 
  et construction de la classe correspondante*/
  if(type == "exchange" || type == "basket"){
    P->extract("payoff coefficients", coeffs, size);
    opt = new Basket(T, degree,dates, size, strike, coeffs);
    c = true;
  }
  else if(type == "bestof"){
    P->extract("payoff coefficients", coeffs, size);
    opt = new Performance(T, degree,dates, size, strike, coeffs);
    c = true;
  }
  else if(type=="geometric_put"){
    opt = new Geometric(T, degree,dates, size, strike);
  }
  else{
    std::cerr << "type unrecognized" << std::endl;
    std::exit(0);
  }

  //Création des classes BlackScholes et MonteCarlo
  BlackScholesModel* bs = new BlackScholesModel(size,rho,r, sigma, divid, spot);
  MonteCarlo* mc = new MonteCarlo(bs, opt, rng, n_samples);



  //************methode JL***************
  //Calcule du prix
  mc->price(prixi);
  std::cout << PricingResults(prixi) << std::endl;


  //************methode Article***************
  
  std::cout << "***********methode Article******" << std::endl;
  bs = new BlackScholesModel(size,rho,r, sigma, divid, spot);
  mc = new MonteCarlo(bs, opt, rng, n_samples);
  tools* tool = new tools(mc, bs, opt,n_samples);
  tool->optim();
  //std::cout << PricingResults(prix) << std::endl;



  pnl_vect_free(&spot);
  pnl_vect_free(&sigma);
  pnl_vect_free(&divid);
  if(c){
    pnl_vect_free(&coeffs);
  }
  delete P;
  delete opt;
  delete bs;
  delete mc;

  exit(0);
}