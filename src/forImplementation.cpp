// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "../inst/include/sp.scRNAseq.h"
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' costFor
//'
//' Wrapper for deconvolution cost function using the Armadillo C++ library.
//'
//' @param oneMultiplet Integer; a integer vector of rounded counts per million
//' for one multiplet.
//' @param singletSubset Matrix; Numeric matrix with the preallocated singlets.
//' Each of the n synthetic multiplets should be stacked, i.e. rbind.
//' @param fractions Numeric; a numeric vector with length equal to
//' ncol(singlets) indicating the fractions that each column should be
//' multiplied with.
//' @param n Integer; length 1 integer indicating the number of synthetic
//' multiplets to generate.
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

double costFor(
  const arma::vec& oneMultiplet,
  const arma::mat& singletSubset,
  const arma::vec& fractions,
  int n
){
  //check all 0 fractions
  if(accu(fractions) == 0) {
    return 999999999;
  }
  
  //normalize fractions
  arma::vec normFractions = fractions / accu(fractions);
  
  //optimize column-major
  arma::mat tSinglets = singletSubset.t();
  
  //loop over genes
  double cost = 0;
  int ci = 0;
  // Loop over genes
  for (int i = 0; i != (tSinglets.n_cols / n); i++){
    // Loop over synthetic multiplets
    double pSums = 0;
    for(int k = 0; k != n; k++) {
      double adjustedSums = 0;
      //Loop over fractions
      for (int j = 0; j != normFractions.n_elem; j++) {
        adjustedSums += normFractions(j) * tSinglets(j, ci);
      }
      //double rp = Rcpp::rpois(1, std::round(adjustedSums))[0];
      pSums += R::dpois(oneMultiplet(i), std::round(adjustedSums), false);
      ++ci;
    }
    
    pSums /= double(n);
    double lpSums = std::log10(pSums);
    if(lpSums == -std::numeric_limits<double>::infinity()) {
      lpSums = -323.0052;
    }
    cost += lpSums;
  }
  return -cost;
}
