// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "../inst/include/sp.scRNAseq.h"
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


////////////////////////////////////////////////////////////////////////////////
//          ARMADILLO FUNCTIONS TO GENERATE SYNTHETIC MULTIPLETS              //
////////////////////////////////////////////////////////////////////////////////

//' sampleSinglets
//' 
//' This function takes a character vector of classes/cell with the same order
//' as the cells in the counts matrix. It returns one random index per unique
//' class and returns them as a integer vector.
//'
//' @param classes Character; a character vector of classes with length equal to
//' the number of cells for which counts exist.
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

arma::uvec sampleSinglets(
    CharacterVector classes
){
  CharacterVector uClasses = unique(classes).sort();
  arma::uvec idxToSubset(uClasses.size());
  
  //get indices for each cell type
  for (int y = 0; y < uClasses.size(); y++) {
    IntegerVector idxs;
    for (int j = 0; j < classes.size(); j++) {
      if(classes[j] == uClasses[y]) {idxs.push_back(j);}
    }
    
    //sample with length 1
    idxToSubset(y) = sample(idxs, 1)[0];
  }
  return idxToSubset;
}

//' subsetSinglets
//' 
//' This function takes a counts matrix and subsets the columns according to the
//' indices provided by the idxToSubset argument.
//'
//' @param singlets Matrix; a counts matrix with cells/samples as columns and 
//' genes as rows.
//' @param idxToSubset Integer; the indexes of cells/samples to be subset. 
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

arma::mat subsetSinglets(
    const arma::mat& singlets,
    const arma::uvec& idxToSubset
){
  return singlets.cols(idxToSubset);
}

//' calculateCost
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

double calculateCost(
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
