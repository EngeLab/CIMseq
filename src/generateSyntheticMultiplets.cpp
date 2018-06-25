// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "../inst/include/sp.scRNAseq.h"
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


////////////////////////////////////////////////////////////////////////////////
//          ARMADILLO FUNCTIONS TO GENERATE SYNTHETIC MULTIPLETS              //
////////////////////////////////////////////////////////////////////////////////

//' sampleSingletsArma
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

arma::uvec sampleSingletsArma(
    CharacterVector classes
){
  CharacterVector uClasses = unique(classes);
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

//' subsetSingletsArma
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

arma::mat subsetSingletsArma(
    const arma::mat& singlets,
    const arma::uvec& idxToSubset
){
  return singlets.cols(idxToSubset);
}
