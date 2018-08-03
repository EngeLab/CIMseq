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

////////////////////////////////////////////////////////////////////////////////
//                        FUNCTIONS FOR TESTING ETC.                          //
////////////////////////////////////////////////////////////////////////////////

//' normalizeFractions
//'
//' Takes a numeric vector and scales to [0, 1] by dividing with its sum.
//'
//' @param fractions Numeric; a numeric vector.
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

arma::vec normalizeFractions(
    const arma::vec& fractions
){
  return fractions / accu(fractions);
}

//' adjustAccordingToFractions
//'
//' This function takes a counts matrix and subsets the columns according to the
//' indices provided by the idxToSubset argument.
//'
//' @param fractions Numeric; a numeric vector with length equal to
//' ncol(singlets) indicating the fractions that each column should be
//' multiplied with.
//' @param singlets Matrix; a counts matrix with cells/samples as columns and
//' genes as rows.
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

arma::mat adjustAccordingToFractions(
    const arma::vec& fractions,
    const arma::mat& singlets
){
  return singlets * arma::diagmat(fractions);
}

//' multipletSums
//'
//' This function takes a counts matrix and calculates the row sums. Output is
//' subsequently rounded for integration downstream.
//'
//' @param adjusted Matrix; a counts matrix with cells/samples as columns and
//' genes as rows.
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

arma::mat multipletSums(
    const arma::mat& adjusted
){
  return arma::round(arma::sum(adjusted, 1));
}

//' vecToMat
//'
//' This function takes a vector and reformats it to a matrix in a column-wise
//' fashion. Since Armadillo "builds" the matrix column-wise, some code tricks
//' are employed to automagically deliever the expected result.
//'
//' @param vec Numeric; the vector to reformat.
//' @param nr integer; length 1 integer indicating the number of matrix rows.
//' @param nc integer; length 1 integer indicating the number of matrix columns.
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

arma::mat vecToMat(arma::vec vec, int nr, int nc){
  arma::mat X(vec.memptr(), nc, nr, false, false);
  return X.t();
}

//' calculateCostDensity
//'
//' This function takes a vector of gene counts per million for one multiplet
//' that is being deconvoluted and a matrix of synthetic multiplets. It then
//' calculates the poisson density for each gene count and each corresponding
//' value in the matrix row.
//'
//' @param oneMultiplet Numeric; a numeric vector of counts per million for one
//' multiplet.
//' @param syntheticMultiplets Matrix; a numeric matrix of synthetic multiplets
//' with samples as columns and genes as rows.
//' @author Jason T. Serviss
// [[Rcpp::export]]

arma::mat calculateCostDensity(
    arma::vec oneMultiplet,
    arma::mat syntheticMultiplets
){
  int nr = syntheticMultiplets.n_rows;
  int nc = syntheticMultiplets.n_cols;
  arma::mat densities(nr, nc);
  //oneMultiplet = round(oneMultiplet); should be done upstream
  
  for(int i = 0; i < nr; i++) {
    for(int j = 0; j < nc; j++) {
      densities(i, j) = R::dpois(oneMultiplet(i), syntheticMultiplets(i, j), false);
    }
  }
  
  return densities;
}

//' calculateLogRowMeans
//'
//' This function takes a matrix of poisson densities and calculates the row
//' means and subsequently their log10 values.
//'
//' @param densities Numeric; a numeric vector of densities.
//' @author Jason T. Serviss
// [[Rcpp::export]]

arma::vec calculateLogRowMeans(
    const arma::mat& densities
){
  arma::vec means = mean(densities, 1);
  return log10(means);
}

//' fixNegInf
//'
//' This function takes a numeric vector and replaces -Inf values with
//' -323.0052. Note: since log10(10^-324) gives -Inf but log10(10^-323)
//' gives -323.0052
//'
//' @param means Numeric; a numeric vector of log10 row means.
//' @author Jason T. Serviss
// [[Rcpp::export]]

arma::vec fixNegInf(
    arma::vec& means
){
  arma::vec noInfMeans(means.n_elem);
  means.elem(find_nonfinite(means)).fill(-323.0052);
  return means;
}

//' costNegSum
//'
//' This function takes a numeric vector and calculates the negative sum.
//'
//' @param means Numeric; a numeric vector of log10 row means.
//' @author Jason T. Serviss
// [[Rcpp::export]]

double costNegSum(
    arma::vec means
){
  double cost = accu(means) *  -1;
  return cost;
}

//' costCalc
//'
//' This function takes a vector of gene counts per million for one multiplet
//' that is being deconvoluted and a matrix of synthetic multiplets and
//' calculates the cost which is returned to the optimization algorithm during
//' deconvolution.
//'
//' @param oneMultiplet Numeric; a numeric vector of counts per million for one
//' multiplet.
//' @param syntheticMultiplets Matrix; a numeric matrix of synthetic multiplets
//' with samples as columns and genes as rows.
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

double costCalc(
    const arma::vec oneMultiplet,
    const arma::mat& syntheticMultiplets
){
  
  //calculate densities
  arma::mat ds = calculateCostDensity(oneMultiplet, syntheticMultiplets);
  
  //calculate log10 row means
  arma::vec means = calculateLogRowMeans(ds);
  
  //Replace -Inf with -323.0052
  arma::vec noInfMeans = fixNegInf(means);
  
  //calculate negative sum and return
  return costNegSum(noInfMeans);
}

//' cost
//'
//' Wrapper for deconvolution cost function using the Armadillo C++ library.
//'
//' @param oneMultiplet Integer; a integer vector of rounded counts per million
//' for one multiplet.
//' @param singletSubset Matrix; a counts matrix or pre-subsetted singlets. 
//' Typically from the .subsetSinglets function.
//' @param fractions Numeric; a numeric vector with length equal to
//' ncol(singlets) indicating the fractions that each column should be
//' multiplied with.
//' @param n Integer; length 1 integer indicating the number of synthetic
//' multiplets to generate.
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

double cost(
    const arma::vec& oneMultiplet,
    const arma::mat& singletSubset,
    const arma::vec& fractions,
    const int n
){
  
  //check all 0 fractions
  if(accu(fractions) == 0) {
    return 999999999;
  }
  
  //normalize fractions
  arma::vec f = normalizeFractions(fractions); //note: returns a matrix(?)
  
  //adjust according to fractions
  arma::mat adjusted = adjustAccordingToFractions(f, singletSubset);
  
  //rowSums
  arma::mat rs = multipletSums(adjusted);
  
  //reformat to wide
  arma::mat sm = vecToMat(rs, oneMultiplet.n_elem, n);
  
  //calculate cost
  double cost = costCalc(oneMultiplet, sm);
  
  return cost;
}

