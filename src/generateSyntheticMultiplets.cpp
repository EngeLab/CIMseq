// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

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

//' normalizeFractionsArma
//' 
//' Takes a numeric vector and scales to [0, 1] by dividing with its sum.
//'
//' @param fractions Numeric; a numeric vector.
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

arma::vec normalizeFractionsArma(
    const arma::vec& fractions
){
  return fractions / accu(fractions);
}

//' adjustAccordingToFractionsArma
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

arma::mat adjustAccordingToFractionsArma(
    const arma::vec& fractions, 
    const arma::mat& singlets
){
  return singlets * arma::diagmat(fractions);
}

//' multipletSumsArma
//' 
//' This function takes a counts matrix and calculates the row sums. Output is
//' subsequently rounded for integration downstream.
//'
//' @param adjusted Matrix; a counts matrix with cells/samples as columns and 
//' genes as rows.
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

arma::mat multipletSumsArma(
    const arma::mat& adjusted
){
  return arma::round(arma::sum(adjusted, 1));
}

//' poissonSampleArma
//' 
//' This function takes the rowSums calculated by the 
//' \code{\link{multipletSumsArma}}.function and randomly samples 1 value for 
//' each input using the Poisson distribution and the input value as lambda.
//' @param rsRcpp Integer; vector of rounded row sums.
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

arma::vec poissonSampleArma(
    const arma::mat& rsRcpp
){
  arma::vec poissonSamp(rsRcpp.n_rows);
  for(int y = 0; y < rsRcpp.n_rows; y++) {
    poissonSamp[y] = Rcpp::rpois(1, rsRcpp(y))[0];
  }
  return poissonSamp;
}

//' vecToMatArma
//' 
//' This function takes a vector and reformats it to a matrix in a column-wise
//' fashion.
//'
//' @param vec Numeric; the vector to reformat.
//' @param nr integer; length 1 integer indicating the number of matrix rows.
//' @param nc integer; length 1 integer indicating the number of matrix columns.
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

arma::mat vecToMatArma(arma::vec vec, int nr, int nc){
  arma::mat X(vec.memptr(), nr, nc, false, false);
  return X;
}

//' cpmArma
//' 
//' Takes a numeric matrix of counts and calculates counts per million.
//'
//' @param counts Matrix; a counts matrix. Must be formated as numeric.
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

arma::mat cpmArma(
    const arma::mat& counts
){
  arma::rowvec colsums = sum(counts, 0);
  arma::mat cNorm = counts.each_row() / colsums;
  return (cNorm * 1000000) + 1;
}

////////////////////////////////////////////////////////////////////////////////
//                        FUNCTIONS TO CALCULATE COST ARMA                    //
////////////////////////////////////////////////////////////////////////////////

//' calculateCostDensityArma
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

arma::mat calculateCostDensityArma(
    arma::vec oneMultiplet,
    arma::mat syntheticMultiplets
){
  int nr = syntheticMultiplets.n_rows;
  int nc = syntheticMultiplets.n_cols;
  arma::mat densities(nr, nc);
  oneMultiplet = round(oneMultiplet);
  
  for(int i = 0; i < nr; i++) {
    for(int j = 0; j < nc; j++) {
      densities(i, j) = R::dpois(oneMultiplet(i), syntheticMultiplets(i, j), false);
    }
  }
  
  return densities;
}

//' calculateLogRowMeansArma
//' 
//' This function takes a matrix of poisson densities and calculates the row 
//' means and subsequently their log10 values.
//'
//' @param densities Numeric; a numeric vector of densities.
//' @author Jason T. Serviss
// [[Rcpp::export]]

arma::vec calculateLogRowMeansArma(
    const arma::mat& densities
){
  arma::vec means = mean(densities, 1);
  return log10(means);
}

//' fixNegInfArma
//' 
//' This function takes a numeric vector and replaces -Inf values with 
//' -323.0052. Note: since log10(10^-324) gives -Inf but log10(10^-323) 
//' gives -323.0052
//'
//' @param means Numeric; a numeric vector of log10 row means.
//' @author Jason T. Serviss
// [[Rcpp::export]]

arma::vec fixNegInfArma(
    arma::vec& means
){
  arma::vec noInfMeans(means.n_elem);
  means.elem(find_nonfinite(means)).fill(-323.0052);
  return means;
}

//' costNegSumArma
//' 
//' This function takes a numeric vector and calculates the negative sum.
//'
//' @param means Numeric; a numeric vector of log10 row means.
//' @author Jason T. Serviss
// [[Rcpp::export]]

double costNegSumArma(
    arma::vec means
){
  double cost = accu(means) *  -1;
  return cost;
}

//' calculateCostArma
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

double calculateCostArma(
    const arma::vec oneMultiplet,
    const arma::mat& syntheticMultiplets
){
  
  //calculate densities
  arma::mat ds = calculateCostDensityArma(oneMultiplet, syntheticMultiplets);
  
  //calculate log10 row means
  arma::vec means = calculateLogRowMeansArma(ds);
  
  //Replace -Inf with -323.0052
  arma::vec noInfMeans = fixNegInfArma(means);
  
  //calculate negative sum and return
  return costNegSumArma(noInfMeans);
}

//' preallocCost
//' 
//' Calculates cost with a preallocated matrix of subsetted singlets.
//'
//' @param oneMultiplet Numeric; a numeric vector of counts per million for one
//' multiplet.
//' @param singletSubset Matrix; Numeric matrix with the preallocated singlets. 
//' Each of the n synthetic multiplets should be stacked, i.e. rbind.
//' @param fractions Numeric; a numeric vector with length equal to 
//' ncol(singlets) indicating the fractions that each column should be 
//' multiplied with.
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

double preallocCost(
    const arma::vec& oneMultiplet,
    const arma::mat& singletSubset,
    const arma::vec& fractions
) {
  //check all 0 fractions
  if(accu(fractions) == 0) {
    return 999999999;
  }
  
  arma::vec f = normalizeFractionsArma(fractions);
  
  //adjust according to fractions
  arma::mat adjusted = adjustAccordingToFractionsArma(f, singletSubset);
  
  //rowSums
  arma::mat rs = multipletSumsArma(adjusted);
  
  //poisson sample
  arma::vec ps = poissonSampleArma(rs);
  
  //reformat into matrix
  arma::mat sm = vecToMatArma(ps, oneMultiplet.n_elem, f.n_elem);
  
  //calculate cpm
  arma::mat cpm = cpmArma(sm);
  
  //calculate cost
  double cost = calculateCostArma(oneMultiplet, cpm);
  return cost;
}
