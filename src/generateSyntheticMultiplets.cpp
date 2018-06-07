// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// FUNCTIONS TO GENERATE SYNTHETIC MULTIPLETS

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

IntegerVector sampleSinglets(
    CharacterVector classes
){
  CharacterVector uClasses = unique(classes).sort();
  IntegerVector idxToSubset(uClasses.size());
  //get indices for each cell type
  for (int y = 0; y < uClasses.size(); y++) {
    IntegerVector idxs;
    for (int j = 0; j < classes.size(); j++) {
      if(classes[j] == uClasses[y]) {idxs.push_back(j);}
    }
    //sample with length 1
    idxToSubset[y] = sample(idxs, 1)[0];
  }
  return idxToSubset;
}

//' subsetSingletsEigen
//' 
//' This function takes a counts matrix and subsets the columns according to the
//' indices provided by the idxToSubset argument.
//'
//' @param singlets Matrix; a counts matrix with cells/samples as columns and 
//' genes as rows.
//' @param idxToSubset Integer; the indexes of cells/samples to be subset. 
//' Typically generared using the \link{\code{subsetSinglets}} function.
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

Eigen::Map<Eigen::MatrixXd> subsetSingletsEigen(
    const arma::mat& singlets,
    Rcpp::NumericVector idxToSubset
){
  arma::uvec v = Rcpp::as<arma::uvec>(idxToSubset);
  arma::mat subMat = singlets.cols(v);
  Eigen::Map<Eigen::MatrixXd> eigenMat = Eigen::Map<Eigen::MatrixXd>(
    subMat.memptr(), subMat.n_rows, subMat.n_cols
  );
  return eigenMat;
}

//' adjustAccordingToFractionsEigen
//' 
//' This function takes a counts matrix and subsets the columns according to the
//' indices provided by the idxToSubset argument.
//'
//' @param fractions Numeric; a numeric vector with length equal to 
//' ncol(singlets) indicating the fractions that each column should be 
//' multiplied with.
//' @param singlets Matrix; a counts matrix with cells/samples as columns and 
//' genes as rows. This matrix should have been previously subset with the 
//' \link{\code{sampleSinglets}} and \code{subsetSingletsEigen} functions so that 
//' only one singlet per class/cell type is present.
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

SEXP adjustAccordingToFractionsEigen(
  Eigen::Map<Eigen::VectorXd> fractions,
  const Eigen::Map<Eigen::MatrixXd> singlets
){
  return Rcpp::wrap(singlets * fractions.asDiagonal());
}

//' multipletSumsEigen
//' 
//' This function takes a counts matrix and calculates the row sums.
//'
//' @param adjusted matrix; A counts matrix with cells/samples as columns and 
//' genes as rows. This matrix should have been previously subset with the 
//' \link{\code{sampleSinglets}} and \code{subsetSingletsEigen} functions so that
//' only one singlet per class/cell type is present. Furthermore, it should have
//' already been asjusted by the fractions using the 
//' \link{\code{adjustAccordingToFractionsEigen}} function.
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

SEXP multipletSumsEigen(
  Eigen::MatrixXd singlets
){
  return Rcpp::wrap(singlets.rowwise().sum().array().round());
}

//' poissonSample
//' 
//' This function takes the rowSums calculated by the 
//' \link{\code{multipletSumsArma}}.or \link{\code{multipletSumsEigen}} 
//' functions and randomly samples 1 value for each input using the Poisson 
//' distribution and the input value as lambda.
//'
//' @param rsRcpp Integer; vector of rounded row sums.
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

NumericVector poissonSample(
    arma::mat rsRcpp
){
  NumericVector poissonSamp(rsRcpp.n_rows);
  for(int y = 0; y < rsRcpp.n_rows; y++) {
    poissonSamp[y] = Rcpp::rpois(1, rsRcpp(y))[0];
  }
  return poissonSamp;
}

//the funciton below is faster but the input will need to be converted
std::vector<double> poissonSample5(
    std::vector<double> rsRcpp
){
  std::vector<double> poissonSamp;
  poissonSamp.reserve(rsRcpp.size());
  for(int y = 0; y < rsRcpp.size(); y++) {
    poissonSamp[y] = Rcpp::rpois(1, rsRcpp[y])[0];
  }
  return poissonSamp;
}

//' cpmC
//' 
//' Takes a numeric matrix of counts and calculates counts per million.
//'
//' @param counts Matrix; a counts matrix. Must be formated as numeric.
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

SEXP cpmC(
    const Eigen::Map<Eigen::MatrixXd> counts
){
  Eigen::VectorXd colsums = counts.colwise().sum();
  Eigen::MatrixXd dCounts = counts.transpose().array().colwise() / colsums.array();
  return wrap((dCounts.transpose().array() * 1000000) + 1);
}

//' generateSyntheticMultipletsEigen
//' 
//' Wrapper for deconvolution cost function using the Eigen C++ library.
//'
//' @param singlets Matrix; a counts matrix with cells/samples as columns and 
//' genes as rows.
//' @param classes Character; a character vector of classes with length equal to
//' the number of cells for which counts exist.
//' @param fractions Numeric; a numeric vector with length equal to 
//' ncol(singlets) indicating the fractions that each column should be 
//' multiplied with.
//' @param n Integer; length 1 integer indicating the number of synthetic 
//' multiplets to generate.
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

NumericMatrix generateSyntheticMultipletsEigen(
  const arma::mat& singlets,
  CharacterVector classes,
  Eigen::Map<Eigen::VectorXd> fractions,
  int n
){
  //arma::arma_rng::set_seed_random();
  NumericMatrix syntheticMultiplets(singlets.n_rows, n);
  int n0 = n - 1;
  
  for (int o = 0; o <= n0; o++) {
    
    //get indices for each cell type and sample for subsequent singlet matrix subsetting
    CharacterVector uClasses = unique(classes).sort();
    IntegerVector idxToSubset;
    //get indices for each cell type
    for (int y = 0; y < uClasses.size(); y++) {
      IntegerVector idxs;
      for (int j = 0; j < classes.size(); j++) {
        if(classes[j] == uClasses[y]) {idxs.push_back(j);}
      }
      //sample with length 1
      //int samp = sample(idxs, 1)[0];
      idxToSubset.push_back(sample(idxs, 1)[0]);
    }
    
    //subset singlets matrix
    arma::uvec v = Rcpp::as<arma::uvec>(idxToSubset);
    arma::mat subMat = singlets.cols(v);
    Eigen::MatrixXd eigenMat = Eigen::Map<Eigen::MatrixXd>(
      subMat.memptr(), subMat.n_rows, subMat.n_cols
    );
    
    //adjust according to fractions
    Eigen::MatrixXd adjusted = eigenMat * fractions.asDiagonal();
    
    //rowSums
    Eigen::VectorXd rs = adjusted.rowwise().sum().array().round();
    
    //poisson sample
    NumericVector ps(rs.rows());
    for(int y = 0; y < rs.rows(); y++) {
      ps[y] = Rcpp::rpois(1, rs(y))[0];
    }
    
    syntheticMultiplets(_, o) = ps;
  }
  
  //calculate cpm
  Eigen::Map<Eigen::MatrixXd> sm_eigen = as<Eigen::Map<Eigen::MatrixXd> >(syntheticMultiplets);
  Eigen::VectorXd colsums = sm_eigen.colwise().sum();
  Eigen::MatrixXd dCounts = sm_eigen.transpose().array().colwise() / colsums.array();
  return wrap((dCounts.transpose().array() * 1000000) + 1);
}

// FUNCTIONS TO CALCULATE COST

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

SEXP calculateCostDensity(
    Eigen::Map<Eigen::VectorXi> oneMultiplet,
    const Eigen::Map<Eigen::MatrixXd> syntheticMultiplets
){
  int nr = syntheticMultiplets.rows();
  int nc = syntheticMultiplets.cols();
  Eigen::MatrixXd densities(nr, nc);
  oneMultiplet = oneMultiplet.array().round();
  
  for(int i = 0; i < nr; i++) {
    for(int j = 0; j < nc; j++) {
      densities(i, j) = R::dpois(oneMultiplet(i), syntheticMultiplets(i, j), false);
    }
  }
  
  return Rcpp::wrap(densities);
}

//' calculateLogRowMeans
//' 
//' This function takes a matrix of poisson densities and calculates the row 
//' means and subsequently their log10 values.
//'
//' @param densities Numeric; a numeric vector of densities. Typically 
//' calculated with the \link{\code{calculateCostDensity}}.function.
//' @author Jason T. Serviss
// [[Rcpp::export]]

SEXP calculateLogRowMeans(
    const Eigen::Map<Eigen::MatrixXd> densities
){
  Eigen::VectorXd means = densities.rowwise().mean().array().log10();
  return wrap(means);
}

//' fixNegInf
//' 
//' This function takes a numeric vector and replaces -Inf values with 
//' -323.0052. Note: since log10(10^-324) gives -Inf but log10(10^-323) 
//' gives -323.0052
//'
//' @param means Numeric; a numeric vector of log10 row means. Typically 
//' calculated with the \link{\code{calculateLogRowMeans}}.function.
//' @author Jason T. Serviss
// [[Rcpp::export]]

SEXP fixNegInf(
    Eigen::VectorXd means
){
  int ms = means.size();
  Eigen::VectorXd noInfMeans(ms);
  Eigen::Array<bool, 1, Eigen::Dynamic> infBool;
  infBool = means.array().isInf();
  
  for(int k = 0; k < ms; k++) {
    if(infBool(k)) {
      noInfMeans(k) = -323.0052;
    } else {
      noInfMeans(k) = means(k);
    }
  }
  return wrap(noInfMeans);
}

//' costNegSum
//' 
//' This function takes a numeric vector and calculates the negative sum.
//'
//' @param means Numeric; a numeric vector of log10 row means. Typically 
//' calculated with the \link{\code{calculateLogRowMeans}}.function with the 
//' -Inf values replaced by the \link{\code{fixNegInf}} function.
//' @author Jason T. Serviss
// [[Rcpp::export]]

SEXP costNegSum(
    Eigen::VectorXd means
){
  double cost = means.sum() *  -1;
  return wrap(cost);
}

//' calculateCostEigen
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

SEXP calculateCostEigen(
    Eigen::Map<Eigen::VectorXi> oneMultiplet,
    const Eigen::Map<Eigen::MatrixXd> syntheticMultiplets
){
  
  //calculate densities
  int nr = syntheticMultiplets.rows();
  int nc = syntheticMultiplets.cols();
  Eigen::MatrixXd ds(nr, nc);
  oneMultiplet = oneMultiplet.array().round();
  
  for(int i = 0; i < nr; i++) {
    for(int j = 0; j < nc; j++) {
      ds(i, j) = R::dpois(oneMultiplet[i], syntheticMultiplets(i, j), false);
    }
  }
  
  //calculate log10 row means
  Eigen::VectorXd means = ds.rowwise().mean().array().log10();
  
  //Replace -Inf with -323.0052
  Eigen::VectorXd noInfMeans(nr);
  Eigen::Array<bool, 1, Eigen::Dynamic> infBool;
  infBool = means.array().isInf();
  
  for(int k = 0; k < nr; k++) {
    if(infBool(k)) {
      noInfMeans(k) = -323.0052;
    } else {
      noInfMeans(k) = means(k);
    }
  }
  
  //calculate negative sum
  double cost = noInfMeans.sum() * -1;
  
  return wrap(cost);
}

// WRAPPER FUNCTION FOR SYNTHETIC MULTIPLET GENERATION AND COST CALCULATION

//' calculateCostC
//' 
//' Wrapper for deconvolution cost function using the Eigen C++ library.
//'
//' @param singlets Matrix; a counts matrix with cells/samples as columns and 
//' genes as rows.
//' @param classes Character; a character vector of classes with length equal to
//' the number of cells for which counts exist.
//' @param fractions Numeric; a numeric vector with length equal to 
//' ncol(singlets) indicating the fractions that each column should be 
//' multiplied with.
//' @param n Integer; length 1 integer indicating the number of synthetic 
//' multiplets to generate.
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

double calculateCostC(
    Eigen::Map<Eigen::VectorXd> oneMultiplet,
    const arma::mat& singlets,
    CharacterVector classes,
    Eigen::Map<Eigen::VectorXd> fractions,
    int n
){
  //arma::arma_rng::set_seed_random();
  NumericMatrix syntheticMultiplets(singlets.n_rows, n);
  int n0 = n - 1;
  
  for (int o = 0; o <= n0; o++) {
    
    //get indices for each cell type and sample for subsequent singlet matrix subsetting
    CharacterVector uClasses = unique(classes).sort();
    IntegerVector idxToSubset;
    //get indices for each cell type
    for (int y = 0; y < uClasses.size(); y++) {
      IntegerVector idxs;
      for (int j = 0; j < classes.size(); j++) {
        if(classes[j] == uClasses[y]) {idxs.push_back(j);}
      }
      //sample with length 1
      //int samp = sample(idxs, 1)[0];
      idxToSubset.push_back(sample(idxs, 1)[0]);
    }
    
    //subset singlets matrix
    arma::uvec v = Rcpp::as<arma::uvec>(idxToSubset);
    arma::mat subMat = singlets.cols(v);
    Eigen::MatrixXd eigenMat = Eigen::Map<Eigen::MatrixXd>(
      subMat.memptr(), subMat.n_rows, subMat.n_cols
    );
    
    //adjust according to fractions
    Eigen::MatrixXd adjusted = eigenMat * fractions.asDiagonal();
    
    //rowSums
    Eigen::VectorXd rs = adjusted.rowwise().sum().array().round();
    
    //poisson sample
    NumericVector ps(rs.rows());
    for(int y = 0; y < rs.rows(); y++) {
      ps[y] = Rcpp::rpois(1, rs(y))[0];
    }
    
    syntheticMultiplets(_, o) = ps;
  }
  
  //calculate cpm
  Eigen::Map<Eigen::MatrixXd> sm_eigen = as<Eigen::Map<Eigen::MatrixXd> >(syntheticMultiplets);
  Eigen::VectorXd colsums = sm_eigen.colwise().sum();
  Eigen::MatrixXd dCounts = sm_eigen.transpose().array().colwise() / colsums.array();
  Eigen::MatrixXd normSyntheticMultipliet = (dCounts.transpose().array() * 1000000) + 1;
  
  //calculate cost
  //calculate densities
  int nr = normSyntheticMultipliet.rows();
  int nc = normSyntheticMultipliet.cols();
  Eigen::MatrixXd ds(nr, nc);
  oneMultiplet = oneMultiplet.array().round();
  
  for(int i = 0; i < nr; i++) {
    for(int j = 0; j < nc; j++) {
      ds(i, j) = R::dpois(oneMultiplet[i], normSyntheticMultipliet(i, j), false);
    }
  }
  
  //calculate log10 row means
  Eigen::VectorXd means = ds.rowwise().mean().array().log10();
  
  //Replace -Inf with -323.0052
  Eigen::VectorXd noInfMeans(nr);
  Eigen::Array<bool, 1, Eigen::Dynamic> infBool;
  infBool = means.array().isInf();
  
  for(int k = 0; k < nr; k++) {
    if(infBool(k)) {
      noInfMeans(k) = -323.0052;
    } else {
      noInfMeans(k) = means(k);
    }
  }
  
  //calculate negative sum
  double cost = noInfMeans.sum() * -1;
  return cost;
}
