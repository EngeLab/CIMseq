// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

////////////////////////////////////////////////////////////////////////////////
//          EIGEN FUNCTIONS TO GENERATE SYNTHETIC MULTIPLETS                  //
////////////////////////////////////////////////////////////////////////////////

//' normalizeFractionsEigen
//' 
//' Takes a numeric vector and scales to [0, 1] by dividing with its sum.
//'
//' @param fractions Numeric; a numeric vector.
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

Eigen::VectorXd normalizeFractionsEigen(
    const Eigen::VectorXd fractions
){
  Eigen::VectorXd normFractions = fractions.array() / fractions.sum();
  return normFractions;
}

//' sampleSingletsEigen
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

Eigen::VectorXi sampleSingletsEigen(
    CharacterVector classes
){
  CharacterVector uClasses = unique(classes).sort();
  Eigen::VectorXi idxToSubset(uClasses.size());
  
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

//' subsetSingletsEigen
//' 
//' This function takes a counts matrix and subsets the columns according to the
//' indices provided by the idxToSubset argument.
//'
//' @param singlets Matrix; a counts matrix with cells/samples as columns and 
//' genes as rows.
//' @param idxToSubset Integer; the indexes of cells/samples to be subset. 
//' Typically generared using the \code{\link{subsetSingletsEigen}} function.
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

Eigen::MatrixXd subsetSingletsEigen(
    Eigen::MatrixXd singlets,
    Eigen::VectorXi idxToSubset
){
  //Currently subsetting matrix by columns is not implemented in Eigen.
  //Functionality for this is in the devel version. As soon as it is implemented
  //change to that. Until them we convert to Armadillo, subset and reconvert to
  //Eigen
  
  //lots of conversions of idxToSubset to get right type for subsetting (painful)
  Eigen::VectorXd eigenD = idxToSubset.cast<double>();
  arma::vec armaVec = arma::vec(eigenD.data(), eigenD.size(), false, false);
  arma::uvec armaIdx = arma::conv_to<arma::uvec>::from(armaVec);
  
  //setup output
  arma::mat singletsArma = arma::mat(singlets.data(), singlets.rows(), singlets.cols(), false, false);
  
  //subset
  arma::mat subMat = singletsArma.cols(armaIdx);
  
  //convert back to Eigen
  Eigen::MatrixXd eigenMat = Eigen::Map<Eigen::MatrixXd>(subMat.memptr(), subMat.n_rows, subMat.n_cols);
  
  return eigenMat;
}

//' adjustAccordingToFractionsEigen
//' 
//' This function takes a counts matrix and subsets the columns according to the
//' indices provided by the idxToSubset argument.
//'
//' @param fractions Numeric; a numeric vector with length equal to 
//' \code{ncol(singlets)} indicating the fractions that each column should be 
//' multiplied with.
//' @param subMat Matrix; a counts matrix with cells/samples as columns and 
//' genes as rows. This matrix should have been previously subset with the 
//' \code{\link{sampleSingletsEigen}} function so that 
//' only one singlet per class/cell type is present.
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

Eigen::MatrixXd adjustAccordingToFractionsEigen(
    const Eigen::VectorXd fractions,
    const Eigen::MatrixXd subMat
){
  Eigen::MatrixXd adjusted = subMat * fractions.asDiagonal();
  return adjusted;
}

//' multipletSumsEigen
//' 
//' This function takes a counts matrix and calculates the row sums.
//'
//' @param singlets matrix; A counts matrix with cells/samples as columns and 
//' genes as rows. This matrix should have been previously subset with the 
//' \code{\link{subsetSingletsEigen}} function so that
//' only one singlet per class/cell type is present. Furthermore, it should have
//' already been asjusted by the fractions using the 
//' \code{\link{adjustAccordingToFractionsEigen}} function.
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

Eigen::VectorXd multipletSumsEigen(
    const Eigen::MatrixXd singlets
){
  Eigen::VectorXd sums = singlets.rowwise().sum().array().round();
  return sums;
}

//' poissonSampleEigen
//' 
//' This function takes the rowSums calculated by the 
//' \code{\link{multipletSumsEigen}} function and randomly samples 1 value for 
//' each input using the Poisson distribution and the input value as lambda.
//'
//' @param rsRcpp Integer; vector of rounded row sums.
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

Eigen::Map<Eigen::VectorXd> poissonSampleEigen(
    const Eigen::VectorXd rsRcpp
){
  NumericVector poissonSamp(rsRcpp.rows());
  for(int y = 0; y < rsRcpp.rows(); y++) {
    poissonSamp[y] = Rcpp::rpois(1, rsRcpp(y))[0];
  }
  Eigen::Map<Eigen::VectorXd> rcppPS = Rcpp::as<Eigen::Map<Eigen::VectorXd> >(poissonSamp);
  return rcppPS;
}

//' cpmEigen
//' 
//' Takes a numeric matrix of counts and calculates counts per million.
//'
//' @param counts Matrix; a counts matrix. Must be formated as numeric.
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

Eigen::MatrixXd cpmEigen(
    const Eigen::MatrixXd counts
){
  Eigen::VectorXd colsums = counts.colwise().sum();
  Eigen::MatrixXd dCounts = counts.transpose().array().colwise() / colsums.array();
  Eigen::MatrixXd cpm = (dCounts.transpose().array() * 1000000) + 1;
  return cpm;
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
//' \code{ncol(singlets)} indicating the fractions that each column should be 
//' multiplied with.
//' @param n Integer; length 1 integer indicating the number of synthetic 
//' multiplets to generate.
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

Eigen::MatrixXd generateSyntheticMultipletsEigen(
    Eigen::MatrixXd singlets,
    CharacterVector classes,
    Eigen::VectorXd fractions,
    int n
){
  
  //normalize fractions
  Eigen::VectorXd f = normalizeFractionsEigen(fractions);
  
  //setup output variables
  Eigen::MatrixXd syntheticMultiplets(singlets.rows(), n);
  
  //setup iterator
  int n0 = n - 1;
  
  //generate synthetic multiplets
  for (int o = 0; o <= n0; o++) {
    
    //get indices for each cell type and sample for subsequent singlet matrix subsetting
    Eigen::VectorXi idxToSubset = sampleSingletsEigen(classes);
    
    //subset singlets matrix
    Eigen::MatrixXd subMat = subsetSingletsEigen(singlets, idxToSubset);
    
    //adjust according to fractions
    Eigen::MatrixXd adjusted = adjustAccordingToFractionsEigen(f, subMat);
    
    //rowSums
    Eigen::VectorXd rs = multipletSumsEigen(adjusted);
    
    //poisson sample
    Eigen::Map<Eigen::VectorXd> ps = poissonSampleEigen(rs);
    
    //add to output
    syntheticMultiplets.col(o) = ps;
    
  }
  
  //calculate cpm and return
  return cpmEigen(syntheticMultiplets);
}

////////////////////////////////////////////////////////////////////////////////
//                        FUNCTIONS TO CALCULATE COST EIGEN                   //
////////////////////////////////////////////////////////////////////////////////

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

Eigen::MatrixXd calculateCostDensity(
    Eigen::Map<Eigen::VectorXd> oneMultiplet,
    const Eigen::MatrixXd syntheticMultiplets
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
  
  return densities;
}

//' calculateLogRowMeans
//' 
//' This function takes a matrix of poisson densities and calculates the row 
//' means and subsequently their log10 values.
//'
//' @param densities Numeric; a numeric vector of densities. Typically 
//' calculated with the \code{\link{calculateCostDensity}}.function.
//' @author Jason T. Serviss
// [[Rcpp::export]]

Eigen::VectorXd calculateLogRowMeans(
    const Eigen::MatrixXd densities
){
  Eigen::VectorXd means = densities.rowwise().mean().array().log10();
  return means;
}

//' fixNegInf
//' 
//' This function takes a numeric vector and replaces -Inf values with 
//' -323.0052. Note: since log10(10^-324) gives -Inf but log10(10^-323) 
//' gives -323.0052
//'
//' @param means Numeric; a numeric vector of log10 row means. Typically 
//' calculated with the \code{\link{calculateLogRowMeans}}.function.
//' @author Jason T. Serviss
// [[Rcpp::export]]

Eigen::VectorXd fixNegInf(
    const Eigen::VectorXd means
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
  return noInfMeans;
}

//' costNegSum
//' 
//' This function takes a numeric vector and calculates the negative sum.
//'
//' @param means Numeric; a numeric vector of log10 row means. Typically 
//' calculated with the \code{\link{calculateLogRowMeans}}.function with the 
//' -Inf values replaced by the \code{\link{fixNegInf}} function.
//' @author Jason T. Serviss
// [[Rcpp::export]]

double costNegSum(
    Eigen::VectorXd means
){
  double cost = means.sum() *  -1;
  return cost;
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

double calculateCostEigen(
    const Eigen::Map<Eigen::VectorXd> oneMultiplet,
    const Eigen::MatrixXd syntheticMultiplets
){
  
  //calculate densities
  Eigen::MatrixXd ds = calculateCostDensity(oneMultiplet, syntheticMultiplets);
  
  //calculate log10 row means
  Eigen::VectorXd means = calculateLogRowMeans(ds);
  
  //Replace -Inf with -323.0052
  Eigen::VectorXd noInfMeans = fixNegInf(means);
  
  //calculate negative sum and return
  return costNegSum(noInfMeans);
}

// WRAPPER FUNCTION FOR SYNTHETIC MULTIPLET GENERATION AND COST CALCULATION

//' calculateCostC
//' 
//' Wrapper for deconvolution cost function using the Eigen C++ library.
//'
//' @param oneMultiplet Numeric; vector with the gene expression for the 
//' multiplet currently being deconvoluted.
//' @param singlets Matrix; a counts matrix with cells/samples as columns and 
//' genes as rows.
//' @param classes Character; a character vector of classes with length equal to
//' the number of cells for which counts exist.
//' @param fractions Numeric; a numeric vector with length equal to 
//' \code{ncol(singlets)} indicating the fractions that each column should be 
//' multiplied with.
//' @param n Integer; length 1 integer indicating the number of synthetic 
//' multiplets to generate.
//' @author Jason T. Serviss
//' @export
// [[Rcpp::export]]

double calculateCostC(
    const Eigen::Map<Eigen::VectorXd> oneMultiplet,
    const Eigen::Map<Eigen::MatrixXd> singlets,
    const CharacterVector classes,
    const Eigen::VectorXd fractions,
    const int n
){
  
  //check all 0 fractions
  if(fractions.sum() == 0) {
    return 999999999;
  }
  
  //normalize fractions
  Eigen::VectorXd f = normalizeFractionsEigen(fractions);
  
  //setup output variables
  Eigen::MatrixXd syntheticMultiplets(singlets.rows(), n);
  
  //setup iterator
  int n0 = n - 1;
  
  //generate synthetic multiplets
  for (int o = 0; o <= n0; o++) {
    
    //get indices for each cell type and sample for subsequent singlet matrix subsetting
    Eigen::VectorXi idxToSubset = sampleSingletsEigen(classes);
    
    //subset singlets matrix
    Eigen::MatrixXd subMat = subsetSingletsEigen(singlets, idxToSubset);
    
    //adjust according to fractions
    Eigen::MatrixXd adjusted = adjustAccordingToFractionsEigen(f, subMat);
    
    //rowSums
    Eigen::VectorXd rs = multipletSumsEigen(adjusted);
    
    //poisson sample
    Eigen::Map<Eigen::VectorXd> ps = poissonSampleEigen(rs);
    
    //add to output
    syntheticMultiplets.col(o) = ps;
    
  }
  
  //calculate cpm and return
  Eigen::MatrixXd cpmSyntheticMultiplets = cpmEigen(syntheticMultiplets);
  
  //calculate cost
  double cost = calculateCostEigen(oneMultiplet, cpmSyntheticMultiplets);
  
  return cost;
}

