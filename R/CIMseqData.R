#' @include All-classes.R
NULL

#' CIMseqSinglets
#'
#' Subtitle
#'
#' Imports count, count.ercc, dimensionality reduction, and classification data
#' to a CIMseqSinglets object for sequenced singlets.
#'
#' @name CIMseqSinglets
#' @rdname CIMseqSinglets
#' @param counts matrix; Counts matrix with samples as columns and genes as rows.
#' @param counts.ercc matrix; A matrix containing ercc spike-in reads.
#' @param counts.cpm matrix; Normalized counts per million.
#' @param counts.log matrix; Log2 normalized counts per million.
#' @param dim.red matrix; Dimensionality reduced representation of the data.
#' @param classification character; Sample classes.
#' @param x CIMseqSinglets; CIMseqSinglets object.
#' @param object CIMseqSinglets; A CIMseqSinglets to show.
#' @param value Replacment value.
#' @param n Data to extract from CIMseqSinglets object.
#' @param .Object Internal object.
#' @param ... additional arguments to pass on
#' @return CIMseqSinglets object.
#' @author Jason T. Serviss
#' @examples
#'
#' #setup test data
#' s <- grepl("^s", colnames(testCounts))
#' ercc <- grepl("^ERCC\\-[0-9]*$", rownames(testCounts))
#' singlets <- testCounts[!ercc, s]
#' singletsERCC <- testCounts[ercc, s]
#' 
#' #dimensionality reduction
#' select <- selectTopMax(singlets, 100)
#' pdist <- pearsonsDist(singlets, select)
#' dim.red <- runTsne(pdist)
#' 
#' #classes
#' classes <- testMeta$cellTypes[match(colnames(singlets), testMeta$sample)]
#' #setup CIMseqSinglets object
#' cObj <- CIMseqSinglets(singlets, singletsERCC, dim.red, classes)
#'
NULL

#' @rdname CIMseqSinglets
#' @export

setGeneric("CIMseqSinglets", function(
  counts, ...
){
  standardGeneric("CIMseqSinglets")
})

#' @rdname CIMseqSinglets
#' @export

setMethod("CIMseqSinglets", "missing", function(
  ...
){
  new(
    "CIMseqSinglets",
    counts = matrix(nrow = 0, ncol = 0),
    counts.log = .norm.log.counts,
    counts.cpm = .norm.counts,
    counts.ercc = matrix(nrow=0, ncol = 0),
    dim.red = matrix(nrow = 0, ncol = 0),
    classification = character(),
    ...
  )
})

#' @rdname CIMseqSinglets
#' @export

setMethod("CIMseqSinglets", "matrix", function(
  counts, counts.ercc, dim.red, classification, ...
){
  .inputCheckCounts(counts, counts.ercc) #also test nrow(dim.red) and length(classification)
  new(
    "CIMseqSinglets",
    counts = counts,
    counts.log = .norm.log.counts,
    counts.cpm = .norm.counts,
    counts.ercc = counts.ercc,
    dim.red = dim.red,
    classification = classification,
    ...
  )
})

.inputCheckCounts <- function(counts, counts.ercc) {
    if((dim(counts)[2]) != (dim(counts.ercc)[2])) {
        message("ncol(counts) != ncol(counts.ercc).")
    }
    if(any(is.na(c(counts, counts.ercc)))) {
        message("is.na(c(counts, counts.ercc) returned TRUE")
    }
}

.norm.log.counts <- function(counts) {
  log2(.norm.counts(counts) + 1)
}

.norm.counts <- function(counts) {
  t(t(counts) / colSums(counts) * 10^6)
}

#' CIMseqMultiplets
#'
#' Subtitle
#'
#' Imports count, count.ercc, and feature data
#' to a CIMseqSinglets object for sequenced multiplets.
#'
#' @name CIMseqMultiplets
#' @rdname CIMseqMultiplets
#' @param counts matrix; Counts matrix with samples as columns and genes as rows.
#' @param counts.ercc matrix; A matrix containing ercc spike-in reads.
#' @param counts.cpm matrix; Normalized counts per million.
#' @param counts.log matrix; Log2 normalized counts per million.
#' @param features numeric; The indexes of the features/genes for use in 
#' deconvolution.
#' @param x CIMseqMultiplets; CIMseqMultiplets object.
#' @param object CIMseqMultiplets; A CIMseqMultiplets to show.
#' @param value Replacment value.
#' @param n Data to extract from CIMseqMultiplets object.
#' @param .Object Internal object.
#' @param ... additional arguments to pass on
#' @return CIMseqMultiplets object.
#' @author Jason T. Serviss
#' @examples
#'
#' s <- grepl("^s", colnames(testCounts))
#' ercc <- grepl("^ERCC\\-[0-9]*$", rownames(testCounts))
#' features <- selectTopMax(testCounts[!ercc, s], 2000)
#' cObj <- CIMseqMultiplets(testCounts[!ercc, !s], testCounts[ercc, !s], features)
#'
NULL

#' @rdname CIMseqMultiplets
#' @export

setGeneric("CIMseqMultiplets", function(
  counts, ...
){
  standardGeneric("CIMseqMultiplets")
})

#' @rdname CIMseqMultiplets
#' @export

setMethod("CIMseqMultiplets", "missing", function(
  ...
){
  new(
    "CIMseqMultiplets",
    counts = matrix(nrow = 0, ncol = 0),
    counts.log = .norm.log.counts,
    counts.cpm = .norm.counts,
    counts.ercc = matrix(nrow=0, ncol = 0),
    features = numeric(),
    ...
  )
})

#' @rdname CIMseqMultiplets
#' @export

setMethod("CIMseqMultiplets", "matrix", function(
  counts, counts.ercc, features, ...
){
  .inputCheckCounts(counts, counts.ercc)
  new(
    "CIMseqMultiplets",
    counts = counts,
    counts.log = .norm.log.counts,
    counts.cpm = .norm.counts,
    counts.ercc = counts.ercc,
    features = features,
    ...
  )
})

#' estimateCells
#'
#' Subtitle
#'
#' Uses ERCC data to calculate the fraction of ERCC reads in the samples. In
#' addition, this function utilizes ERCC data to estimate the cell number
#' in each sample.
#'
#' @name estimateCells
#' @rdname estimateCells
#' @param singlets A CIMseqSinglets object.
#' @param multiplets A CIMseqMultiplets object.
#' @param warning logical; Indicates if a warning should be issued when all ERCC
#'  counts for a sample are equal to 0. If this warning is issued it can be
#'  satisfactorily resolved by setting the ERCC reads for these samples to NA.
#' @param ... additional arguments to pass on
#' @return A data frame including the fraction of ercc reads and cell counts for
#'    each sample.
#' @author Jason T. Serviss
#' @examples
#'
#' output <- estimateCells(CIMseqSinglets_test, CIMseqMultiplets_test)
#'
NULL

#' @export

setGeneric("estimateCells", function(
  singlets,
  multiplets,
  ...
){
  standardGeneric("estimateCells")
})

#' @rdname estimateCells
#' @import tibble tibble
#' @importFrom stats median quantile
#' @importFrom rlang .data
#' @importFrom dplyr pull
#' @export

setMethod("estimateCells", "CIMseqSinglets", function(
  singlets, multiplets, warning = TRUE, ...
){
  
  counts <- cbind(
    getData(singlets, "counts"),
    getData(multiplets, "counts")
  )
  
  counts.ercc <- cbind(
    getData(singlets, "counts.ercc"),
    getData(multiplets, "counts.ercc")
  )
  
  #check if any samples have ERCC are all 0
  all0 <- apply(counts.ercc, 2, function(x) all(x == 0))
  if(any(all0, na.rm = TRUE) & warning) {
    zeroIDs <- colnames(counts.ercc)[which(all0)]
    if(length(zeroIDs) > 5) {
      zeroIDs <- paste0(paste(zeroIDs[1:5], collapse = ", "), ", ...")
    } else {
      zeroIDs <- paste(zeroIDs, collapse = ", ")
    }
    warning(paste0(
      "Results will not be accurate. ", 
      "These samples ERCC reads are all 0's: ", zeroIDs
    ))
  }
  
  d <- tibble(
    sampleName = c(
      colnames(getData(singlets, "counts")),
      colnames(getData(multiplets, "counts"))
    ),
    sampleType = c(
      rep("Singlet", ncol(getData(singlets, "counts"))),
      rep("Multiplet", ncol(getData(multiplets, "counts")))
    ),
    frac.ercc = colSums(counts.ercc) / (colSums(counts.ercc) + colSums(counts))
  )
  
  d$cellNumberMin <-
    d %>%
      filter(.data$sampleType == "Singlet") %>%
      pull(.data$frac.ercc) %>%
      quantile(., na.rm = TRUE) %>%
      `[` (2) %>%
      `/` (d$frac.ercc)
  
  d$cellNumberMedian <-
    d %>%
      filter(.data$sampleType == "Singlet") %>%
      pull(.data$frac.ercc) %>%
      median(., na.rm = TRUE) %>%
      `/` (d$frac.ercc)
  
  d$cellNumberMax <-
    d %>%
      filter(.data$sampleType == "Singlet") %>%
      pull(.data$frac.ercc) %>%
      quantile(., na.rm = TRUE) %>%
      `[` (4) %>%
      `/` (d$frac.ercc)
  
  return(d)
})
