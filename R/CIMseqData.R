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
    counts.ercc = matrix(nrow = 0, ncol = 0),
    classification = character(),
    norm.to=1E6,
    dim.red = matrix(nrow = 0, ncol = 0)
  )
})

#' @rdname CIMseqSinglets
#' @export

setMethod("CIMseqSinglets", "matrix", function(counts, counts.ercc=matrix(nrow=0, ncol=0),
                                               classification, norm.to = 1E6, dim.red=matrix(nrow=0, ncol=0))
{
    .inputCheckSinglets(counts, counts.ercc, dim.red, classification)
    new(
        "CIMseqSinglets",
        counts = counts,
        counts.log = .norm.log.counts,
        counts.cpm = .norm.counts,
        counts.ercc = counts.ercc,
        dim.red = dim.red,
        classification = classification,
        norm.to = norm.to
        )
})

.inputCheckSinglets <- function(counts, counts.ercc, dim.red, classification) {
  if(!length(counts.ercc)) {
    message("No ERCC controls given, will use counts for size estimation")
  }
  if(length(counts.ercc) & (dim(counts)[2]) != (dim(counts.ercc)[2])) {
    message("ncol(counts) != ncol(counts.ercc).")
  }
  if(any(is.na(c(counts, counts.ercc)))) {
    message("is.na(c(counts, counts.ercc) returned TRUE")
  }
  if(ncol(counts) != length(classification)) {
    message("length(classification) != ncol(counts)")
  }
  if(nrow(dim.red) != 0 & nrow(dim.red) != ncol(counts)) {
      message("nrow(dim.red) != ncol(counts)")
  }
}

.norm.log.counts <- function(counts, norm.to) {
  log2(.norm.counts(counts, norm.to) + 1)
}

.norm.counts <- function(counts, norm.to) {
  t(t(counts) / colSums(counts) * norm.to)
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
  counts, counts.ercc=matrix(nrow=0, ncol=0), features, norm.to = 1E6, ...
){
  featInd <- features
  if(is.character(features)) {
      featInd <- match(features, rownames(counts))
  }
  .inputCheckMultiplets(counts, counts.ercc, featInd)
  new(
    "CIMseqMultiplets",
    counts = counts,
    counts.log = .norm.log.counts,
    counts.cpm = .norm.counts,
    counts.ercc = counts.ercc,
    features = featInd,
    norm.to = norm.to,
    ...
  )
})

.inputCheckMultiplets <- function(counts, counts.ercc, features) {
  if(!length(counts.ercc)) {
    message("No ERCC controls given, will use counts for size estimation")
  }
  if(length(counts.ercc) & (dim(counts)[2]) != (dim(counts.ercc)[2])) {
    message("ncol(counts) != ncol(counts.ercc).")
  }
  if(any(is.na(c(counts, counts.ercc)))) {
    message("is.na(c(counts, counts.ercc) returned TRUE")
  }
  if(any(is.na(features))) {
      message("features contains non-existing genes or NA")
  }
}

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
#' @importFrom stats median
#' @importFrom rlang .data
#' @export

setMethod("estimateCells", "CIMseqSinglets", function(
                                                      singlets, multiplets, warning = TRUE, maxCellsPerMultiplet = Inf, multiplet.factor=NA) {
  frac.ercc <- NULL
  counts <- cbind(
    getData(singlets, "counts"), 
    getData(multiplets, "counts")
  )
  counts.ercc <- cbind(
    getData(singlets, "counts.ercc"), 
    getData(multiplets, "counts.ercc")
  )
  if(length(counts.ercc)) {
      n.sng <- ncol(getData(singlets, "counts"))
      n.mul <- ncol(getData(multiplets, "counts"))
  
  #check if any samples have ERCC that are all 0
      if(warning) .checkEstimateCellsInput(counts.ercc)
      fe <- colSums(counts.ercc) / (colSums(counts.ercc) + colSums(counts))
      ecn <- median(fe[1:n.sng]) / fe
      ecn[ecn > maxCellsPerMultiplet] <- maxCellsPerMultiplet
  
      return( tibble(
          sample = colnames(counts),
          sampleType = c(rep("Singlet", n.sng), rep("Multiplet", n.mul)),
          frac.ercc = fe, estimatedCellNumber = ecn
      ) )
  } else {
      # No ERCC controls, use mRNA counts directly instead
  n.sng <- ncol(getData(singlets, "counts"))
  n.mul <- ncol(getData(multiplets, "counts"))

  fe.sng <- colSums(counts)[1:n.sng]
  fe.mul <- colSums(counts)[-1:-n.sng]
  calcFrac <- function(x, fe.sng, fe.mul) {
      fe.mul <- fe.mul*x
      ecn.sng <- round(fe.sng / median(fe.sng))
      ecn.mul <- round(fe.mul / median(fe.sng))
      frac.sng <- sum(ecn.sng==0)/sum(ecn.sng == 1)
      frac.mul <- sum(ecn.mul==0)/(sum(ecn.mul==1))*(1+frac.sng*sum(ecn.mul==1)/sum(ecn.mul == 2))
      log(frac.sng/frac.mul)
  }
  if(is.na(multiplet.factor)) {
      m.factors <- seq(0.05, 5, by=0.05)
      m.scores <- sapply(m.factors, calcFrac, fe.sng=fe.sng, fe.mul=fe.mul)
      multiplet.factor <- m.factors[which.min(abs(m.scores))] # Pick the lowest score. min(abs) works since m.scores is logged
  }
  fe <- c(fe.sng, fe.mul*multiplet.factor)
  ecn <- fe / median(fe.sng)
  ecn[ecn > maxCellsPerMultiplet] <- maxCellsPerMultiplet
  
  return( tibble(
    sample = colnames(counts),
    sampleType = c(rep("Singlet", n.sng), rep("Multiplet", n.mul)),
    estimatedCellNumber = ecn
  ) )
  }
})
  
.checkEstimateCellsInput <- function(counts.ercc) {
  all0 <- apply(counts.ercc, 2, function(x) all(x == 0))
  if(any(all0, na.rm = TRUE)) {
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
}
