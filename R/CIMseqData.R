#' @include All-classes.R
NULL

#' CIMseqSinglets
#'
#' Subtitle
#'
#' Imports count, dimensionality reduction, and classification data
#' to a CIMseqSinglets object for sequenced singlets.
#'
#' @name CIMseqSinglets
#' @rdname CIMseqSinglets
#' @param counts matrix; Counts matrix with samples as columns and genes as rows.
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
    dim.red = matrix(nrow = 0, ncol = 0),
    classification = character(),
    ...
  )
})

#' @rdname CIMseqSinglets
#' @export

setMethod("CIMseqSinglets", "matrix", function(
  counts, dim.red, classification, ...
  ){
  .inputCheckSinglets(counts, dim.red, classification)
  new(
    "CIMseqSinglets",
    counts = counts,
    counts.log = .norm.log.counts,
    counts.cpm = .norm.counts,
    dim.red = dim.red,
    classification = classification,
    ...
  )
})

.inputCheckSinglets <- function(counts, dim.red, classification) {
  if(ncol(counts) != length(classification)) {
    message("length(classification) != ncol(counts)")
  }
  if(nrow(dim.red) != ncol(counts)) {
    message("nrow(dim.red) != ncol(counts)")
  }
}

.norm.log.counts <- function(counts) {
  log2(.norm.counts(counts) + 1)
}

# Hack, to cater for lower numbers in 10x genomics.
.norm.counts <- function(counts, norm.factor=10000) { 
    if(is.na(norm.factor)) {
        norm.factor <- mean(colSums(counts))
    }
  t(t(counts) / colSums(counts) * norm.factor)
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
    features = numeric(),
    ...
  )
})

#' @rdname CIMseqMultiplets
#' @export

setMethod("CIMseqMultiplets", "matrix", function(
  counts, features, ...
  ){
        
  .inputCheckMultiplets(counts)
  new(
    "CIMseqMultiplets",
    counts = counts,
    counts.log = .norm.log.counts,
    counts.cpm = .norm.counts,
    features = features,
    ...
  )
})

.inputCheckMultiplets <- function(counts) {
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

# Changed to use UMI counts.
setMethod("estimateCells", "CIMseqSinglets", function(
  singlets, multiplets, warning = TRUE, maxCellsPerMultiplet = Inf, multiplet.factor=NA){
  counts <- cbind(
    getData(singlets, "counts"), 
    getData(multiplets, "counts")*median(colSums(getData(singlets, "counts")))/median(colSums(getData(multiplets, "counts")))
  )
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
#  calcFrac <- function(x, fe.sng, fe.mul) {
#      fe.mul <- fe.mul*x
#      ecn.sng <- round(fe.sng / median(fe.sng))
#      ecn.mul <- round(fe.mul / median(fe.sng))
#      log((sum(ecn.sng==0)/sum(ecn.sng!=0))/(sum(ecn.mul==0)/sum(ecn.mul!=0)))
#  }
  if(is.na(multiplet.factor)) {
      m.factors <- seq(0.05, 5, by=0.05)
      m.scores <- sapply(m.factors, calcFrac, fe.sng=fe.sng, fe.mul=fe.mul)
      multiplet.factor <- m.factors[which.min(abs(m.scores))] # Pick the lowest score. min(abs) works since m.scores is logged
  }
  fe <- c(fe.sng, fe.mul*multiplet.factor)
  ecn <- fe / median(fe.sng)
  ecn[ecn > maxCellsPerMultiplet] <- maxCellsPerMultiplet
  
  tibble(
    sample = colnames(counts),
    sampleType = c(rep("Singlet", n.sng), rep("Multiplet", n.mul)),
    estimatedCellNumber = ecn
  )

# Older version
#    counts <- cbind(
#    getData(singlets, "counts"), 
#    getData(multiplets, "counts")*2*median(colSums(getData(singlets, "counts")))/median(colSums(getData(multiplets, "counts"))) # Normalize by multiplet/singlet batch, adjust since median multiplet count is 2
#  )
#  n.sng <- ncol(getData(singlets, "counts"))
#  n.mul <- ncol(getData(multiplets, "counts"))
  
#  fe <- colSums(counts)
#  ecn <- fe / median(fe[1:n.sng])
#  ecn[ecn > maxCellsPerMultiplet] <- maxCellsPerMultiplet
  
#  tibble(
#    sample = colnames(counts),
#    sampleType = c(rep("Singlet", n.sng), rep("Multiplet", n.mul)),
#    estimatedCellNumber = ecn
#  )
})

