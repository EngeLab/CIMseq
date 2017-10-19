
#'@include All-classes.R
NULL

#' spCounts
#'
#' Subtitle
#'
#' Imports count, sampleType, and count.ercc data to a sp.scRNAseq object.
#'
#' @name spCounts
#' @rdname spCounts
#' @aliases spCounts
#' @param counts Counts matrix with samples as columns and genes as rows.
#' @param counts.ercc A matrix containing ercc spike-in reads and their counts.
#' @param counts.cpm Normalized counts per million.
#' @param counts.log Log2 normalized counts per million.
#' @param object spCounts object.
#' @param n Data to extract from spCounts object.
#' @param .Object Internal object.
#' @param ... additional arguments to pass on
#' @return spCounts object.
#' @author Jason T. Serviss
#' @keywords spCounts
#' @examples
#'
#' s <- grepl("^s", colnames(testCounts))
#' cObj <- spCounts(testCounts[, s], testErcc[, s])
#'
NULL

#' @rdname spCounts
#' @export

setGeneric("spCounts", function(
    counts,
    ...
){
    standardGeneric("spCounts")
})

#' @rdname spCounts
#' @export
setMethod("spCounts", "matrix", function(
    counts,
    counts.ercc,
    ...
){
    .inputCheckCounts(counts, counts.ercc)
    new("spCounts",
        counts = counts,
        counts.log = .norm.log.counts(counts),
        counts.cpm = .norm.counts(counts),
        counts.ercc = counts.ercc
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
    norm.fact <- colSums(counts)
    counts.norm <- t(apply(counts, 1, .norm, n = norm.fact))
    counts.log <- log2(counts.norm)
}

.norm.counts <- function(counts) {
    norm.fact <- colSums(counts)
    counts.cpm <- t(apply(counts, 1, .norm, n = norm.fact))
}

.norm <- function(x, n) {
    x / n * 1000000 + 1
}

#.deconv <- function(counts, counts.ercc) {
#    D <- rbind(counts, counts.ercc)
#    sce_deconv<-newSCESet(countData=data.frame(D))
#    sce_deconv<-calculateQCMetrics(sce_deconv,feature_controls=list(MySpikes=c(rep(FALSE, nrow(counts)), rep(TRUE, nrow(counts.ercc)))))
#    setSpike(sce_deconv) <- "MySpikes"
#    #spikes(sce_deconv)[1:5,1:5] #extracts matrix of spikes
#
#    clusters <- quickCluster(sce_deconv,get.spikes=FALSE,min.size=10)
#    sce_deconv <- computeSumFactors(sce_deconv)
#    summary(sizeFactors(sce_deconv))
#
#    sce_deconv <- computeSpikeFactors(sce_deconv,general.use=FALSE)
#    summary(sizeFactors(sce_deconv))
#
#    sce_deconv<-normalize(sce_deconv)
#
#    as.data.frame(t(apply(D[grepl("ERCC_",rownames(D))==FALSE,],1,"/",sce_deconv@phenoData@data$size_factor)))
#    #expr_deconv_ercc<-as.data.frame(t(apply(counts.ercc,1,"/",sce_deconv@phenoData@data$size_factor_MySpikes)))
#}

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
#' @aliases estimateCells
#' @param spCountsSng A spCounts object with singlets.
#' @param spCountsMul A spCounts object with multiplets.
#' @param ... additional arguments to pass on
#' @return A data frame including the fraction of ercc reads and cell counts for
#'    each sample.
#' @author Jason T. Serviss
#' @keywords spCounts
#' @examples
#'
#' #use test data
#' s <- grepl("^s", colnames(testCounts))
#' cObjSng <- spCounts(testCounts[, s], testErcc[, s])
#' cObjMul <- spCounts(testCounts[, !s], testErcc[, !s])
#'
#' #run function
#' output <- estimateCells(cObjSng, cObjMul)
#'
NULL

#' @export

setGeneric("estimateCells", function(
    spCountsSng,
    spCountsMul,
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
setMethod("estimateCells", "spCounts", function(
    spCountsSng,
    spCountsMul,
    ...
){
    counts <- cbind(
        getData(spCountsSng, "counts"),
        getData(spCountsMul, "counts")
    )
    
    counts.ercc <- cbind(
        getData(spCountsSng, "counts.ercc"),
        getData(spCountsMul, "counts.ercc")
    )
    
    d <- tibble(
        sampleName = c(
            colnames(getData(spCountsSng, "counts")),
            colnames(getData(spCountsMul, "counts"))
        ),
        sampleType = c(
            rep(
                "Singlet",
                ncol(getData(spCountsSng, "counts"))
            ),
            rep(
                "Multiplet",
                ncol(getData(spCountsMul, "counts"))
            )
        ),
        frac.ercc = colSums(counts.ercc) /
            (colSums(counts.ercc) + colSums(counts))
    )
    
    d$cellNumberMin <-
        d %>%
            filter(.data$sampleType == "Singlet") %>%
            pull(.data$frac.ercc) %>%
            quantile %>%
            `[` (2) %>%
            `/` (d$frac.ercc)
    
    d$cellNumberMedian <-
        d %>%
            filter(.data$sampleType == "Singlet") %>%
            pull(.data$frac.ercc) %>%
            median %>%
            `/` (d$frac.ercc)
    
    d$cellNumberMax <-
        d %>%
            filter(.data$sampleType == "Singlet") %>%
            pull(.data$frac.ercc) %>%
            quantile %>%
            `[` (4) %>%
            `/` (d$frac.ercc)
    
    return(d)
})
