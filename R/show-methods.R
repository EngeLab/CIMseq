#'@include All-classes.R
NULL

#' @rdname CIMseqSinglets
#' @importFrom utils head
#' @export

setMethod("show", "CIMseqSinglets", function(object){
  .showCIMseqData(object)
})

#' @rdname CIMseqMultiplets
#' @importFrom S4Vectors show
#' @importFrom utils head
#' @export

setMethod("show", "CIMseqMultiplets", function(object){
  .showCIMseqData(object)
})

#internal show function
.showCIMseqData <- function(object){
  names <- slotNames(object)
  cat("Class:", class(object)[[1]], "\n")
  cat("Contains: \n\n")
  for(i in 1:length(names)){
    cat(paste(i, ". ", names[i], "\n", sep = ""))
    obj <- slot(object, names[i])
    if(is.matrix(obj)) .showMatrix(obj)
    if(is.function(obj)) .showFunction(obj)
    if((is.character(obj) | is.numeric(obj)) & !is.matrix(obj)) .showBasics(obj)
  }
}

#' @rdname CIMseqSwarm
#' @importFrom S4Vectors show
#' @importFrom utils head
#' @export

setMethod("show", "CIMseqSwarm", function(object){
  .showCIMseqSwarm(object)
})

#internal show function
.showCIMseqSwarm <- function(object){
  names <- slotNames(object)
  cat("Class:", class(object)[[1]], "\n")
  cat("Contains: \n\n")
  for(i in 1:length(names)){
    cat(paste(i,". ", names[i], "\n",sep=""))
    mat <- slot(object, names[i])
    if(is.list(mat) & !is_tibble(mat)) .showList(mat)
    if(is.matrix(mat)) .showMatrix(mat)
    if(is_tibble(mat)) .showTibble(mat)
    if((is.character(mat) | is.numeric(mat)) & !is.matrix(mat)) .showBasics(mat)
  }
}

################################################################################
#                                                                              #
#                                   Units                                      #
#                                                                              #
################################################################################

.showBasics <- function(obj){
  if(class(obj) == "call") {return("")}
  if(length(obj) <= 5) {
    cat(head(obj))
    cat("\n-----------\n\n")
  } else {
    additionalElements <- length(obj) - 5
    cat(head(obj, n = 3))
    cat("...\n")
    cat(paste("<", additionalElements, " more elements>", sep=""))
    cat("\n-----------\n\n")
  }
}

.showList <- function(obj){
  print(paste0("List of length ", length(obj)))
  cat("-----------\n\n")
}

.showMatrix <- function(obj){
  if(all(is.na(obj)) == TRUE) {
    print("NA")
    cat("\n-----------\n\n")
  } else {
    nc <- min(3, ncol(obj))
    nr <- min(3, nrow(obj))
    cat(paste("<", nrow(obj), " elements>", sep=""))
    cat(paste("<", ncol(obj), " columns>", sep=""))
    cat("\n")
    print(obj[1:nr, 1:nc])
    cat("\n-----------\n\n")
  }
}

.showTibble <- function(obj) {
  print(obj)
  cat("\n-----------\n\n")
}

.showFunction <- function(obj) {
  print(obj)
  cat("\n-----------\n\n")
}
