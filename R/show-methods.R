#'@include All-classes.R
NULL

#' @rdname spCounts
#' @importFrom S4Vectors show
#' @importFrom utils head
#' @export
setMethod("show", "spCounts", function(
    object
){
    .showCounts(object)
})

#internal show function
.showCounts <- function(
    object
){
    names <- slotNames(object)
    cat("Class:","spCounts\n")
    cat("Contains: \n")
    for(i in 1:length(names)){
        cat(paste(i,". ", names[i], "\n",sep=""))
        mat <- slot(object, names[i])

        if(class(mat) == "matrix") {
            .showMatrix(mat)
        } else {
            .showBasics(mat)
        }
    }
}

#' @rdname spUnsupervised
#' @importFrom S4Vectors show
#' @importFrom utils head
#' @export
setMethod("show", "spUnsupervised", function(
    object
){
    .showUnsupervised(object)
})

#internal show function
.showUnsupervised <- function(
    object
){
    names <- slotNames(object)
    cat("Class:","spUnsupervised\n")
    cat("Contains: \n")
    for(i in 1:length(names)){
        
        cat(paste(i, ". ", names[i], "\n",sep=""))
        mat <- slot(object, names[i])
        
        if(class(mat) == "matrix" | class(mat) == "data.frame") {
            .showMatrix(mat)
        } else if(class(mat) == "list") {
            .showList(mat)
        } else {
            .showBasics(mat)
        }
    }
}

#' @rdname spSwarm
#' @importFrom S4Vectors show
#' @importFrom utils head
#' @export
setMethod("show", "spSwarm", function(
    object
){
    .showSpSwarm(object)
})

#internal show function
.showSpSwarm <- function(
    object
){
    names <- slotNames(object)
    cat("Class:","spSwarm\n")
    cat("Contains: \n")
    for(i in 1:length(names)){
        
        cat(paste(i,". ", names[i], "\n",sep=""))
        mat <- slot(object, names[i])
        
        if(class(mat) == "data.frame") {
            .showMatrix(mat)
        } else if(class(mat) == "list") {
            .showList(mat)
        } else if(class(mat) == "spCounts" | class(mat) == "spUnsupervised") {
            show(mat)
        } else {
            .showBasics(mat)
        }
    }
}

################################################################################
#                                                                              #
#                                   Units                                      #
#                                                                              #
################################################################################

.showBasics <- function(
    obj
){
    
    if(class(obj) == "call") {return("")}
    
    if(length(obj) <= 5) {
        cat(head(obj))
        cat("\n-----------\n\n")
    } else {
        additionalElements <- length(obj) - 5
        cat(head(obj))
        cat("...\n")
        cat(paste("<", additionalElements, " more elements>", sep=""))
        cat("\n-----------\n\n")
    }
}


.showList <- function(
    obj
){
    show(as(obj, "List"))
}

.showMatrix <- function(
    obj
){
    if(all(is.na(obj)) == TRUE) {
        
        print("NA")
        cat("\n-----------\n\n")
        
    } else {
        
        cat(paste("<", nrow(obj), " elements>", sep=""))
        cat(paste("<", ncol(obj), " columns>", sep=""))
        cat("\n-----------\n\n")
    }

}