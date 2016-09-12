#'@include All-classes.R
NULL

#' @rdname spCounts
#' @export
setMethod("show", "spCounts", function(object){ .showCounts(object) })

#internal show function
.showCounts <- function(object
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
#' @export
setMethod("show", "spUnsupervised", function(object){ .showUnsupervised(object) })

#internal show function
.showUnsupervised <- function(object
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
#' @export
setMethod("show", "spSwarm", function(object){ .showSpSwarm(object) })

#internal show function
.showSpSwarm <- function(object
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

####################################################################################
#                                                                                  #
#                                   Units                                          #
#                                                                                  #
####################################################################################

.showBasics <- function(obj) {
    
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


.showList <- function(obj) {
    S4Vectors::show(as(obj, "List"))
}


.showMatrix <- function(obj) {
    
    if(all(is.na(obj)) == TRUE) {
        
        print("NA")
        cat("\n-----------\n\n")
        
    } else {
        
        if(nrow(obj) <= 10) {
            additionalRows <- 0
        } else {
            additionalRows <- nrow(obj) - 10
        }
        
        if(ncol(obj) <= 5) {
            additionalColumns <- 0
        } else {
            additionalColumns <- ncol(obj) - 5
        }
        
        p <- S4Vectors:::makePrettyMatrixForCompactPrinting(
            obj,
            function(x){
                x[,1:2]
            }
        )
        print(p)
        cat("...\n")
        cat(paste("<", additionalRows, " more elements>", sep=""))
        cat(paste("<", additionalColumns, " more columns>", sep=""))
        cat("\n-----------\n\n")
    }
}

