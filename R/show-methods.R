#'@include All-classes.R
NULL

#' @rdname spCounts
#' @export
setMethod("show", "spCounts", function(object){ .showCounts(object) })

#internal show function
.showCounts <- function(object
){
    names <- slotNames(object)
    cat("Class:","Counts\n")
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
    cat("Class:","Counts\n")
    cat("Contains: \n")
    for(i in 1:length(names)){
        
        cat(paste(i,". ", names[i], "\n",sep=""))
        mat <- slot(object, names[i])
        
        if(class(mat) == "matrix") {
            .showMatrix(mat)
        } else if(class(mat) == "list") {
            .showList(mat)
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

#.showList <- function(obj) {
#    for(oo in 1:length(obj)) {
#        curr <- obj[[oo]]
#        name <- names(obj[oo])
#        cat(paste(name, ": ", sep=""))
#
#        if(class(curr) == "matrix" | class(curr) == "mclustBIC") {
#            .showMatrix(curr)
#        } else if(class(curr) == "list") {
#            .showList(curr)
#        } else {
#            .showBasics(curr)
#       }
#    }
#}

.showList <- function(obj) {
    leng <- length(obj)
    cat(paste("List with ", leng, " elements", sep=""))
}

.showMatrix <- function(obj) {
    
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
