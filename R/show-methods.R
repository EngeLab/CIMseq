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
            additionalRows <- nrow(mat) - 10
            additionalColumns <- ncol(mat) - 5
            p <- S4Vectors:::makePrettyMatrixForCompactPrinting(
                mat,
                function(x){
                    x[,1:5]
                }
            )
            print(p)
            cat("...\n")
            cat(paste("<", additionalRows, " more elements>", sep=""))
            cat(paste("<", additionalColumns, " more columns>", sep=""))
            cat("\n-----------\n\n")
        } else {
            additionalElements <- length(mat) - 5
            print(head(mat))
            cat("...\n")
            cat(paste("<", additionalElements, " more elements>", sep=""))
            cat("\n-----------\n\n")
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
            additionalRows <- nrow(mat) - 10
            additionalColumns <- ncol(mat) - 5
            p <- S4Vectors:::makePrettyMatrixForCompactPrinting(
            mat,
            function(x){
                x[,1:5]
            }
            )
            print(p)
            cat("...\n")
            cat(paste("<", additionalRows, " more elements>", sep=""))
            cat(paste("<", additionalColumns, " more columns>", sep=""))
            cat("\n-----------\n\n")
        } else {
            additionalElements <- length(mat) - 5
            print(head(mat))
            cat("...\n")
            cat(paste("<", additionalElements, " more elements>", sep=""))
            cat("\n-----------\n\n")
        }
    }
}
