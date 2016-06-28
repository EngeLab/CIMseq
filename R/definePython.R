
#'@include All-classes.R
NULL

#' pySwarm
#'
#' Subtitle
#'
#' Imports count and count.ercc data to a sp.scRNAseq object.
#'
#' @name pySwarm
#' @rdname pySwarm
#' @aliases pySwarm
#' @param counts Counts matrix with samples as columns and genes as rows.
#' @param ... additional arguments to pass on
#' @return The spCounts function returns an object of class spCounts.
#' @author Jason T. Serviss
#' @keywords pySwarm
#' @examples
#'
#' #use demo data
#'
#'
#' #run function
#'
#'
NULL

#' @rdname pySwarm
#' @export

setGeneric("pySwarm", function(spCounts, ...
){ standardGeneric("spCounts") })

#' @importFrom rPython python.exec python.assign python.get
#' @rdname spCounts
#' @export
setMethod("spCounts", "matrix",
function(
    counts,
    ...
){
    cellTypes <-
    slice <-
    fractions <-
    
    ##run pySwarm
    .defineImport(cellTypes, slice, fractions)
    .defineCost()
    .definePySwarm()
    
    cmd1 <- 'result = optimize(cellTypes, slice, fractions)'
    python.exec(cmd1)
    result <- python.get('result')
    return(result)
})



#python.exec('print type(py_dict1)')

.defineImport <- function(cellTypes, slice, fractions) {
    python.exec('import pandas as pd')
    python.exec('import numpy as np')
    
    python.assign('cellTypesDictionary', as.data.frame(cellTypes))
    python.assign('sliceDictionary', as.data.frame(slice))
    python.assign('fractions', fractions)
    
    python.exec('cellTypes = pd.DataFrame(cellTypesDictionary)')
    python.exec('slice = pd.DataFrame(sliceDictionary)')
    python.exec('fractions = np.asarray(fractions)')
}

.defineCost <- function() {
    python.exec('import pandas as pd')
    python.exec('import numpy as np')
    python.exec('import math')
    
    cmd1 <- 'def makeSyntheticSlice(cellTypes, fractions):
        s = sum(fractions)
    
        for i, f in enumerate(fractions):
            fractions[i] = np.float64(fractions[i]) / s
    
        func = lambda x: sum(np.asarray(x) * np.asarray(fractions))
        return cellTypes.apply(func, axis=1)'
    
    cmd2 <- 'def distToSlice(fractions, *args):
        cellTypes, slice = args
        a = makeSyntheticSlice(cellTypes, fractions)
        for i in a:
            if math.isnan(i):
                print "NaN in make synthetic slice!"
                print a
    
        diff = []
        for index, row in cellTypes.iterrows():
            d = float(a.iloc[[index]]) - float(slice.iloc[[index]])
            diff.append(abs(d))
        cost = sum(diff)
        return cost'
    
    python.exec(cmd1)
    python.exec(cmd2)
}


.definePySwarm <- function() {
    python.exec('from pyswarm import pso')
    
    #cmd1 <- 'def constraints(fractions, args*):'
    
    cmd2 <- 'def optimize(cellTypes, slice, fractions):
        lb = [0, 0, 0]
        ub = [1, 1, 1]
        args = (cellTypes, slice)
        xopt, fopt = pso(distToSlice, lb, ub, args=args)
        ##!!!add return statment'
    
    #python.exec(cmd1)
    python.exec(cmd2)
}


##set up example data

s1 = [x for x in range(1,11)]
s2 = [x*10 for x in range(1,11)]
s3 = [x*100 for x in range(1,11)]

exp = {
    's1': s1,
    's2': s2,
    's3': s3
}

cellTypes = pd.DataFrame(exp, columns=exp.keys())
slice = pd.Series(exp['s1'])
fractions = [0.25, 0.5, 0.25]








