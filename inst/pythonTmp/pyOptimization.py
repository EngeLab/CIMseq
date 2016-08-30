import pandas as pd
import numpy as np
import math
from pyswarm import pso
from collections import Counter
from random import randint

#import data
with open('/Users/jasonserviss/Github/sp.scRNAseq/inst/pythonTmp/cellTypes.txt') as f:
    cellTypes = pd.read_table(f, sep='\t', index_col=0, header=0, lineterminator='\n')

with open('/Users/jasonserviss/Github/sp.scRNAseq/inst/pythonTmp/slice.txt') as f:
    slice = pd.read_table(f, sep='\t', index_col=0, header=0, lineterminator='\n')

fractions = []
for l in range(len(cellTypes.columns)):
    fractions.append(1/np.float64(len(cellTypes.columns)))

#define function to minimize

def makeSyntheticSlice(cellTypes, fractions):
    func = lambda x: sum(np.asarray(x) * np.asarray(fractions))
    return cellTypes.apply(func, axis=1)

def distToSlice(fractions, *args):
    cellTypes, slice, col = args
    normFractions = fractions / sum(fractions)
    a = makeSyntheticSlice(cellTypes, normFractions)
    
    for i in a:
        if math.isnan(i):
            print "NaN in make synthetic slice!"
            print a
    
    diff = []
    for index in range(cellTypes.shape[0]):
        d = a.iloc[index,] - slice.iloc[index, col]
        diff.append(abs(d))
    cost = sum(diff)
    return cost

#define constraints
"""f_ieqcons : function
    Returns a 1-D array in which each element must be greater or equal to 0.0 in a successfully optimized problem. If f_ieqcons is specified, ieqcons is ignored (Default: None)"""
#fractions sum to 1
def con1(fractions, *args):
    if sum(fractions) > 0.1:
        return [0]
    else:
        return [-1]

#with real data you would want to add a constraint on non-null values, i.e. a multuplet should be at least x cells

#define optimization function
def optimize(fractions, cellTypes, slice, col):
    lb = np.asarray([0] * len(fractions))
    ub = np.asarray([1] * len(fractions))
    args = (cellTypes, slice, col)
    xopt, fopt = pso(distToSlice, lb, ub, args=args, f_ieqcons=con1, maxiter=10, swarmsize=150)
    dictionary = dict(zip(list(cellTypes.columns.values), xopt.tolist()))
    return { 'xopt':dictionary, 'fopt':fopt }

#test
"""
#fourSeven
result0 = optimize(fractions, cellTypes, slice, 0)
#fourNine
result1 = optimize(fractions, cellTypes, slice, 1)
#sevenNine
result2 = optimize(fractions, cellTypes, slice, 2)
#fourSevenNine
result3 = optimize(fractions, cellTypes, slice, 3)
"""

#run optimization
sampleNames = slice.columns.values
df = pd.DataFrame(columns=['sample', 'xopt', 'fopt'])
r = slice.shape[1]

for y in range(10):
    s = randint(0,r)
    currentSample = sampleNames[s]
    result = optimize(fractions, cellTypes, slice, s)
    tmp = pd.DataFrame({'sample':currentSample, 'xopt': [result['xopt']], 'fopt': result['fopt']})
    df = df.append(tmp)

df.to_csv("~/Desktop/seqOptResults.txt", sep='\t')
