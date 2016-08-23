import numpy as np
from pyswarm import pso

fractions = []
for l in range(0, 3):
    fractions.append(1/np.float64(3))
fractions = np.asarray(fractions)

cellTypes = np.asarray([50, 50, 0])
slice = [100]


#define function to minimize
def makeSyntheticSlice(cellTypes, fractions):
    return sum(cellTypes * fractions)

def distToSlice(fractions, *args):
    cellTypes, slice = args
    a = makeSyntheticSlice(cellTypes, fractions)
    if math.isnan(a):
        print "NaN in make synthetic slice!"
        print a
    
    cost = abs(a - slice)
    return cost


#define constraints
"""f_ieqcons : function
    Returns a 1-D array in which each element must be greater or equal to 0.0 in a successfully optimized problem. If f_ieqcons is specified, ieqcons is ignored (Default: None)"""

#fractions sum to 1
def con1(fractions, *args):
    if sum(fractions) == 1:
        return 0
    else:
        return -1

#constrains each fraction to be 0 or > the fractions if all cell types were divided equally among the multuplet
def con2(fractions, *args):
    limit = 1/np.float64(len(fractions))
    oc1 = list(set([i for i in fractions if i < limit and i != 0.0]))
    if sum([Counter(fractions)[k] for k in oc1]) == 0:
        return 0
    else:
        return -1

def constraints(fractions, *args):
    cons1 = con1(fractions, *args)
    cons2 = con2(fractions, *args)
    return [cons1, cons2]

#define optimization function
def optimize(fractions, cellTypes, slice, col):
    lb = np.asarray([0] * len(fractions))
    ub = np.asarray([1] * len(fractions))
    args = (cellTypes, slice, col)
    xopt, fopt = pso(distToSlice, lb, ub, args=args, f_ieqcons=constraints, maxiter=500)
    dictionary = dict(zip(list(cellTypes.columns.values), xopt.tolist()))
    return { 'xopt':dictionary, 'fopt':fopt }


result0 = optimize(fractions, cellTypes, slice, 0)

