import numpy as np
from pyswarm import pso
import math

"""
   A really simple example with 3 possible cell types and only one gene. 
   The "multuplet" (or slice) value for the gene is the mean of the 2 contributing cell types.
"""
#define fractions, cellTypes, and slice
fractions = [1/np.float64(3)] * 3
cellTypes = np.asarray([25, 75, 0])
slice = [50]

#define function to minimize
def distToSlice(fractions, *args):
    return abs(slice - sum(fractions * cellTypes))

#define constraint; fractions sum to 1
def con1(fractions, *args):
    if sum(fractions) == 1:
        return [0]
    else:
        return [-1]

#test distToSlice with fractions expected after optimization; expect cost to be 0
print distToSlice(np.asarray([0.5, 0.5, 0]), (cellTypes, slice))

#define optimization function
def optimize(fractions, cellTypes, slice):
    lb = np.asarray([0] * len(fractions))
    ub = np.asarray([1] * len(fractions))
    args = (cellTypes, slice)
    xopt, fopt = pso(distToSlice, lb, ub, args=args, f_ieqcons=con1, maxiter=10000, swarmsize=1000)
    return { 'xopt':xopt, 'fopt':fopt }

#run optimization; expected to result in fractions = [0.5, 0.5, 0.0] and cost = 0
result = optimize(fractions, cellTypes, slice)
print result

"""
    This is typically resulting as one of the fractions = 1 and the remaining = 0 even though the cost for this solution is higher than the cost for the desired solution.
"""