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
    norm_fracts = fractions / sum(fractions) # Normalize instead of constrain
    return abs(slice - sum(norm_fracts * cellTypes))


#define constraint; fractions do not sum to 0
def con0(fractions, *args):
    if sum(fractions) > 0.01: # Any small value > 0
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
    xopt, fopt = pso(distToSlice, lb, ub, f_ieqcons=con0, args=args, maxiter=1000, swarmsize=250, minstep=1e-16, minfunc=1e-16)
    xopt = xopt / sum(xopt); # Normalize result.
    return { 'xopt':xopt, 'fopt':fopt }


#run optimization; expected to result in fractions = [0.5, 0.5, 0.0] and cost = 0
result0 = optimize(fractions, cellTypes, slice)
print result0;

"""
    result0 should end up with the right solution most of the time.
    Small minstep and minfunc are for robustness (also, there are many 'decent' solutions to this optimization).
"""
