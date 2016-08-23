import pandas as pd
import numpy as np
import math
from pyswarm import pso
from collections import Counter

#import data
with open('/Users/jasonserviss/Desktop/cellTypes.txt') as f:
    cellTypes = pd.read_table(f, sep='\t', index_col=0, header=0, lineterminator='\n')

with open('/Users/jasonserviss/Desktop/slice.txt') as f:
    slice = pd.read_table(f, sep='\t', index_col=0, header=0, lineterminator='\n')

fractions = []
for l in range(len(cellTypes.columns)):
    fractions.append(1/np.float64(len(cellTypes.columns)))

#define function to minimize
def makeSyntheticSlice(cellTypes, fractions):
    """s = sum(fractions)
    
    for i, f in enumerate(fractions):
        fractions[i] = np.float64(fractions[i]) / s"""
    
    func = lambda x: sum(np.asarray(x) * np.asarray(fractions))
    return cellTypes.apply(func, axis=1)

def distToSlice(fractions, *args):
    cellTypes, slice, col = args
    a = makeSyntheticSlice(cellTypes, fractions)
    for i in a:
        if math.isnan(i):
            print "NaN in make synthetic slice!"
            print a
    
    diff = []
    for index in range(cellTypes.shape[0]):
        d = a[index] - slice.iloc[index, col]
        diff.append(abs(d))
    cost = sum(diff)
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

#with real data you would want to add a constraint on non-null values, i.e. a multuplet should be at least x cells

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

#test
#fourSeven
result0 = optimize(fractions, cellTypes, slice, 0)
#fourNine
result1 = optimize(fractions, cellTypes, slice, 1)
#sevenNine
result2 = optimize(fractions, cellTypes, slice, 2)
#fourSevenNine
result3 = optimize(fractions, cellTypes, slice, 3)

#run optimization
"""xopt = {}
fopt = {}
for col in range(len(slice.columns)):
    print "analyzing multuplet number: " + str(col)
    result = optimize(fractions, cellTypes, slice, col)
    xopt = {col: result['xopt']}
    fopt = {col: result['fopt']}"""


##################################
#try with minimize
from scipy.optimize import minimize

cons = (
{
    'type': 'eq',
    'fun': lambda fractions: 1 - sum(fractions)
},
{
    'type': 'eq',
    'fun': lambda fractions: sum([Counter(fractions)[k] for k in list(set([i for i in fractions if i < 1/np.float64(len(fractions)) and i != 0.0]))])
})

args = (cellTypes, slice, 5)
bnds = ((0, 1), (0, 1), (0, 1))
res = minimize(distToSlice, [0, 0, 0], args=args, constraints=cons, method='SLSQP', options={'disp': True, 'maxiter': 100}, bounds=bnds, tol=1e-10)

##################################
#try with pyOpt
import pyOpt

def objfunc(fractions, *args):
    f = distToSlice(fractions, *args)
    g = constraints(fractions, *args)
    fail = 0
    return f,g, fail

opt_prob = pyOpt.Optimization('test', objfunc)

#Define Problem Design Variables (3 Continuous Variables)
lb = [0] * len(fractions)
up = [1] * len(fractions)
opt_prob.addVarGroup('fractions', len(fractions), type='c', lower=lb, upper=up)

#add name of variable in objective function which includes the function to be minimized's output
opt_prob.addObj('f', optimum=0.0)

#add constraints
opt_prob.addConGroup('name', len(fractions), type='i', lower=lb, upper=up, equal='i')

#initialize optimizers
slsqp = pyOpt.SLSQP()

#optimize
args = (cellTypes, slice, 5)
[fstr, xstr, inform] = slsqp(opt_prob, *args)



