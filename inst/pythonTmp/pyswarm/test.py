from pyswarm import pso
import numpy as np

def mydist(x):
    k = len(x)/2
    clusts = np.round(x[0:k], decimals=0).astype(int)
    fracs = x[k:len(x)]
    sumfracs = sum(fracs)
    fracs = fracs / sum(fracs)
    sumsig = np.repeat(0.0, len(celltypes[:,clusts[0]]))
    for i in range(0, k):
        sumsig += celltypes[:,clusts[i]] * fracs[i]
    mydist = np.linalg.norm(dbl1 - sumsig)
    return mydist
    
def banana(x):
    x1 = x[0]
    x2 = x[1]
    return x1**4 - 2*x2*x1**2 + x2**2 + x1**2 - 2*x1 + 5

def con(x):
    x1 = x[0]
    x2 = x[1]
    return [-(x1 + 0.25)**2 + 0.75*x2]

k = 4
nclusts = 12
ngenes = 500
maxiter = 10000
swarmsize = 1000

#f = open ( 'test_celltypes.txt' , 'r')
f = open ( 'top2000_genes_in_celltypes.txt' , 'r')
celltypes = []
celltypes = np.array([ map(float, line.split("\t")) for line in f])
celltypes = celltypes[0:ngenes,:]
print celltypes

f = open ( 'top2000_genes_in_doubles.txt' , 'r')
doubles = []
doubles = np.array([ map(float, line.split("\t")) for line in f])
doubles = doubles[0:ngenes,:]
print doubles
fout  = open('result.csv', 'w')
dbl1 = doubles[:,1]

clustRes = range(0, k)
fracRes = range(0, k)

#for i in range(0, 5):
for i in range(0, len(doubles[0,:])):
    dbl1 = doubles[:,i]
    lb = np.repeat([0,.0001], k)
    ub = np.repeat([nclusts-1, 1], k)
#    print lb, "\n", ub, "\n"
    xopt, fopt = pso(mydist, lb, ub, f_ieqcons=con, swarmsize=swarmsize, minstep=1e-10, minfunc=1e-10, maxiter=maxiter)
    print xopt, "\n"
    print fopt, "\n"
    k = len(xopt)/2
    clusts = np.round(xopt[0:k], decimals=0).astype(int)
    fracs = xopt[k:len(xopt)]
    fracs /= sum(fracs)

    clustRes = np.vstack([clustRes, clusts])
    for clust in clusts:
        print >>fout, clust, "\t"
        
    fracRes = np.vstack([fracRes, fracs])
    for frac in fracs:
        print >>fout, frac, "\t"

    print >>fout, "\n"
    print clusts, "\n"
    print fracs, "\n\n"

np.savetxt("clusts.csv", clustRes, delimiter=",", fmt='%i')
np.savetxt("fracs.csv", fracRes, delimiter=",", fmt='%.3f')

#xopt, fopt = pso(banana, lb, ub, f_ieqcons=con, swarmsize=10000, processes=10)





