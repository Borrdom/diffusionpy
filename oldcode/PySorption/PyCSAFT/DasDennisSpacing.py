import numpy as np
from scipy.special import comb
import itertools



def DasDennis(p,dim):
    co=(p+dim-2)/p
    nco=int(round(co*p,0))
    covec=np.linspace(0,co,nco+1)
    complexity=comb(dim+p-1,p,exact=True)
    combinations=np.asarray(list(itertools.combinations(covec, dim-1))).T
    for i in range(dim-1):
        for j in range(complexity):
            combinations[i,j]=combinations[i,j]-i/p
    Nlinspaces=complexity//p
    spacevec=np.zeros((dim,complexity))
    spacevec[0,:]=combinations[0,:]
    for i in range(1,dim-1):
        spacevec[i,:]=combinations[i,:]-combinations[i-1,:]
    spacevec[dim-1,:]=1-combinations[dim-2,:]    

    return spacevec

if __name__=="__main__":
    p=30
    dim=3
    spacevec=DasDennis(p, dim)