import numpy as np
from math import comb
import itertools
import pandas as pd


def DasDennis(p,dim):
    """create a equidistance n dimensional spacing which satisfies the mass balance constraint
    
    Examples:
        >>> p=30
        >>> dim=3
        >>> spacevec=DasDennis(p, dim)
        >>> pd.DataFrame(spacevec.T).to_excel("test.xlsx")
    """
    co=(p+dim-2)/p
    nco=int(round(co*p,0))
    covec=np.linspace(0,co,nco+1)
    complexity=comb(dim+p-1,p)
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
    spacevec[spacevec<=1E-8]=0
    return spacevec
