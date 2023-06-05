import numpy as np
from numba import njit
import matplotlib.pyplot as plt
def MEOS_mode(ode,EJ, etaJ, exponent,RV,v2):
    def MEOS_ode(t,x,THFaktor,dmuext,rhoiB,drhovdtB):
        _,nz=dmuext.shape
        nTH=drhovdtB.shape[0]
        nJ=len(EJ)
        rhov=np.zeros((nTH,nz+1))
        for i in range(nTH):
            rhovtemp=x[(nz+1)*(i):(nz+1)*(1+i)]
            rhov[i,:]=rhovtemp
        sigmaJ=np.zeros((nz+1,nJ))
        for J in range(nJ):
            sigmaJtemp=x[(nz+1)*(nTH+J):(nz+1)*(nTH+1+J)]#*np.atleast_1d(EJ)[J]
            sigmaJ[:,J]=sigmaJtemp
        X2II=rhov/np.sum(rhoiB)
        w2II=X2II/(X2II+1)
        WL=np.prod(np.exp(-w2II*np.expand_dims(exponent,-1)),0)
        etaWL=np.expand_dims(WL,1)*np.expand_dims(etaJ,0)
        rhoiB[0]=rhov[:,-1]
        rhov=np.ascontiguousarray(rhov)
        dmuext=MDF(sigmaJ,EJ,RV)
        drhovdt=ode(t,rhov,THFaktor,dmuext,rhoiB,drhovdtB)
        dsigmaJdt=stress(etaWL,EJ,sigmaJ,drhovdt,v2)
        drhovdt[:,-1]=drhovdtB
        fvec=np.hstack((drhovdt.flatten(),dsigmaJdt.flatten()))
        return fvec
    return MEOS_ode

#@njit
def MDF(sigmaJ,EJ,RV):
    sigma=np.sum(sigmaJ*EJ,axis=1)
    dsigma=np.diff(sigma)
    dmuext=1/np.expand_dims(RV,1)*np.expand_dims(dsigma,0)
    return dmuext

#@njit
def stress(etaWL,EJ,sigmaJ,drhodtNF,v2):
    nJ=len(EJ)
    nz_1,_=sigmaJ.shape
    dsigmaJdt=np.zeros((nz_1,nJ))
    for i in range(nJ):
        dsigmaJdt[:,i]=-1/etaWL[:,i]*sigmaJ[:,i]*EJ[i]+np.sum(drhodtNF,axis=0)*v2
    return dsigmaJdt

def boundary(rhov,etaWL,EJ,sigmaJ,RV,THFaktor,v2):
    b=np.sum(1/(etaWL[-1,:])*sigmaJ[-1,:]*EJ**2,axis=0)/RV
    A=(THFaktor/rhov[:,-1]+np.sum(EJ)*v2/RV)
    drhovdtB=np.linalg.solve(A,b) if len(A)>1 else b/A
    return drhovdtB

def initialboundary(RV,EJ,v2,THFaktors,rho,mobiles,wi0,wi8):
    from scipy.special import lambertw
    B=RV/EJ/v2*np.diag(np.mean(THFaktors,axis=0))
    rhovB=B*lambertw(wi8*rho*(1/B)).real
    sigmaJB=(rhovB-wi0*rho)*v2
    rhoiB=rho*wi8
    rhoiB[mobiles]=rhovB
    return rhoiB,sigmaJB
    
    
