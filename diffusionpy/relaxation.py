import numpy as np
from numba import njit

def MEOS_mode(ode,EJ, etaJ, exponent,RV,v2,mobiles):
    def MEOS_ode(t,x,THFaktor,dmuext,rhoiB,drhovdtB):
        nc,nz=dmuext.shape
        nTH=drhovdtB.shape[0]
        nJ=len(EJ)
        rhov=np.zeros((nTH,nz+1))
        for i in range(nTH):
            rhovtemp=x[(nz+1)*(i):(nz+1)*(1+i)]
            rhov[i,:]=rhovtemp
        sigmaJ=np.zeros((nz+1,nJ))
        for J in range(nJ):
            sigmaJtemp=x[(nz+1)*(nTH+J):(nz+1)*(nTH+1+J)]*np.atleast_1d(EJ)[J]
            sigmaJ[:,J]=sigmaJtemp
        rhov=np.ascontiguousarray(rhov)
        dmuext=MDF(sigmaJ,EJ,RV)
        drhovdt=ode(t,rhov,THFaktor,dmuext,rhoiB,drhovdtB)
        dsigmaJdt,drhovdtB=stress_and_boundary(rhov,drhovdt,sigmaJ,THFaktor,EJ,etaJ,exponent,RV,v2,mobiles)
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


def stress_and_boundary(rhov,drhodtNF,sigmaJ,THFaktor,EJ,etaJ,exponent,RV,v2,mobiles):
    nJ=len(EJ)
    rho=1200
    nz_1,_=sigmaJ.shape
    X2II=rhov/rho
    w2II=X2II/(X2II+1)
    WL=np.prod(np.exp(-w2II*np.expand_dims(exponent,-1)),0)
    etaWL=np.expand_dims(WL,1)*np.expand_dims(etaJ,0)
    b=np.sum(1/(etaWL[-1,:])*sigmaJ[-1,:]*EJ**2,axis=0)/RV
    #b=1/etaWL[-1,:].T*sigmaJ[-1,:].T*EJ**2/RV
    A=(THFaktor[np.ix_(mobiles,mobiles)]/rhov[:,-1]+np.sum(EJ)*v2/RV)
    drhovdtB=np.linalg.solve(A,b) if len(A)>1 else b/A
    dsigmaJdt=np.zeros((nz_1,nJ))
    for i in range(nJ):
        dsigmaJdt[:,i]=-1/etaWL[:,i]*sigmaJ[:,i]*EJ[i]+np.sum(drhodtNF,axis=0)/rho
    return dsigmaJdt,drhovdtB



    
    

  


    
    
