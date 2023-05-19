import numpy as np

def MEOS(t,rhoiII,sigmaJ,tint,THFaktor,EJ,tauJ,exponent,nz,L0,mobiles,rho0,drhodtNF,Mi):
    nJ=len(EJ)
    nc=len(memoryview)
    rho=np.sum(rho0)
    etaJ=EJ*tauJ
    deltaz=L0/nz
    nz=20
    THcorr=np.ones((nc,nc))
    for i in range(nc):
        for j in range(nc):
            THcorr[i,j]=np.interp(t,tint,THFaktor[j,i,:])    
    #sigmaJ=np.zeros((nz+1,nJ))
    sigma=np.sum((sigmaJ.T*EJ).T,axis=1)
    dsigma=np.diff(sigma)
    v2=1/rho0[mobiles]
    rho2II=rhoiII[mobiles]
    X2II=rho2II/rho
    w2II=X2II/(X2II+1)
    WL=np.exp(-w2II*exponent)

    R=8.145
    T=298.15
    M2=Mi[mobiles]
    RV=R*T*1/M2*1/v2
    MDF=1/RV*dsigma/deltaz 
    drhodtNF+=np.diff(MDF)/deltaz
    etaWL=(etaJ@WL.T).T
    drho2dt_hist=np.sum(1/etaWL[-1,:].T*sigmaJ[-1,:].T*EJ**2/RV,axis=0)/(THcorr[mobiles,mobiles]/rho2II[-1]+np.sum(EJ)*v2/RV)
    drhodtNF[-1]=drho2dt_hist
    dsigmaJdt=np.zeros(nz+1,nJ)
    for i in range(nJ):
        dsigmaJdt[:,i]=-1/etaWL[:,i]*sigmaJ[:,i]*EJ[i]+drhodtNF*v2
    return drhodtNF,dsigmaJdt