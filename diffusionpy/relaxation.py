import numpy as np


def MEOS(t,rhov,sigmaJ,tint,THFaktor,EJ,tauJ,exponent,nz,L0,mobiles,rho0,drhodtNF,Mi):
    nJ=len(EJ)
    nc=len(Mi)
    rho=np.sum(rho0)
    etaJ=EJ*tauJ
    deltaz=L0/nz
    nz=20
    THcorr=np.ones((nc,nc))
    for i in range(nc):
        for j in range(nc):
            THcorr[i,j]=np.interp(t,tint,THFaktor[j,i,:])    
    #sigmaJ=np.zeros((nz+1,nJ))
    rhoi=np.zeros((nc,nz+1))
    rhoi[mobiles,:]=rhov
 
    sigma=np.sum(sigmaJ*EJ,axis=1)
    dsigma=np.diff(sigma)
    v2=1/rho0[mobiles]
    rhov=rhoi[mobiles,:]
    X2II=rhov/rho
    w2II=X2II/(X2II+1)
    WL=np.exp(-np.outer(w2II,exponent))

    R=8.145
    T=298.15
    M2=Mi[mobiles]
    RV=R*T*1/M2*1/v2
    MDF=np.hstack((0.,1/RV*dsigma/deltaz)) 
    drhodtNF[:-1]=drhodtNF[:-1]+np.diff(MDF)/deltaz
    etaWL=etaJ*WL
    b=np.sum(1/etaWL[-1,:].T*sigmaJ[-1,:].T*EJ**2/RV,axis=0)
    A=(THcorr[mobiles,mobiles]/rhov[:,-1]+np.sum(EJ)*v2/RV)
    drho2dt_hist=np.linalg.solve(A,b) if len(A)>1 else b/A
    drhodtNF[-1]=drho2dt_hist

    dsigmaJdt=np.zeros((nz+1,nJ))
    for i in range(nJ):
        dsigmaJdt[:,i]=-1/etaWL[:,i]*sigmaJ[:,i]*EJ[i]+np.sum(drhodtNF,axis=0)/rho
    return drhodtNF,dsigmaJdt
if __name__=="__main__":
    nt=500
    nJ=2
    nz=20
    nc=3
    t=40.
    tint=np.linspace(0,100,nt)
    rhoi=np.ones((nc,nz+1))
    sigmaJ=np.zeros((nz+1,nJ))
    THFaktor=np.zeros((nc,nc,nt))
    EJ=np.asarray([1.E9,2.E9])
    tauJ=np.asarray([1.,0.5])
    exponent=np.asarray([0.2,0.2])
    L0=1E-6
    mobiles=np.where(np.asarray([True,False,False]))[0]
    rho0=np.asarray([997.,1200.,1190.])
    drhodtNF=np.zeros(nz+1)
    Mi=np.asarray([18.15,25700.,350.])
    print(MEOS(t,rhoi,sigmaJ,tint,THFaktor,EJ,tauJ,exponent,nz,L0,mobiles,rho0,drhodtNF,Mi))
   


    
    
