import numpy as np
from numba import njit
import matplotlib.pyplot as plt
from .FEM_collocation import collocation,collocation_space

def relaxation_mode(rhovinit,ode,EJ, etaJ, exponent,M2,v2,tint,THFaktor,mobiles,immobiles,Mi,D,allflux,swelling,rho,wi0,dmuext,rhoiB,drhovdtB):
    """alter the ode function in diffusionpy.Diffusion_MS, to also solve the relaxation 

    Args:
        rhovinit (array_like): vector of the partial densities of the mobile components
        ode (array_like): ode fuinction which is modified by the function
        EJ (array_like): Elasticity Moduli
        etaJ (array_like): Damper constants
        exponent (array_like): Plasticization factors
        M2 (array_like): Molar masses of mobile compnents
        v2 (array_like): specific volumes of mobile compnents

    Returns:
        array_like: new modified ode function with the same format as the input ode function
    """
    R=8.31448
    T=298.15
    RV=R*T*1/(M2/1000.)*1/v2
    nJ=len(EJ)
    _,nz_1=rhovinit.shape
    def relaxation_ode(t,x,tint,THFaktor,mobiles,immobiles,Mi,D,allflux,swelling,rho,wi0,dmuext,rhoiB,drhovdtB):
        """solves the genralized Maxwell model for relaxation"""
        _,nz_1=dmuext.shape
        nTH=drhovdtB.shape[0]
        nJ=len(EJ)
        rhov=np.zeros((nTH,nz_1))
        for i in range(nTH):
            rhovtemp=x[(nz_1)*(i):(nz_1)*(1+i)]
            rhov[i,:]=rhovtemp
        sigmaJ=np.zeros((nz_1,nJ))
        for J in range(nJ):
            sigmaJtemp=x[(nz_1)*(nTH+J):(nz_1)*(nTH+1+J)]
            sigmaJ[:,J]=sigmaJtemp
        X2II=rhov/np.sum(rhoiB)
        w2II=X2II/(X2II+1)
        WL=np.prod(np.exp(-w2II*np.expand_dims(exponent,-1)),0)
        etaWL=np.expand_dims(WL,1)*np.expand_dims(etaJ,0)
        #rhoiB[0]=rhov[:,-1]
        rhov=np.ascontiguousarray(rhov)
        dmuext=MDF(sigmaJ,EJ,RV)
        #float64, array(float64, 1d, C), array(float64, 1d, C), pyobject, array(int64, 1d, C), array(int64, 1d, C), array(float64, 1d, C), array(float64, 2d, F), bool, bool, float64, array(float64, 1d, C), array(float64, 2d, C), array(float64, 2d, C), array(float64, 1d, C)
        drhovdt=ode(t,np.ascontiguousarray(rhov.flatten()),tint,THFaktor,mobiles,immobiles,Mi,D,allflux,swelling,rho,wi0,dmuext,rhoiB,drhovdtB)
        dsigmaJdt=stress(etaWL,EJ,sigmaJ,drhovdt,v2)
        # drhovdt[:,-1]=drhovdtB
        fvec=np.hstack((drhovdt.flatten(),dsigmaJdt.flatten()))
        return fvec
    
    sigmaJ0=np.zeros((nz_1,nJ))
    sigmaJB=np.zeros((nJ))
    sigmaJ0[-1,:]=sigmaJB
    xinit=np.hstack((rhovinit.flatten(),sigmaJ0.flatten()))
    return xinit,relaxation_ode

@njit
def MDF(sigmaJ,EJ,RV):
    """the mechanical driving force for the stress gradient"""
    sigma=np.sum(sigmaJ*EJ,axis=1)
    nz_1=sigma.shape[0]
    dsigma=collocation(sigma,nz_1,True)
    dsigma[0]=0
    dmuext=1/np.expand_dims(RV,1)*np.expand_dims(dsigma,0)
    return dmuext

@njit
def stress(etaWL,EJ,sigmaJ,drhodtNF,v2):
    """calculate the change in the stresses of the maxwell elements"""
    nJ=len(EJ)
    nz_1,_=sigmaJ.shape
    dsigmaJdt=np.zeros((nz_1,nJ))
    for i in range(nJ):
        dsigmaJdt[:,i]=-1/etaWL[:,i]*sigmaJ[:,i]*EJ[i]+np.sum(drhodtNF*np.expand_dims(v2,1),axis=0)
    return dsigmaJdt

# def boundary(rhov,etaWL,EJ,sigmaJ,RV,THFaktor,v2):
#     """calculate the mass fraction vector of the surface elment as afunction of time"""
#     b=np.sum(1/(etaWL[-1,:])*sigmaJ[-1,:]*EJ**2,axis=0)/RV
#     A=(THFaktor/rhov[:,-1]+np.sum(EJ)*v2/RV)
#     drhovdtB=np.linalg.solve(A,b) if len(A)>1 else b/A
#     return drhovdtB

# def initialboundary(RV,EJ,v2,THFaktors,rho,mobiles,wi0,wi8):
#     """calculate the boudnary of the surface concnetration and surface stresses as a function time"""
#     from scipy.special import lambertw
#     B=RV/EJ/v2*np.diag(np.mean(THFaktors,axis=0))
#     rhovB=B*lambertw(wi8*rho*(1/B)).real
#     sigmaJB=(rhovB-wi0*rho)*v2
#     rhoiB=rho*wi8
#     rhoiB[mobiles]=rhovB
#     return rhoiB,sigmaJB
    
    
