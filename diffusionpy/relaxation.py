import numpy as np
from numba import njit

def relaxation_mode(ode,EJ,etaJ,exponent,T=298.15,Mv=18.015,vv=0.001):
    """alter the ode function in diffusionpy.Diffusion_MS, to also solve the relaxation of a component
    Args:
        ode (array_like): ode function which is modified by this function
        EJ (array_like): Elasticity Moduli
        etaJ (array_like): Damper constants
        exponent (array_like): Plasticization factors
        T (float, optional) : temperature /K
        Mv (float, optional) : solvent molar mass g/mol
        vv (float, optional) : solvent specific volume kg/m3
    Returns:
        array_like: new modified ode function with the same format as the input ode function
    """
    R=8.31448
    RV=R*T*1/(Mv/1000.)*1/vv
    nJ=len(EJ)
    #t,x,tint,THFaktor,mobiles,immobiles,D,allflux,wi0_relax,dmuext,wiB
    def relaxation_ode(*args):
        """solves the generalized Maxwell model for relaxation"""
        x=args[1]
        wi0=args[8]
        immobiles=args[5]
        nTH,nz_1=args[-2].shape
        nJ=len(EJ)
        wv=np.zeros((nTH,nz_1))
        for i in range(nTH):  wv[i,:]=x[(nz_1)*(i):(nz_1)*(1+i)]
        sigmaJ=np.zeros((nz_1,nJ))
        for J in range(nJ):   sigmaJ[:,J]=x[(nz_1)*(nTH+J):(nz_1)*(nTH+1+J)]
        wv=np.ascontiguousarray(wv)
        WL=np.prod(np.exp(-wv*np.expand_dims(exponent,-1)),0)
        etaWL=np.expand_dims(WL,1)*np.expand_dims(etaJ,0)
        dmuext=MDF(sigmaJ,EJ,RV)
        argsmod=list(args)
        argsmod[-2]=dmuext
        argsmod[1]=np.ascontiguousarray(wv.flatten())
        dwvdt_=ode(*argsmod)
        dwvdt=np.reshape(dwvdt_,(nTH,nz_1))
        omega=(np.sum(wi0[immobiles,:],axis=0)/(1-np.sum(wv,axis=0)))**-1
        rho=1/vv
        drhovdt=dwvdt*rho/omega-np.sum(dwvdt/np.sum(wi0[immobiles,:]),axis=0)*wv
        dsigmaJdt=stress(etaWL,EJ,sigmaJ,drhovdt,vv)
        fvec=np.hstack((dwvdt.flatten(),dsigmaJdt.flatten()))
        return fvec
    return relaxation_ode

@njit
def MDF(sigmaJ,EJ,RV):
    """the mechanical driving force for the stress gradient"""
    sigma=np.sum(sigmaJ*EJ,axis=1)
    nz_1=sigma.shape[0]
    dsigma=np.zeros_like(sigma)
    dsigma[1:]=np.diff(sigma)
    dsigma[0]=0
    dmuext=1/RV*np.expand_dims(dsigma,0)
    return dmuext

@njit
def stress(etaWL,EJ,sigmaJ,drhodtNF,vv):
    """calculate the change in the stresses of the maxwell elements"""
    nJ=len(EJ)
    nz_1,_=sigmaJ.shape
    dsigmaJdt=np.zeros((nz_1,nJ))
    for i in range(nJ):
        dsigmaJdt[:,i]=-1/etaWL[:,i]*sigmaJ[:,i]*EJ[i]+np.sum(drhodtNF*vv,axis=0)
    return dsigmaJdt

    
    
