import numpy as np
from numba import njit
from .FEM_collocation import collocation,collocation_space

def liquidseperation_mode(wvinit,ode,kappaii,tint,THFaktor,mobiles,immobiles,Mi,D,allflux,wi0,dmuext,wiB):
    def liquidseperatiopn_ode(t,x,tint,THFaktor,mobiles,immobiles,Mi,D,allflux,wi0,dmuext,wiB):
        """solves the diffusion model coupled with the density gradient theory"""
        nTH,nz_1=dmuext.shape
        wv=np.zeros((nTH,nz_1))
        for i in range(nTH):
            wvtemp=x[(nz_1)*(i):(nz_1)*(1+i)]
            wv[i,:]=wvtemp
        for j in mobiles:
            wv[j,-1]=np.interp(t,tint,wiB[:,j])
        wv=np.ascontiguousarray(wv)
        if not np.allclose(THFaktor[0,0,:,:],np.eye(nTH)): dmuext=DGT(wv,kappaii) 
        return ode(t,np.ascontiguousarray(wv.flatten()),tint,THFaktor,mobiles,immobiles,Mi,D,allflux,wi0,dmuext,wiB)
    return wvinit.flatten(),liquidseperatiopn_ode

@njit
def DGT(wv,kappaii):
    """the mechanical driving force for the stress gradient"""
    nTH,nz_1=wv.shape
    dwv=np.zeros((nTH,nz_1))
    for j in range(nTH):
        dwv[j,:]=collocation(wv[j,:],nz_1,True)
    dwv[:,0]=0
    dlnwv=dwv/wv
    ddwv=np.zeros((nTH,nz_1))
    for j in range(nTH):
        ddwv[j,:]=collocation(dwv[j,:],nz_1,False)
        ddwv[j,-1]=0
    dddwv=np.zeros((nTH,nz_1))

    for j in range(nTH):
            dddwv[j,:]=collocation(ddwv[j,:],nz_1,True)
    dmuext=dddwv*np.expand_dims(kappaii,1)
    return dmuext

