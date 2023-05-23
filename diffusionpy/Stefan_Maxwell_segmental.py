import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from numba import njit,config
from .PyCSAFT_nue import dlnai_dlnxi,vpure,lngi
import time
#config.DISABLE_JIT = True

@njit(['f8[:,:](f8, f8[:,:], f8[::1], f8[:,:,:], f8[::1], i8[::1], i8[::1], i8, f8[::1], f8[:,:], b1, b1, i8, f8, f8[:], f8[:])'],cache=True)
def drhodt(t,rhov,tint,THFaktor,taui,mobiles,immobiles,nc,ri,D,allflux,swelling,nz,rho,wi0,wi8):
    def averaging(a): return (a[1:]+a[:-1])/2.
    def BIJ_Matrix(D,wi,mobiles,allflux):
        nc=wi.shape[0]
        B=np.zeros((nc,nc))
        for i in range(nc):
            Din=wi[-1]/D[i,-1] if (i+1)!=nc else 0
            Dii=0
            for j in range(nc):
                if j!=i:
                    Dij=-wi[i]/D[i,j] 
                    if allflux: Dij+=Din
                    B[i,j]=Dij
                    Dii+=wi[j]/D[i,j]
            if allflux: Dii+=Din
            B[i,i]=Dii
        return B[mobiles,:][:,mobiles] if not allflux else B[:-1,:-1]  #8E-13
    
    
    rhoi=np.zeros((nc,nz+1))
    rhoi[mobiles,:]=rhov
    rhoi[immobiles,:]=rho*wi0[immobiles]
    wi0_immobiles=wi0[immobiles]/np.sum(wi0[immobiles])
    if taui[0]!=0.: rhoi[mobiles,-1]=(wi8[mobiles]+(wi0[mobiles]-wi8[mobiles])*np.exp(-t/taui))*rho
    rhoi[immobiles,-1]=(rho-np.sum(rhoi[mobiles,-1],axis=0))*wi0_immobiles if taui[0]!=0. else rho*wi8[immobiles]
    ji=np.zeros((nc,nz))
    dlnwi=np.zeros((nc,nz))
    wibar=np.zeros((nc,nz))
    rhoibar=np.zeros((nc,nz))
    drhoidt=np.zeros((nc,nz+1))
    wi=np.zeros((nc,nz+1))
    for i in range(nc):
        wi[i,:]=rhoi[i,:]/np.sum(rhoi,axis=0)
    for i in range(nc):
        dlnwis=np.zeros(nz)
        for j in range(nc):
            THcorr=np.interp(t,tint,THFaktor[j,i,:])
            dlnwis+=THcorr*np.diff(np.log(np.fmax(wi[j,:],1E-8)))
        dlnwi[i,:]=dlnwis
        wibar[i,:]= averaging(wi[i,:])
        rhoibar[i,:]= averaging(rhoi[i,:])
    for z in range(nz):
        B=BIJ_Matrix(D,wibar[:,z],mobiles,allflux) 
        dmui=(dlnwi[:,z])
        omega=rho/np.sum(rhoibar[:,z])
        di=rho*wibar[:,z]*dmui/ri*omega if swelling else rhoibar[:,z]*dmui/ri
        if not allflux:
            ji[mobiles,z]=np.linalg.solve(B,di[mobiles])
        else:
            ji[:-1,z]=np.linalg.solve(B,di[:-1])

    for i in range(nc):
        dji=np.diff(np.hstack((np.asarray([0.]),ji[i,:])))
        djib=np.hstack((dji,np.asarray([0.]))) 
        drhoidt[i,:]=djib
    drhovdt=drhoidt[mobiles,:]
    return  drhovdt

def Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,full_output=False,dlnai_dlnwi=None,swelling=False,taui=None,rho0i=None,**kwargs):
    """
    Method that computes the multi-component diffusion kinetics 
    Inputs
    ----------
    t  :    array_like  time                                /s
    L  :    float       dry thickness                       /m
    Dvec:   array_like  Vector of diffusion coefficients. 
    The length of this vector is nd=(nc-1)*n/2, where nc 
    is the number of components                             /m^2/s
    wi0:    array_like   Mass fractions at t=0               /-
    wi8:    array_like   Mass fraction at t=infinity         /-
    Mi:    array_like   Molar mass of components nc         /g/mol
    dlnai_dlnwi: array_like estimate for DlnaiDlnx at t          /-
    Returns
    -------
    wt:    array_like   Matrix of mass fractions at t       /-

    Full output
    -------
    wt:    array_like   Matrix of mass fractions at t       /-
    wtz:    array_like  Matrix of mass fractions at t,z     /-
    """
    nc=len(wi0)
    nz=20
    dz=L/nz
    D=D_Matrix(Dvec/dz**2,nc)
    zvec=np.linspace(0,L,nz+1)
    nf=np.sum(mobile)
    nt=len(t)
    rho=1200. if rho0i is None else np.sum(rho0i*wi0)
    refsegment=np.argmin(Mi)
    ri= Mi/Mi[refsegment]
    allflux=nc==nf
    nTH=nf if not allflux else nc-1
    mobiles=np.where(mobile)[0] if not allflux else np.arange(0,nc-1,dtype=np.int64)
    immobiles=np.where(~mobile)[0] if not allflux else np.asarray([-1],dtype=np.int64)
    
    #initial conditions
    rhoiinit=np.broadcast_to(rho*wi0,(nz+1,nc)).T.copy()
    wi0_immobiles=wi0[immobiles]/np.sum(wi0[immobiles])
    rhoiinit[:,-1]=rho*wi8 if taui is None else rho*wi0
    if taui is None: taui=np.zeros_like(mobiles,dtype=float)
    #rhoiinit[immobiles,-1]=rho*wi8[immobiles]
    rhovinit=rhoiinit[mobiles,:]
    xinit=np.reshape(rhovinit,np.multiply(*rhovinit.shape))
    THFaktor=np.asarray([np.eye(nc)]*nt).T
    if dlnai_dlnwi is not None:
        if len(dlnai_dlnwi.shape)<3: dlnai_dlnwi=dlnai_dlnwi.T
        dlnai_dlnwi=dlnai_dlnwi.reshape((nc,nc,nt))
        THij = dlnai_dlnwi[mobiles,:,:][:,mobiles,:]
        massbalancecorrection=np.stack([np.sum((dlnai_dlnwi[mobiles,:,:][:,immobiles,:])*wi0_immobiles.reshape((1,nc-nTH,1)),axis=1)]*nTH)
        for i in range(nTH):
            for j in range(nTH):
                    THFaktor[mobiles[i],mobiles[j],:]= THij[i,j,:]-massbalancecorrection[j,i,:]

    def wrapdrhodt(t,rhov,tint,THFaktor,taui,mobiles,immobiles,nc,ri,D,allflux,swelling,nz,rho,wi0,wi8):
        rhov=np.ascontiguousarray(np.reshape(rhov,(nTH,nz+1)))
        drhovdt=drhodt(t,rhov,tint,np.ascontiguousarray(THFaktor),taui,mobiles,immobiles,nc,ri,D,allflux,swelling,nz,rho,wi0,wi8)
        return np.reshape(drhovdt,np.multiply(*drhovdt.shape))

    if "EJ" or "tauJ" or "exponent" in kwargs: 
        EJ=kwargs["EJ"]
        exponent=kwargs["exponent"]
        tauJ=kwargs["tauJ"]
        nJ=len(EJ)
        def wrapdrhodt(t,xvec,tint,THFaktor,taui,mobiles,immobiles,nc,ri,D,allflux,swelling,nz,rho,wi0,wi8,EJ,tauJ,exponent):
            rhov=np.zeros((nTH,nz+1))
            for i in range(nTH):
                rhovtemp=xvec[(nz+1)*(i):(nz+1)*(1+i)]
                rhov[:,i]=rhovtemp
            sigmaJ=np.zeros((nz+1,nJ))
            for J in range(nJ):
                sigmaJtemp=xvec[(nz+1)*(nTH+J):(nz+1)*(nTH+1+J)]*np.atleast_1d(EJ)[J]
                sigmaJ[:,J]=sigmaJtemp
            rhov=np.ascontiguousarray(rhov)
            drhovdt=drhodt(t,rhov,tint,np.ascontiguousarray(THFaktor),taui,mobiles,immobiles,nc,ri,D,allflux,swelling,nz,rho,wi0,wi8)
            from .relaxation import MEOS
            if rho0i is None: rho0i=rho*np.ones(nc)
            drhovdt,dsigmaJdt=MEOS(t,rhov,sigmaJ,tint,THFaktor,EJ,tauJ,exponent,nz,L,mobiles,rho0i,drhovdt,Mi)
            dsigmaJdtvec=dsigmaJdt.flatten()
            #xvec=np.hstack(rhov.flatten(),sigmaJ.flatten())
            fvec=np.hstack(drhovdt.flatten(),dsigmaJdtvec.flatten())
            return fvec

    
    print("------------- Start diffusion modeling ----------------")
    start=time.time_ns()
    sol=solve_ivp(wrapdrhodt,(t[0],t[-1]),xinit,args=(t,THFaktor,taui,mobiles,immobiles,nc,ri,D,allflux,swelling,nz,rho,wi0,wi8),method="Radau",t_eval=t)
    end=time.time_ns()
    print("------------- Diffusion modeling took "+str((end-start)/1E9)+" seconds ----------------")
    if not sol["success"]: print("------------- Modeling failed the initial conditions are returned instead ----------------"); return wi0*np.ones((nc,nt)).T 
    x_sol=sol["y"] 
    nt=len(t)
    wt=np.zeros((nc,nt))
    wik=np.zeros((nt,nc,nz+1))
    rhoik=np.zeros((nt,nc,nz+1))
    rhok=np.zeros(nt)
    for k in range(nt):
        rhovk=np.reshape(x_sol[:,k],(nTH,nz+1))
        rhoik[k,:,:]=rhoiinit
        rhoik[k,mobiles,:]=rhovk
        for i in range(nc):
            wik[k,i,:]=rhoik[k,i,:]/np.sum(rhoik[k,:,:],axis=0)
        wik[k,mobiles,-1]=(wi8[mobiles]+(wi0[mobiles]-wi8[mobiles])*np.exp(-t[k]/taui))
        wik[k,immobiles,-1]=(1-np.sum(wik[k,mobiles,-1],axis=0))*wi0_immobiles
        rhok[k]=np.sum(np.sum(rhoik[k,:,:-1]/nz,axis=1),axis=0)
        wt[:,k]=np.sum(wik[k,:,:-1]/nz,axis=1)
    Lt=rhok/rho*L
    return wt.T if not full_output else (wt.T,wik,zvec,Lt)


def D_Matrix(Dvec,nc):
    nd=(nc-1)*nc//2 
    if len(Dvec)!=nd: 
        raise Exception("Wrong number of diffusion coefficients. Provide array with "+str(nd)+" entries")
    else:
        D=np.zeros((nc,nc))
        D[np.triu_indices_from(D,k=1)]=Dvec
        D[np.tril_indices_from(D,k=-1)]=Dvec
    return D

def Diffusion_MS_iter(t,L,Dvec,wi0,wi8,Mi,mobile,full_output=False,swelling=False,taui=None,rho0i=None,T=298.15,par={}):
    nt=len(t)
    nc=len(wi0)
    dlnai_dlnwi=np.stack([dlnai_dlnxi(T,wi8*0.5+wi0*0.5,**par)]*nt).T
    wt_old=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,dlnai_dlnwi=dlnai_dlnwi,swelling=swelling,taui=taui,rho0i=rho0i)
    def wt_obj(wt_old):
        wt=wt_old.reshape((nt,nc))
        dlnai_dlnwi=np.asarray([dlnai_dlnxi(T,wt[i,:],**par) for i in range(nt)]).T
        return (wt-Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,dlnai_dlnwi=dlnai_dlnwi,swelling=swelling,taui=taui,rho0i=rho0i)).flatten()
    wtopt=root(wt_obj,wt_old.flatten(),method="df-sane",tol=1E-4)["x"].reshape((nt,nc))

    if not full_output:
        return wtopt
    else: 
        dlnai_dlnwiopt=np.asarray([dlnai_dlnxi(T,wtopt[i,:],**par) for i,val in enumerate(wtopt[:,0])]).T
        return Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,full_output=True,dlnai_dlnwi=dlnai_dlnwiopt,swelling=swelling,taui=taui,rho0i=rho0i)
    
def convert(x,M,axis=0):
    y=np.empty(x.shape,M.dtype)
    np.copyto(y.T,M)
    return x*y/np.sum(x*y,axis=axis)

if __name__=="__main__":
    import matplotlib.pyplot as plt
    nt=500
    t=np.linspace(0,429.8518201,nt)
    Dvec=np.asarray([1E-6,2.3E-11,1.7E-11])
    nc=3
    L=2E-5
    wi0=np.asarray([0.333333333,0.333333333,0.333333333])
    wi8=np.asarray([0.00001,0.127606346,0.872393654])
    Mi=np.asarray([32.042,92.142,90000.])
    mobile=np.asarray([True,True,False])
    wt=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile)
    plt.plot(t,wt[:,0])
    plt.plot(t,wt[:,1])
    plt.plot(t,wt[:,2])
    from .PyCSAFT_nue import dlnai_dlnxi,vpure,lngi
    T=298.15
    p=1E5
    kij=D_Matrix(np.asarray([0.029,-0.05855362,0.027776682]),nc)
    par={"mi" :np.asarray([1.5255, 2.8149, 2889.9]),
    "ui" : np.asarray([188.9, 285.69, 204.65]),
    "si" : np.asarray([3.2300, 3.7169, 3.3972]),
    "kAi": np.asarray([0.035176, 0., 0.]),
    "eAi": np.asarray([2899.5, 0., 0.]),
    "NAi": np.asarray([1., 0., 1047.]),
    "Mi" : np.asarray([32.042,  92.142, 90000.]),
    "kij": kij}
    vpures=vpure(p,T,**par)
    par["vpure"]=vpures
    
    wt=Diffusion_MS_iter(t,L,Dvec,wi0,wi8,Mi,mobile,T=T,par=par)
    plt.plot(t,wt[:,0])
    plt.plot(t,wt[:,1])
    plt.plot(t,wt[:,2])
    plt.show()
