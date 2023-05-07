import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from numba import njit,config

@njit(['f8[:,:](f8, f8[:,:], f8[::1], f8[:,:,:], f8[::1], i8[::1], i8[::1], i8, f8[::1], f8[:,:], b1, b1, i8, f8, f8[:,:], f8[:], f8[:])'],cache=True)
def drhodt(t,rhov,tint,Gammai,taui,volatiles,nonvolatiles,nc,ri,D,allflux,swelling,nz,rho,rhoiinit,wi0,wi8):
    def averaging(a): return (a[1:]+a[:-1])/2.
    def BIJ_Matrix(D,wi,volatiles,allflux):
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
        return B[volatiles,:][:,volatiles] if not allflux else B[:-1,:-1]  #8E-13
    
    
    rhoi=np.zeros((nc,nz+1))
    rhoi[volatiles,:]=rhov
    rhoi[nonvolatiles,:]=rhoiinit[nonvolatiles,:]
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
            Gammacorr=np.interp(t,tint,Gammai[j,i,:])
            dlnwis+=Gammacorr*np.diff(np.log(np.fmax(wi[j,:],1E-6)))
        dlnwi[i,:]=dlnwis
        wibar[i,:]= averaging(wi[i,:])
        rhoibar[i,:]= averaging(rhoi[i,:])
    for z in range(nz):
        B=BIJ_Matrix(D,wibar[:,z],volatiles,allflux) 
        dmui=(dlnwi[:,z])
        omega=rho/np.sum(rhoibar[:,z])
        di=rho*wibar[:,z]*dmui/ri*omega if swelling else rhoibar[:,z]*dmui/ri
        if not allflux:
            ji[volatiles,z]=np.linalg.solve(B,di[volatiles])
        else:
            ji[:-1,z]=np.linalg.solve(B,di[:-1])

    for i in range(nc):
        dji=np.diff(np.hstack((np.asarray([0.]),ji[i,:])))
        djib=np.hstack((dji,np.asarray([0.]))) 
        drhoidt[i,:]=djib
    if taui[0]!=0.: drhoidt[volatiles,-1]=1./taui*np.exp(-t/taui)*(wi8[volatiles]-wi0[volatiles])*rho
    drhovdt=drhoidt[volatiles,:]
    return  drhovdt

def Diffusion_MS(t,L,Dvec,wi0,wi8,Mw,volatile,full_output=False,Gammai=None,swelling=False,taui=None):
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
    Mw:    array_like   Molar mass of components nc         /g/mol
    Gammai: array_like estimate for DlnaiDlnx at t          /-
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
    nf=np.sum(volatile)
    nt=len(t)
    rho=1200.
    refsegment=np.argmin(Mw)
    ri= Mw/Mw[refsegment]
    allflux=nc==nf
    nTH=nf if not allflux else nc-1
    volatiles=np.where(volatile)[0] if not allflux else np.arange(0,nc-1,dtype=np.int64)
    nonvolatiles=np.where(~volatile)[0] if not allflux else np.asarray([-1],dtype=np.int64)
    
    #initial conditions
    rhoiinit=np.zeros((nc,nz+1))
    for z in range(nz+1):
        rhoiinit[:,z]=wi0*rho
    rhoiinit[:,-1]=rho*wi8 if taui is None else rho*wi0
    if taui is None: taui=np.zeros_like(volatiles,dtype=float)
    rhoiinit[nonvolatiles,-1]=rho*wi8[nonvolatiles]
    rhovinit=rhoiinit[volatiles,:]
    xinit=np.reshape(rhovinit,np.multiply(*rhovinit.shape))
    GammaiT=np.asarray([np.eye(nc)]*nt).T
    if Gammai is not None:
        
        Gammai=Gammai.reshape((nc,nc,nt))
        wi0_nonvolatiles=wi0[nonvolatiles]/np.sum(wi0[nonvolatiles])
        THij = Gammai[volatiles,:,:][:,volatiles,:]
        massbalancecorrection=np.stack([np.sum((Gammai[volatiles,:,:][:,nonvolatiles,:])*wi0_nonvolatiles.reshape((1,nc-nTH,1)),axis=1)]*nTH)
        for i in range(nTH):
            for j in range(nTH):
                    GammaiT[volatiles[i],volatiles[j],:]= THij[i,j,:]-massbalancecorrection[j,i,:]

    def wrapdrhodt(t,rhov,tint,GammaiT,taui,volatiles,nonvolatiles,nc,ri,D,allflux,swelling,nz,rho,rhoiinit,wi0,wi8):
        rhov=np.ascontiguousarray(np.reshape(rhov,(nTH,nz+1)))
        drhovdt=drhodt(t,rhov,tint,np.ascontiguousarray(GammaiT),taui,volatiles,nonvolatiles,nc,ri,D,allflux,swelling,nz,rho,rhoiinit,wi0,wi8)
        return np.reshape(drhovdt,np.multiply(*drhovdt.shape))
    
    print("------------- Start diffusion modeling ----------------")
    start=time.time_ns()
    sol=solve_ivp(wrapdrhodt,(t[0],t[-1]),xinit,args=(t,GammaiT,taui,volatiles,nonvolatiles,nc,ri,D,allflux,swelling,nz,rho,rhoiinit,wi0,wi8),method="Radau",t_eval=t)
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
        rhoik[k,volatiles,:]=rhovk
        for i in range(nc):
            wik[k,i,:]=rhoik[k,i,:]/np.sum(rhoik[k,:,:],axis=0)
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

def Diffusion_MS_iter(t,L,Dvec,wi0,wi8,Mw,volatile,full_output=False,swelling=False,taui=None,par={}):
    Gammaiave=np.stack([dlnai_dlnxi(T,vpures,wi8*0.5+wi0*0.5,**par)]*nt).T
    wt_old=Diffusion_MS(t,L,Dvec,wi0,wi8,Mw,volatile,Gammai=Gammaiave,swelling=swelling,taui=taui)
    def wt_obj(wt_old):
        wt=wt_old.reshape((nt,nc))
        Gammai=np.asarray([dlnai_dlnxi(T,vpures,wt[i,:],**par) for i in range(nt)]).T
        return (wt-Diffusion_MS(t,L,Dvec,wi0,wi8,Mw,volatile,Gammai=Gammai,swelling=swelling,taui=taui)).flatten()
    return root(wt_obj,wt_old.flatten(),method="df-sane")["x"].reshape((nt,nc))

if __name__=="__main__":
    import matplotlib.pyplot as plt
    nt=500
    t=np.linspace(0,429.8518201,nt)
    Dvec=np.asarray([1E-6,2.3E-11,1.7E-11])
    nc=3
    L=2E-5
    wi0=np.asarray([0.333333333,0.333333333,0.333333333])
    wi8=np.asarray([0.00001,0.127606346,0.872393654])
    Mw=np.asarray([32.042,92.142,90000.])
    volatile=np.asarray([True,True,False])
    wt=Diffusion_MS(t,L,Dvec,wi0,wi8,Mw,volatile)
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
    "Mw" : np.asarray([32.042,  92.142, 90000.]),
    "kij": kij}
    vpures=vpure(p,T,**par)
    wt=Diffusion_MS_iter(t,L,Dvec,wi0,wi8,Mw,volatile,full_output=False,swelling=False,taui=None,par=par)
    plt.plot(t,wt[:,0])
    plt.plot(t,wt[:,1])
    plt.plot(t,wt[:,2])
    plt.show()