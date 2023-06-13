import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from scipy.interpolate import interp1d
from numba import njit,config
from .PyCSAFT_nue import dlnai_dlnxi,vpure,lngi
import time
config.DISABLE_JIT = True

@njit(['f8[:,:](f8, f8[:,::1], f8[:,::1], i8[::1], i8[::1], f8[::1], f8[:,:], b1, b1, f8, f8[::1],f8[:,:],f8[::1],f8[::1])'],cache=True)
def drhodt(t,rhov,THFaktor,mobiles,immobiles,Mi,D,allflux,swelling,rho,wi0,dmuext,rhoiB,drhovdtB):
    def averaging(a): return (a[1:,:]+a[:-1,:])/2.
    def np_linalg_solve(A,b):
        ret = np.empty_like(b)
        for i in range(b.shape[1]):
            ret[:, i] = np.linalg.solve(A[:,:,i], b[:,i])
        return ret
    def BIJ_Matrix(D,wi,mobiles,allflux):
        nc,nz=wi.shape
        B=np.zeros((nc,nc,nz))
        for i in range(nc):
            Din=wi[-1,:]/D[i,-1] if (i+1)!=nc else np.zeros(nz)
            Dii=np.zeros(nz)
            for j in range(nc):
                if j!=i:
                    Dij=-wi[i,:]/D[i,j] 
                    if allflux: Dij+=Din
                    B[i,j,:]=Dij
                    Dii+=wi[j,:]/D[i,j]
            if allflux: Dii+=Din
            B[i,i,:]=Dii
        return B[mobiles,:,:][:,mobiles,:] if not allflux else B[:-1,:-1,:]
    refsegment=np.argmin(Mi)
    ri= Mi[mobiles]/Mi[refsegment]
    nc=len(Mi)
    nTH,nz_1=rhov.shape
    rhoi=np.zeros((nc,nz_1))
    rhoi[mobiles,:]=rhov
    rhoi[immobiles,:]=rho*wi0[immobiles]
    rhoi[immobiles,-1]=rhoiB[immobiles]
    if not np.any(drhovdtB): rhoi[mobiles,-1]=rhoiB[mobiles]
    wi=rhoi/np.sum(rhoi,axis=0)
    wibar = averaging(wi.T).T
    rhoibar= averaging(rhoi.T).T
    dlnai=THFaktor@np.diff(np.log(np.fmax(wi[mobiles,:],1E-8)))
    B=BIJ_Matrix(D,wibar,mobiles,allflux) 
    dmui=dlnai+dmuext
    omega=rho/np.sum(rhoibar,axis=0)
    di=rho*wibar[mobiles,:]*dmui/np.atleast_2d(ri).T*omega if swelling else rhoibar[mobiles,:]*dmui/np.atleast_2d(ri).T
    ji=np_linalg_solve(B,di)
    jiB=np.zeros((nTH,1))
    dji=np.diff(np.hstack((jiB,ji)))
    drhovdt=np.hstack((dji,np.atleast_2d(drhovdtB).T))
    return  drhovdt


def Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,full_output=False,dlnai_dlnwi=None,swelling=False,**kwargs):
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
    rho=1200. if "rho0i" not in kwargs else np.sum(kwargs["rho0i"]*wi0)
    if "rho0i" not in kwargs  : rho0i=rho*np.ones(nc)
    allflux=nc==nf
    nTH=nf if not allflux else nc-1
    mobiles=np.where(mobile)[0] if not allflux else np.arange(0,nc-1,dtype=np.int64)
    immobiles=np.where(~mobile)[0] if not allflux else np.asarray([-1],dtype=np.int64)
    wi0_immobiles=wi0[immobiles]/np.sum(wi0[immobiles])
    rhoiinit=(rho*wi0*np.ones((nz+1,nc))).T
    rhovinit=rhoiinit[mobiles,:]
    #Construct TH Factor
    THFaktor=np.eye(nTH)
    if dlnai_dlnwi is not None:
        slc1=np.ix_(np.asarray(range(nt)),mobiles, mobiles)
        slc2=np.ix_(np.asarray(range(nt)),immobiles, mobiles)
        massbalancecorrection=np.sum(dlnai_dlnwi[slc2]*wi0_immobiles,axis=1)
        #THFaktor[slc1]=dlnai_dlnwi[slc1]-massbalancecorrection[:,None,:]
        THFaktor=dlnai_dlnwi[slc1]-massbalancecorrection[:,None,:]
        THFaktor= interp1d(t,THFaktor,axis=0,bounds_error=False,fill_value="extrapolate")
    
    
    #Set up the PDEs boundary,initial conditions and additional driving forces.
    #____________________________________
    #def default_mode():
    xinit=rhovinit.flatten()
    dmuext=np.zeros((nTH,nz))
    rhoiB=wi8*rho
    drhovdtB=np.zeros(nTH)
    def ode(t,x,THFaktor,dmuext,rhoiB,drhovdtB):
        rhov=np.ascontiguousarray(np.reshape(x,(nTH,nz+1))) 
        return drhodt(t,rhov,THFaktor,mobiles,immobiles,Mi,D,allflux,swelling,rho,wi0,dmuext,rhoiB,drhovdtB)
    #override these variables to alter the PDE
    #_____________________________________
    # Time-dependant surface concentration 
    # if  taui is not None: 
    #     from .surface_activity import time_dep_surface,gassided
    #     #rhoiB=lambda t: time_dep_surface(t,wi0,wi8,mobiles,immobiles,taui,rho)
    #     drhovdtB=lambda tvar,rhov: gassided(tvar,rhov,wi0[mobiles],wi8[mobiles],taui,THFaktor,t)
    #     def rhoiB(tvar,rhov):
    #         rhoiBtemp=wi0*rho
    #         rhoiBtemp[immobiles]=(rho-np.sum(np.reshape(rhov,rhovinit.shape),axis=0))[-1]*wi0[immobiles]/np.sum(wi0[immobiles])
    #         return rhoiBtemp
    #_____________________________________
    # Mechanical equation of state (MEOS)      
    if "EJ" in kwargs or "etaJ" in kwargs or "exponent" in kwargs: 
        from .relaxation import MEOS_mode
        xinit,ode=MEOS_mode(rhovinit,ode,kwargs["EJ"],kwargs["etaJ"],kwargs["exponent"],Mi[mobiles],1/rho0i[mobiles])
    #_____________________________________
    # if "par" in kwargs:
    if "witB" in kwargs: 
        rhoiB=lambda tvar: kwargs['witB'](tvar)*rho

    def wrapode(t,x,ode,THFaktor,dmuext,rhoiB,drhovdtB):
        THFaktor=THFaktor(t)    if callable(THFaktor)   else THFaktor
        rhoiB=rhoiB(t)          if callable(rhoiB)      else rhoiB
        drhovdtB=drhovdtB(t,x)    if callable(drhovdtB)   else drhovdtB
        dxdt=ode(t,x,THFaktor,dmuext,rhoiB,drhovdtB)
        return dxdt.flatten()
    print("------------- Start diffusion modeling ----------------")
    start=time.time_ns()
    sol=solve_ivp(wrapode,(t[0],t[-1]),xinit,args=(ode,THFaktor,dmuext,rhoiB,drhovdtB),method="Radau",t_eval=t)
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
        rhovk=np.reshape(x_sol[:(nz+1)*nTH,k],(nTH,nz+1))
        rhoik[k,:,:]=rhoiinit
        rhoik[k,mobiles,:]=rhovk
        if callable(rhoiB): rhoik[k,:,-1]=rhoiB(t[k])
        for i in range(nc):
            wik[k,i,:]=rhoik[k,i,:]/np.sum(rhoik[k,:,:],axis=0)
        # wik[k,mobiles,-1]=(wi8[mobiles]+(wi0[mobiles]-wi8[mobiles])*np.exp(-t[k]/taui)) if taui is not None else wi8[mobiles]
        # wik[k,mobiles,-1]=wi8[mobiles]
        wik[k,immobiles,-1]=(1-np.sum(wik[k,mobiles,-1],axis=0))*wi0_immobiles
        rhok[k]=np.sum(np.sum(rhoik[k,:,:-1]/nz,axis=1),axis=0)
        wt[:,k]=np.sum(wik[k,:,:-1]/nz,axis=1)
    Lt=rhok/rho*L
    # plt.plot(t,rhoik[:,1,-1])
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

def Diffusion_MS_iter(t,L,Dvec,wi0,wi8,Mi,mobile,full_output=False,swelling=False,T=298.15,par={},**kwargs):
    nt=len(t)
    nc=len(wi0)
    dlnai_dlnwi=np.stack([dlnai_dlnxi(T,wi8*0.5+wi0*0.5,**par)]*nt)
    wt_old=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,dlnai_dlnwi=dlnai_dlnwi,swelling=swelling,**kwargs)
    def wt_obj(wt_old):
        wt=wt_old.reshape((nt,nc))
        dlnai_dlnwi=np.asarray([dlnai_dlnxi(T,wt[i,:],**par) for i in range(nt)])
        return (wt-Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,dlnai_dlnwi=dlnai_dlnwi,swelling=swelling,**kwargs)).flatten()
    wtopt=root(wt_obj,wt_old.flatten(),method="df-sane",tol=1E-3)["x"].reshape((nt,nc))

    if not full_output:
        return wtopt
    else: 
        dlnai_dlnwiopt=np.asarray([dlnai_dlnxi(T,wtopt[i,:],**par) for i,val in enumerate(wtopt[:,0])]).T
        return Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,full_output=True,dlnai_dlnwi=dlnai_dlnwiopt,swelling=swelling)
    
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


    from .PyCSAFT_nue import dlnai_dlnxi,vpure,lngi
    T=298.15
    p=1E5
    Dvec=np.asarray([1E-6,2.3E-10,1.7E-10])
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
    dlnai_dlnwi=np.stack([dlnai_dlnxi(T,wi8*0.5+wi0*0.5,**par)]*nt)
    wt=Diffusion_MS_iter(t,L,Dvec,wi0,wi8,Mi,mobile,T=T,par=par)
    plt.plot(t,wt[:,0])
    plt.plot(t,wt[:,1])
    plt.plot(t,wt[:,2])

    taui=np.asarray([1,70])
    from .surface_activity import time_dep_surface2
    witB=time_dep_surface2(t,wi0,wi8,mobile,taui,par=None)
    wt=Diffusion_MS_iter(t,L,Dvec,wi0,wi8,Mi,mobile,T=T,par=par,witB=witB)
    # wt=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,dlnai_dlnwi=dlnai_dlnwi,taui=np.asarray([1,70]))
    # wt=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,dlnai_dlnwi=dlnai_dlnwi,taui=np.asarray([1,70]),par=par)
    plt.plot(t,wt[:,0])
    plt.plot(t,wt[:,1])
    plt.plot(t,wt[:,2])

    witB=time_dep_surface2(t,wi0,wi8,mobile,taui,par=par)
    wt=Diffusion_MS_iter(t,L,Dvec,wi0,wi8,Mi,mobile,T=T,par=par,witB=witB)
    # wt=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,dlnai_dlnwi=dlnai_dlnwi,taui=np.asarray([1,70]))
    # wt=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,dlnai_dlnwi=dlnai_dlnwi,taui=np.asarray([1,70]),par=par)
    plt.plot(t,wt[:,0])
    plt.plot(t,wt[:,1])
    plt.plot(t,wt[:,2])
    plt.savefig('filename.png', format='png',  transparent=True)
    plt.show()
    # wt=Diffusion_MS_iter(t,L,Dvec,wi0,wi8,Mi,mobile,T=T,par=par)
    # plt.plot(t,wt[:,0])
    # plt.plot(t,wt[:,1])
    # plt.plot(t,wt[:,2])

    # EJ=np.asarray([1E10])
    # etaJ=np.asarray([1E10])
    # exponent=np.asarray([0.,0.])
    # dlnai_dlnwi=np.stack([dlnai_dlnxi(T,wi8*0.5+wi0*0.5,**par)]*nt)
    # wt=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,dlnai_dlnwi=dlnai_dlnwi,EJ=EJ,etaJ=etaJ,exponent=exponent)
    # wt=Diffusion_MS_iter(t,L,Dvec,wi0,wi8,Mi,mobile,T=T,par=par,EJ=EJ,etaJ=etaJ,exponent=exponent)
    # plt.plot(t,wt[:,0])
    # plt.plot(t,wt[:,1])
    # plt.plot(t,wt[:,2])
    # plt.show()
    # t=np.linspace(0.,4000.*60.,nt)
    # Mi=np.asarray([18.015,65000])
    # EJ=np.asarray([1E10])
    # etaJ=np.asarray([1E9])
    # exponent=np.asarray([0.,0.])
    # mobile=np.asarray([True,False])
    # Dvec=np.asarray([5E-14])
    # wi0=np.asarray([1E-4,1-1E-4])
    # wi8=np.asarray([0.1,0.9])
    # #ai0=np.asarray([0.001,1])
    # #ai8=np.asarray([0.25,1])
    # a10=np.asarray([0.001])
    # a18=np.asarray([0.25])
    # L=1E-5
    # dlna1_dlnw1=(np.log(a10)-np.log(a18))/(np.log(wi0[0])-np.log(wi8[0]))
    # #dlnai_dlnwi=np.subtract.outer(np.log(ai0),np.log(ai8))/np.subtract.outer(np.log(wi0),np.log(wi8))
    # dlnai_dlnwi=np.asarray([[[dlna1_dlnw1,0.],[0.,1.]]]*len(t))
    # rho0i=np.asarray([997.,1180])
    # etavec=np.asarray([1E9,1E12,1E14,1E16])
    # for etaJ in etavec:
    #     wt=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,dlnai_dlnwi=dlnai_dlnwi,EJ=EJ,etaJ=etaJ,exponent=exponent,rho0i=rho0i)
    #     plt.plot((t/60)**(1/2),wt[:,0])
    # # plt.plot(t,wt[:,2])
    # plt.show()

