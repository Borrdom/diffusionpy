import xloil as xlo
from Stefan_Maxwell_segmental import Diffusion_MS,Diffusion1D,D_Matrix,averaging,BIJ_Matrix,SolveODEs
import numpy as np
import casadi as cs
import time

@xlo.func
def Diffusion_MS_xloil(t:xlo.Array(float,dims=1),L:float,Dvec:xlo.Array(float,dims=1),w0:xlo.Array(float,dims=1),w8:xlo.Array(float,dims=1),Mi:xlo.Array(float,dims=1),
volatile:xlo.Array(bool,dims=1),full_output:bool,Gammai:xlo.Array(float,dims=2)):
    """
    Method that computes the multi-component diffusion kinetics 
    Inputs
    ----------
    t  :    array_like  time                                /s
    L  :    float       dry thickness                       /m
    Dvec:   array_like  Vector of diffusion coefficients. 
    The length of this vector is nd=(nc-1)*n/2, where nc 
    is the number of components                             /m^2/s
    w0:    array_like   Mass fractions at t=0               /-
    w8:    array_like   Mass fraction at t=infinity         /-
    Mi:    array_like   Molar mass of components nc         /g/mol
    lngammai: array_like estimate for lngammai at t         /t
    Returns
    -------
    wt:    array_like   Matrix of mass fractions at t       /-
    """
    #volatile=volatile.flatten()
    #w0=w0.flatten()
    #w8=w0.flatten()
    #Mi=Mi.flatten()
    #t=t.flatten()
    #Dvec=Dvec.flatten()
    nc=len(w0)
    D=D_Matrix(Dvec,nc)
    nz=200
    nf=np.sum(volatile)
    nt=len(t)
    rho=1200
    refsegment=np.argmin(Mi)
    #initial
    rhoiinit=cs.DM.zeros((nc,nz+1))
    for z in range(nz+1):
        rhoiinit[:,z]=w0*rho
    rhoiinit[:,-1]=w8*rho
    rhovinit=rhoiinit[np.where(volatile)[0],:]

    #decision variable vector
    rhov=cs.MX.sym("rhov",(nf,nz+1))
    Time=cs.MX.sym("Time")
    #lngammai=np.zeros((nc,nt))
    #lngammai=lngammai.T
    #lngamma0i=lngammai[:,0]
    #lngamma8i=lngammai[:,-1]
    #lngammaiT=cs.MX.zeros(nc)
    Gammai=Gammai.T
    GammaiT=cs.MX.zeros(nc)
    for i in range(nc):
        GammaiT[i]=cs.interpolant("Gammai_fun","linear",[t],Gammai[i,:])(Time)
        #GammaiT[i]=cs.gradient(lngammaiT[i],Time)
    rhoi=cs.MX.ones((nc,nz+1))
    rhoi[np.where(volatile)[0],:]=rhov
    rhoi[np.where(~volatile)[0],:]=rhoiinit[np.where(~volatile)[0],:]



    
    ji=cs.MX.zeros((nc,nz))
    dlnwi=cs.MX.ones((nc,nz))
    wibar=cs.MX.ones((nc,nz))
    rhoibar=cs.MX.ones((nc,nz))
    drhoidt=cs.MX.ones((nc,nz+1))
    wi=cs.MX.zeros((nc,nz+1))
    dz=L/nz
    ri= Mi/Mi[refsegment]

    for i in range(nc):
        wi[i,:]=rhoi[i,:]/cs.sum1(rhoi)
        dlnwi[i,:]= cs.diff(cs.log(wi[i,:]))*GammaiT[i]
        wibar[i,:]= averaging(wi[i,:])
        rhoibar[i,:]= averaging(rhoi[i,:])

    
    for z in range(nz):
        B=BIJ_Matrix(D,wibar[:,z],volatile)
        allflux=nc==nf
        Binv=cs.inv(B)
        dmui=(dlnwi[:,z])/dz
        di=rhoibar[:,z]*dmui/ri
        if not allflux:
            ji[np.where(volatile)[0],z]=cs.sum1(di[np.where(volatile)[0]]*Binv)
        else:
            ji[:-1,z]=cs.sum1(di[:-1]*Binv)
    
    ji[-1,:]=-cs.sum1(ji[:-1,:]) if allflux else ji[-1,:]
    for i in range(nc):
        dji=cs.diff(cs.horzcat(0,ji[i,:]))
        drhoidt[i,:]=cs.horzcat(dji/dz,0)
    drhovdt=drhoidt[np.where(volatile)[0],:]
    ode=cs.reshape(drhovdt,(np.multiply(*drhovdt.shape),1))
    x=cs.reshape(rhov,(np.multiply(*rhov.shape),1))
    xinit=cs.reshape(rhovinit,(np.multiply(*rhovinit.shape),1))

    #np.cond(A) 10**16 
    

    opts = {"grid":t,"max_num_steps":10000,"regularity_check":True,"output_t0":True}
    start=time.time_ns()
    x_sol=SolveODEs(x,ode,xinit,opts,Time)
    end=time.time_ns()
    print("------------- Diffusion modeling took "+str((end-start)/1E9)+" seconds ----------------")
    nt=len(t)
    wt=np.ones((nc,nt))
    wik=np.ones((nc,nz+1))
    rhoik=rhoiinit
    for k in range(nt):
        rhovk=cs.reshape(x_sol[:,k],(nf,nz+1))
        rhoik[np.where(volatile)[0],:]=rhovk
        for i in range(nc): 
            wik[i,:]=rhoik[i,:]/cs.sum1(rhoik)
        wt[:,k]=cs.sum2(wik[:,:-1]/nz).full().flatten()
    return wt.T if not full_output else (wt,x_sol)

@xlo.func
def Diffusion1D_xloil(t:xlo.Array(float,dims=1),L0:float,Ds:float,ws0:float,ws8:float):
    return Diffusion1D(t,L0,Ds,ws0,ws8)
    
@xlo.func
def gradient(x:xlo.Array(float,dims=1),y:xlo.Array(float,dims=1)):
    return np.gradient(x,y)