import numpy as np
import casadi as cs
import time
import copy

def SolveODEs(x,f,x0,opts,Time):
    dae={"x":x,"t":Time,"ode":f}
    fint=cs.integrator("integrator","idas",dae,opts)
    F=fint(x0=x0)
    return F["xf"] 

def averaging(a):
    return (a[1:]+a[:-1])/2

def Diffusion_MS(t,L,Dvec,w0,w8,Mi,volatile,full_output=False,Gammai=None):
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
    #Dref=1E-14
    #tau=(L**2)/Dref
    nc=len(w0)
    D=D_Matrix(Dvec,nc)
    nz=20
    zvec=np.linspace(0,L,nz+1)
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
    
    GammaiT=cs.MX.ones(nc)
    if Gammai is not None:
        Gammai_2=copy.deepcopy(Gammai.T)
        Gammai_2[Gammai_2<=0]=1.
        if np.any(Gammai_2!=Gammai.T):
           print("This is a warning")
        for i in range(nc):
            GammaiT[i]=cs.interpolant("Gammai_fun","linear",[t],Gammai_2[i,:])(Time)
            
    rhoi=cs.MX.ones((nc,nz+1))
    rhoi[np.where(volatile)[0],:]=rhov
    rhoi[np.where(~volatile)[0],:]=rhoiinit[np.where(~volatile)[0],:]
    rhoi=cs.fmax(rhoi,1E-8)


    
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
    

    opts = {"grid":t,"max_num_steps":1000,"regularity_check":True,"output_t0":True,"verbose":True}
    opts = {"grid":t,"regularity_check":True,"output_t0":True,"verbose":True}
    start=time.time_ns()
    x_sol=SolveODEs(x,ode,xinit,opts,Time)
    end=time.time_ns()
    print("------------- Diffusion modeling took "+str((end-start)/1E9)+" seconds ----------------")
    nt=len(t)
    wt=np.ones((nc,nt))
    wik=np.ones((nt,nc,nz+1))
    rhoik=rhoiinit
    for k in range(nt):
        rhovk=cs.reshape(x_sol[:,k],(nf,nz+1))
        rhoik[np.where(volatile)[0],:]=rhovk
        for i in range(nc): 
            wik[k,i,:]=rhoik[i,:]/cs.sum1(rhoik)
        wt[:,k]=cs.sum2(wik[k,:,:-1]/nz).full().flatten()
    return wt.T if not full_output else (wt.T,wik,zvec)


def BIJ_Matrix(D,wi,volatile):
    nc=wi.shape[0]
    nf=np.sum(volatile)
    allflux=nc==nf
    B=cs.MX.ones((nc,nc))
    comp=np.arange(0,nc)
    for i in range(nc):
        Din=wi[-1]/D[i,-1] if (i+1)!=nc else 0
        Dij=-wi[i]/D[i,:]
        BIJ= Dij if not allflux else Dij+Din
        Dii=cs.sum1((wi/D[i,:])[comp[comp!=i]])
        BII=Dii if not allflux else Dii+Din
        B[i,:]=BIJ
        B[i,i]=BII
    #alpha=1/np.min(D)
    #B+alpha*wi@wi.T

    return B[np.where(volatile)[0],np.where(volatile)[0]] if not allflux else B[:-1,:-1]  #8E-13

def D_Matrix(Dvec,nc):
    nd=(nc-1)*nc//2 
    if len(Dvec)!=nd: 
        raise Exception("Wrong number of diffusion coefficients. Provide array with "+str(nd)+" entries")
    else:
        D=np.zeros((nc,nc))
        D[np.triu_indices_from(D,k=1)]=Dvec
        D[np.tril_indices_from(D,k=-1)]=Dvec
    return D

def Diffusion1D(t,L0,Ds,ws0,ws8):
    """
    Method that computes the solvent sorption kinetics in a multi-component mixture 
    of non-volatile components
    Inputs
    ----------
    t  :    array_like  time                                /s
    L0 :    float       dry thickness                       /m
    Ds :    float       Solvent diffusion coefficient       /m^2/s
    ws0:    float       Solvent mass fraction at t=0        /-
    ws8:    float       Solvent mass fraction at t=infinity /-
    lnphi_fun
    Returns
    -------
    wst:    array_like  Solvent mass fraction at t          /-
    """
    nz=200
    deltaz=L0/nz
    Xs0=ws0/(1-ws0)
    Xs8=ws8/(1-ws8)
    Xs_ini=np.ones(nz+1)*Xs0
    Xs_ini[nz]=Xs8
    Ds0=Ds/deltaz**2
    Time=cs.MX.sym("Time")
    Xs=cs.MX.sym("Xs",nz+1)
    ws=Xs/(Xs+1)
    wsbar=averaging(ws)
    Xsbar= averaging(Xs)
    Ds=Ds0/(1-wsbar)
    lnfs= np.log(ws)
    TDF= cs.diff(lnfs)
    j_bnd=np.asarray([0.])
    js= cs.vertcat(j_bnd,Xsbar*Ds*TDF)
    dXsdt_bnd=np.asarray([0.])
    dXsdt=cs.vertcat(cs.diff(js),dXsdt_bnd)
    opts = {"grid":t,"max_num_steps":10000,"regularity_check":True,"output_t0":True}
    start=time.time_ns()
    Xs_sol=SolveODEs(Xs,dXsdt,Xs_ini,opts,Time).full()
    end=time.time_ns()
    print((end-start)/1E9)
    Xst=np.sum(Xs_sol[:-1,:],axis=0)/(nz)
    wst=Xst/(1+Xst)
    return wst

if __name__=="__main__":
    nt=100
    t=np.linspace(0,400,nt)*60
    nc=3
    nd=(nc-1)*nc//2 
    Dvec=np.asarray([1E-13,2.3E-13,1.7E-13])
    Dvec=np.asarray([1E-10,2.3E-10,3E-14])
    #np.fill_diagonal(D,np.ones(nc)*1E-30)
    L=0.001
    wi0=np.asarray([0.01,0.495,0.495])
    wi8=np.asarray([0.9,0.005,0.095])
    Mi=np.asarray([18.015,357.57,65000])

    volatile=np.asarray([True,True,False])
    wt,wtz,zvec=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,volatile,True)
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots()
    fig1,ax1=plt.subplots()
    ax.plot(t/60,wt[:,0])
    ax.plot(t/60,wt[:,1])
    ax.plot(t/60,wt[:,2])
    #[ax1.plot(zvec,wtz[i*10,0,:]) for i,val in enumerate(t[::10])]
    [ax1.plot(zvec,wtz[i*10,1,:]) for i,val in enumerate(t[::10])]
    #[ax1.plot(zvec,wtz[i*10,2,:]) for i,val in enumerate(t[::10])]
    ax.set_xlabel("t/min")
    ax.set_ylabel("wi/-")
    ax.legend(["w1","w2","w3"])
    plt.show()

