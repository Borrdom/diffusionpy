import numpy as np
import casadi as cs
import time
import copy
import matplotlib.pyplot as plt

def SolveODEs(x,f,x0,opts,Time):
    dae={"x":x,"t":Time,"ode":f}
    fint=cs.integrator("integrator","idas",dae,opts)
    F=fint(x0=x0)
    return F["xf"] 

def averaging(a):
    return (a[1:]+a[:-1])/2
def diff(a):
    #abar=averaging(a)
    #b=4/3*(a[1:]-a[:-1])
    #c=b[:1]-2/12*(abar[1:]-abar[:-1])
    #c=b[1:]+1/2*a[2:]-1/2*a[:-2]
    #return cs.horzcat(a[1]-abar[0],c)
    return (a[1:]-a[:-1])

def Diffusion_MS(t,L,Dvec,w0,w8,Mi,volatile,full_output=False,Gammai=None,swelling=False,taui=None):
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
    Gammai: array_like estimate for DlnaiDlnx at t          /-
    Returns
    -------
    wt:    array_like   Matrix of mass fractions at t       /-

    Full output
    -------
    wt:    array_like   Matrix of mass fractions at t       /-
    wtz:    array_like  Matrix of mass fractions at t,z     /-
    """

    #Dref=np.average(Dvec)
    #tau=(L**2)/Dref
    nc=len(w0)
    nz=20
    dz=L/nz
    D=D_Matrix(Dvec/dz**2,nc)
    zvec=np.linspace(0,L,nz+1)
    nf=np.sum(volatile)
    nt=len(t)
    rho=1200
    refsegment=np.argmin(Mi)
    allflux=nc==nf
    isvolatile=np.where(volatile)[0]
    notvolatile=np.where(~volatile)[0]
    #initial conditions
    rhoiinit=cs.DM.zeros((nc,nz+1))
    for z in range(nz+1):
        rhoiinit[:,z]=w0*rho
    #phi=w8[1]/w0[1]
    #rhoiinit[:,-1]=(1-phi)/phi*rho*w8
    #swelling=w0II is not None
    #omega8=(w0-w0II)/(w8-w0II) if swelling else 1
    #omega8=1
    rhoiinit[:,-1]=rho*w8 if taui is None else rho*w0
    rhoiinit[notvolatile,-1]=rho*w8[notvolatile]
    rhovinit=rhoiinit[isvolatile,:]

    #decision variable vector
    rhov=cs.MX.sym("rhov",(nf,nz+1))
    Time=cs.MX.sym("Time")
    lambi=1 if taui is None else 1-cs.exp(-Time/taui)
    dlambidt=0 if taui is None else 1/taui*cs.exp(-Time/taui)

    x=cs.reshape(rhov,(np.multiply(*rhov.shape),1))
    
    
    GammaiT=cs.MX.eye(nc)
    if Gammai is not None:
        nTH=nf if not allflux else nc-1
        Gammai=Gammai.reshape((nc,nc,nt))
        comp=np.arange(0,nc)
        volatiles=comp[isvolatile] if not allflux else comp[isvolatile][:-1]
        nonvolatiles=comp[notvolatile] if not allflux else comp[-1,None]
        w0_nonvolatiles=w0[nonvolatiles]/np.sum(w0[nonvolatiles])
        THij = Gammai[volatiles,:,:][:,volatiles,:]
        massbalancecorrection=np.stack([np.sum((Gammai[volatiles,:,:][:,nonvolatiles,:])*w0_nonvolatiles.reshape((1,nc-nTH,1)),axis=1)]*nTH)
        #THij -= massbalancecorrection
        for i in range(nTH):
            for j in range(nTH):
                    GammaiT[volatiles[i],volatiles[j]]=cs.interpolant("Gammai_fun","linear",[t],THij[i,j,:]-massbalancecorrection[j,i,:])(Time)
           
    # if Gammai is not None:
    #     Gammai_2=copy.deepcopy(Gammai.T)
    #     Gammai_2[Gammai_2<=0]=1.
    #     if np.any(Gammai_2!=Gammai.T):
    #        print("This is a warning")

    rhoi=cs.MX.ones((nc,nz+1))
    rhoi[isvolatile,:]=rhov
    rhoi[notvolatile,:]=rhoiinit[notvolatile,:]

    
    ji=cs.MX.zeros((nc,nz))
    dlnwi=cs.MX.ones((nc,nz))
    wibar=cs.MX.ones((nc,nz))
    rhoibar=cs.MX.ones((nc,nz))
    drhoidt=cs.MX.ones((nc,nz+1))
    wi=cs.MX.zeros((nc,nz+1))
    
    ri= Mi/Mi[refsegment]


    for i in range(nc):
        wi[i,:]=rhoi[i,:]/cs.sum1(rhoi)
        
    GammaiT=GammaiT.T
    for i in range(nc):
        #dlnwi[i,:]=cs.diff(cs.log(wi[i,:]))*GammaiT[i,i] # here it would be for every component
        dlnwi[i,:]=cs.sum1(cs.diff(cs.log(wi),1,1)*GammaiT[i,:].T)# here it would be for every component
        wibar[i,:]= averaging(wi[i,:])
        rhoibar[i,:]= averaging(rhoi[i,:])
    for z in range(nz):
        B=BIJ_Matrix(D,wibar[:,z],volatile) 
        Binv=cs.inv(B)
        dmui=(dlnwi[:,z])
        #di=rhoibar[:,z]*dmui/ri # in rhoibar steckt die transformierte riesige Dichte drin. Der Volumenexpansionsfaktor korriegiert diesen
        omega=rho/cs.sum1(rhoibar[:,z])
        di=rho*wibar[:,z]*dmui/ri*omega if swelling else rhoibar[:,z]*dmui/ri
        #di=rhoibar[:,z]*dmui/ri*omega**2
        
        if not allflux:
            ji[isvolatile,z]=di[isvolatile].T@Binv
        else:
            ji[:-1,z]=di[:-1].T@Binv

    for i in range(nc):
        dji=cs.diff(cs.horzcat(0,ji[i,:]))
        djib=cs.horzcat(dji,0) 
        drhoidt[i,:]=djib
    
    X0=w0[isvolatile]/(1-cs.sum1(w0[isvolatile])) if not allflux else w0
    X8=w8[isvolatile]/(1-cs.sum1(w8[isvolatile])) if not allflux else w8
    #drhoidt[isvolatile,-1]=dlambidt*(X8-X0)*rho
    drhoidt[isvolatile,-1]=dlambidt*(w8[isvolatile]-w0[isvolatile])*rho
    # need for boundary when swelling
    # drhoidtmean=cs.sum2(drhoidt[:,:-1])/nz
    # drhoidt[:,-1]=cs.sum1(drhoidtmean)*w8
    #is this it

    drhovdt=drhoidt[isvolatile,:]
    ode=cs.reshape(drhovdt,(np.multiply(*drhovdt.shape),1))
    
    xinit=cs.reshape(rhovinit,(np.multiply(*rhovinit.shape),1))


    opts = {"grid":t,"regularity_check":True,"output_t0":True,"verbose":True}
    start=time.time_ns()
    x_sol=SolveODEs(x,ode,xinit,opts,Time)
    end=time.time_ns()
    print("------------- Diffusion modeling took "+str((end-start)/1E9)+" seconds ----------------")
    nt=len(t)
    wt=np.zeros((nc,nt))
    wik=np.zeros((nt,nc,nz+1))
    rhoik=np.zeros((nt,nc,nz+1))
    rhok=np.zeros(nt)
    for k in range(nt):
        rhovk=cs.reshape(x_sol[:,k],(nf,nz+1))
        rhoik[k,:,:]=rhoiinit
        rhoik[k,isvolatile,:]=rhovk
        for i in range(nc):
            wik[k,i,:]=rhoik[k,i,:]/cs.sum1(rhoik[k,:,:]).full()
        rhok[k]=cs.sum1(cs.sum2(rhoik[k,:,:-1]/nz)).full()
        wt[:,k]=cs.sum2(wik[k,:,:-1]/nz).full().flatten()
    Lt=rhok/rho*L
    #[plt.plot(zvec,rhoik[k,1,:]) for k,val in enumerate(rhoik[:,1,0])]
    #plt.show()
    return wt.T if not full_output else (wt.T,wik,zvec,Lt)


def BIJ_Matrix(D,wi,volatile):
    nc=wi.shape[0]
    nf=np.sum(volatile)
    allflux=nc==nf
    B=cs.MX.ones((nc,nc))
    comp=np.arange(0,nc)
    isvolatile=np.where(volatile)[0]
    notvolatile=np.where(~volatile)[0]
    for i in range(nc):
        Din=wi[-1]/D[i,-1] if (i+1)!=nc else 0
        Dij=-wi[i]/D[i,:]
        BIJ= Dij if not allflux else Dij+Din
        Dii=cs.sum1((wi/D[i,:])[comp[comp!=i]])
        BII=Dii if not allflux else Dii+Din
        B[i,:]=BIJ
        B[i,i]=BII
    #alpha=1/np.min(D) #implement an augmented transport matrix
    #B+alpha*wi@wi.T

    return B[isvolatile,isvolatile] if not allflux else B[:-1,:-1]  #8E-13

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
    t=np.linspace(0,500,nt)*60
    # nc=2
    # nd=(nc-1)*nc//2

    # Dvec=np.asarray([1E-13])
    # #Dvec=np.asarray([1E-10,2.3E-10,3E-14])
    # #np.fill_diagonal(D,np.ones(nc)*1E-30)
    # L=0.0002
    # wi0=np.asarray([0.01,0.99])
    # wi8=np.asarray([0.3,0.7])
    # Mi=np.asarray([18.015,25700])
    # volatile=np.asarray([True,False])
    # wt=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,volatile)
    # from .PyCSAFT_nue import DlnaDlnx,vpure,SAFTSAC
    # from numba import config
    # config.DISABLE_JIT = True
    # T=298.15
    # p=1E5
    # npoint=12000
    # mi=np.asarray([1.20469,1045.6])
    # sigi=np.asarray([2.797059952,2.71])
    # ui=np.asarray([353.95,205.599])
    # epsAiBi=np.asarray([2425.67,0.])
    # kapi=np.asarray([0.04509,0.02])
    # N=np.asarray([1.,231.])
    # vpures=vpure(p,T,mi,sigi,ui,epsAiBi,kapi,N)
    # kij=D_Matrix(np.asarray([-0.148]),nc)
    # for i in range(1):
    #     #plt.plot(wt[:,0])
    #     Gammai=np.asarray([DlnaDlnx(T,vpures,np.ascontiguousarray(wt[i,:]),mi,sigi,ui,epsAiBi,kapi,N,Mi,kij) for i in range(nt)]).T
    #     wt=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,volatile,Gammai=Gammai)
    # ww=np.linspace(0,1,nt)
    # wt=np.stack((ww,1-ww)).T


    # Gammai2=np.asarray([DlnaDlnx(T,vpures,np.ascontiguousarray(wt[i,:]),mi,sigi,ui,epsAiBi,kapi,N,Mi,kij,idx=-1) for i in range(nt)]).T
    # lngammai=np.asarray([SAFTSAC(T,vpures,np.ascontiguousarray(wt[i,:]),mi,sigi,ui,epsAiBi,kapi,N,Mi,kij)[0] for i in range(nt)])
    # plt.plot(ww,Gammai2[0,0,:])

    # THFaktor1ave=np.average(Gammai2[0,0,:])
    # wtid=Diffusion_MS(t,L,Dvec*THFaktor1ave,wi0,wi8,Mi,volatile)
    # #plt.plot(wtid[:,0],'kx')

    # plt.show()


    # nc=3
    # nd=(nc-1)*nc//2

    # Dvec=np.asarray([1E-13,1E-13,1E-16])
    # #Dvec=np.asarray([1E-10,2.3E-10,3E-14])
    # #np.fill_diagonal(D,np.ones(nc)*1E-30)
    # L=0.0002
    # wi0=np.asarray([0.01,0.485,0.485])
    # wi8=np.asarray([0.1,0.45,0.45])
    # Mi=np.asarray([18.015,25700,357.79])
    # volatile=np.asarray([True,False,True])
    # wt=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,volatile)
    # from .PyCSAFT_nue import DlnaDlnx,vpure,SAFTSAC
    # T=298.15
    # p=1E5
    # npoint=12000
    # mi=np.asarray([1.20469,1045.6,14.283])
    # sigi=np.asarray([2.797059952,2.71,3.535])
    # ui=np.asarray([353.95,205.599,262.791])
    # epsAiBi=np.asarray([2425.67,0.,886.4])
    # kapi=np.asarray([0.04509,0.02,0.02])
    # N=np.asarray([1.,231.,3.])
    # vpures=vpure(p,T,mi,sigi,ui,epsAiBi,kapi,N)
    # kij=D_Matrix(np.asarray([-0.148,0.,0.]),nc)
    # for i in range(10):
    #     plt.plot(wt[:,0])
    #     Gammai=np.asarray([DlnaDlnx(T,vpures,np.ascontiguousarray(wt[i,:]),mi,sigi,ui,epsAiBi,kapi,N,Mi,kij) for i in range(nt)]).T
    #     #Gammai=np.stack([np.eye(nc)]*nt).T#*0.0001+Gammai*0.99999
    #     wt=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,volatile,Gammai=Gammai)
    # Gammai2=np.asarray([DlnaDlnx(T,vpures,np.ascontiguousarray(wt[i,:]),mi,sigi,ui,epsAiBi,kapi,N,Mi,kij,idx=-1) for i in range(nt)]).T
    # THFaktor1ave=np.average(Gammai2[0,0,:])#np.asarray([np.average(Gammai2[0,0,:]),np.average(Gammai2[0,1,:]),np.average(Gammai2[1,1,:])])
    # print(THFaktor1ave)
    # wtid=Diffusion_MS(t,L,Dvec*THFaktor1ave,wi0,wi8,Mi,volatile)
    # plt.plot(wtid[:,0],'kx')

    # plt.show()


    from .PyCSAFT_nue import DlnaDlnx,vpure,SAFTSAC
    T=298.15
    p=1E5
    texp=np.asarray([0,4.427792916,14.50035208,23.87257753,33.76909653,45.58674953,58.69408811,71.80142669,91.44521324,120.9089796,147.0930411,177.8656278,212.5493678,264.9136638,429.8518201])

    wpvac=np.asarray([0.333333333,0.397028757,0.559683846,0.704648614,0.763915741,0.791401827,0.801167219,0.822653035,0.830878575,0.841458325,0.846286716,0.857199094,0.851138308,0.858388537,0.8723936549])

    wtol=np.asarray([0.333333333,0.361083897,0.34858624,0.273123698,0.22781745,0.202716913,0.188759833,0.17692216,0.169121425,0.158541675,0.153713284,0.142800906,0.143816828,0.137082518,0.123003508])

    wmet=np.asarray([0.333333333,0.241887346,0.091729914,0.022227687,0.008266808,0.00588126,0.010072947,0.000424805,0,0,0,0,0,0,0])

    mi=np.asarray([1.5255, 2.8149, 2889.9])
    sigi=np.asarray([3.2300, 3.7169, 3.3972])
    ui=np.asarray([188.9, 285.69, 204.65])
    kapi=np.asarray([0.035176, 0., 0.])
    epsAiBi=np.asarray([2899.5, 0., 0.])
    N=np.asarray([1., 0., 1047.])
    Mw=np.asarray([32.042,  92.142, 90000.])
    vpures=vpure(p,T,mi,sigi,ui,epsAiBi,kapi,N)
    Dvec=np.asarray([1E-3,2.3E-11,1.7E-11])
    nt=300
    nc=3
    L=2E-5
    wi0=np.asarray([0.333333333,0.333333333,0.333333333])
    wi8=np.asarray([0.00001,0.127606346,0.872393654])
    Mi=np.asarray([32.042,92.142,90000.])
    t=np.linspace(0,texp[-1],nt)
    volatile=np.asarray([True,True,False])
    wt=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,volatile)
    kij=D_Matrix(np.asarray([0.029,-0.05855362,0.027776682]),nc)

    for i in range(10):
        plt.plot(wt[:,0])
        Gammai=np.asarray([DlnaDlnx(T,vpures,np.ascontiguousarray(wt[i,:]),mi,sigi,ui,epsAiBi,kapi,N,Mw,kij) for i in range(nt)]).T
        Gammai2=np.asarray([DlnaDlnx(T,vpures,np.ascontiguousarray(wt[i,:]),mi,sigi,ui,epsAiBi,kapi,N,Mw,kij,idx=-1) for i in range(nt)]).T #leads to the same
        #Gammai=(Gammai*i/5+np.stack([np.eye(nc)]*nt).T*(5-i)/5)
        wt=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,volatile,Gammai=Gammai)

    plt.show()