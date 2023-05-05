import numpy as np
from scipy.integrate import ode,odeint,solve_ivp
import time
import copy
import matplotlib.pyplot as plt

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
    nf=int(np.sum(volatile))
    nt=len(t)
    rho=1200
    refsegment=np.argmin(Mi)
    allflux=nc==nf
    isvolatile=np.where(volatile)[0]
    notvolatile=np.where(~volatile)[0]
    #initial conditions
    rhoiinit=np.zeros((nc,nz+1))
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
    xinit=np.reshape(rhovinit,(np.multiply(*rhovinit.shape),1))

    #decision variable vector
    #rhov=cs.MX.sym("rhov",(nf,nz+1))
    #Time=cs.MX.sym("Time")
    lambi=1 if taui is None else lambda Time : 1-np.exp(-Time/taui)
    dlambidt=0 if taui is None else lambda Time: 1/taui*np.exp(-Time/taui)

    #x=np.reshape(rhov,(np.multiply(*rhov.shape),1))
    
    
    GammaiT=np.eye(nc, dtype=object)
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
                    GammaiT[volatiles[i],volatiles[j]]= lambda Time,c=i,d=j :np.interp(Time,t,THij[c,d,:]-massbalancecorrection[d,c,:])    
    # if Gammai is not None:
    #     Gammai_2=copy.deepcopy(Gammai.T)
    #     Gammai_2[Gammai_2<=0]=1.
    #     if np.any(Gammai_2!=Gammai.T):
    #        print("This is a warning")
    def drhodt(t,rhov):
        rhoi=np.ones((nc,nz+1))
        rhov=np.reshape(rhov,(nf,nz+1))
        rhoi[isvolatile,:]=rhov
        rhoi[notvolatile,:]=rhoiinit[notvolatile,:]


        ji=np.zeros((nc,nz))
        dlnwi=np.zeros((nc,nz))
        wibar=np.ones((nc,nz))
        rhoibar=np.ones((nc,nz))
        drhoidt=np.ones((nc,nz+1))
        wi=np.zeros((nc,nz+1))

        ri= Mi/Mi[refsegment]


        for i in range(nc):
            wi[i,:]=rhoi[i,:]/np.sum(rhoi,axis=0)
            
    
        for i in range(nc):
            #dlnwi[i,:]=np.diff(np.log(wi[i,:]))*GammaiT(t)
            Gammacorr=GammaiT[i,i](t) if callable(GammaiT[i,i]) else GammaiT[i,i]
            #print(Gammacorr)
            dlnwi[i,:]=np.diff(np.log(wi[i,:]))*Gammacorr # here it would be for every component
            #dlnwi[i,:]=np.sum(np.diff(np.log(wi),1,1)*GammaiT[i,:].T)# here it would be for every component
            wibar[i,:]= averaging(wi[i,:])
            rhoibar[i,:]= averaging(rhoi[i,:])
        for z in range(nz):
            B=BIJ_Matrix(D,wibar[:,z],volatile) 
            Binv=np.linalg.inv(B) if len(B.shape)>1 else B**-1
            dmui=(dlnwi[:,z])
            #di=rhoibar[:,z]*dmui/ri # in rhoibar steckt die transformierte riesige Dichte drin. Der Volumenexpansionsfaktor korriegiert diesen
            omega=rho/np.sum(rhoibar[:,z])
            di=rho*wibar[:,z]*dmui/ri*omega if swelling else rhoibar[:,z]*dmui/ri
            #di=rhoibar[:,z]*dmui/ri*omega**2
            
            if not allflux:
                ji[isvolatile,z]=di[isvolatile].T@Binv
            else:
                ji[:-1,z]=di[:-1].T@Binv

        for i in range(nc):
            dji=np.diff(np.hstack((0,ji[i,:])))
            djib=np.hstack((dji,0)) 
            drhoidt[i,:]=djib

        X0=w0[isvolatile]/(1-np.sum(w0[isvolatile])) if not allflux else w0
        X8=w8[isvolatile]/(1-np.sum(w8[isvolatile])) if not allflux else w8
        #drhoidt[isvolatile,-1]=dlambidt*(X8-X0)*rho
        drhoidt[isvolatile,-1]=dlambidt*(w8[isvolatile]-w0[isvolatile])*rho
        # need for boundary when swelling
        # drhoidtmean=cs.sum2(drhoidt[:,:-1])/nz
        # drhoidt[:,-1]=cs.sum1(drhoidtmean)*w8
        #is this it

        drhovdt=drhoidt[isvolatile,:]
        ode=np.reshape(drhovdt,np.multiply(*drhovdt.shape))
        return ode

    xinit=np.reshape(rhovinit,np.multiply(*rhovinit.shape))
    drhodt(t[0],xinit)
    
    opts = {"grid":t,"regularity_check":True,"output_t0":True,"verbose":True}
    
    r=ode(drhodt).set_integrator('vode', method='adams',order=12)
    
    r.set_initial_value(xinit,0.)
    start=time.time_ns()
    #x_sol=odeint(drhodt,xinit,t,tfirst=True).T
    #x_sol=np.asarray([xinit]+[r.integrate(val) for val in t[1:]]).T
    x_sol=solve_ivp(drhodt,(t[0],t[-1]),xinit,method="Radau",t_eval=t)["y"]
    end=time.time_ns()
    print("------------- Diffusion modeling took "+str((end-start)/1E9)+" seconds ----------------")
    nt=len(t)
    wt=np.zeros((nc,nt))
    wik=np.zeros((nt,nc,nz+1))
    rhoik=np.zeros((nt,nc,nz+1))
    rhok=np.zeros(nt)
    for k in range(nt):
        rhovk=np.reshape(x_sol[:,k],(nf,nz+1))
        rhoik[k,:,:]=rhoiinit
        rhoik[k,isvolatile,:]=rhovk
        for i in range(nc):
            wik[k,i,:]=rhoik[k,i,:]/np.sum(rhoik[k,:,:],axis=0)
        rhok[k]=np.sum(np.sum(rhoik[k,:,:-1]/nz,axis=1),axis=0)
        wt[:,k]=np.sum(wik[k,:,:-1]/nz,axis=1)
    Lt=rhok/rho*L
    #[plt.plot(zvec,rhoik[k,1,:]) for k,val in enumerate(rhoik[:,1,0])]
    #plt.show()
    return wt.T if not full_output else (wt.T,wik,zvec,Lt)


def BIJ_Matrix(D,wi,volatile):
    nc=wi.shape[0]
    nf=np.sum(volatile)
    allflux=nc==nf
    B=np.ones((nc,nc))
    comp=np.arange(0,nc)
    isvolatile=np.where(volatile)[0]
    notvolatile=np.where(~volatile)[0]
    for i in range(nc):
        Din=wi[-1]/D[i,-1] if (i+1)!=nc else 0
        Dij=-wi[i]/D[i,:]
        BIJ= Dij if not allflux else Dij+Din
        Dii=np.sum((wi/D[i,:])[comp[comp!=i]])
        BII=Dii if not allflux else Dii+Din
        B[i,:]=BIJ
        B[i,i]=BII
    #alpha=1/np.min(D) #implement an augmented transport matrix
    #B+alpha*wi@wi.T

    return B[isvolatile,:][:,isvolatile] if not allflux else B[:-1,:-1]  #8E-13

def D_Matrix(Dvec,nc):
    nd=(nc-1)*nc//2 
    if len(Dvec)!=nd: 
        raise Exception("Wrong number of diffusion coefficients. Provide array with "+str(nd)+" entries")
    else:
        D=np.zeros((nc,nc))
        D[np.triu_indices_from(D,k=1)]=Dvec
        D[np.tril_indices_from(D,k=-1)]=Dvec
    return D

if __name__=="__main__":
    nt=100
    t=np.linspace(0,429.8518201,nt)
    Dvec=np.asarray([1E-3,2.3E-11,1.7E-11])
    nc=3
    L=2E-5
    wi0=np.asarray([0.333333333,0.333333333,0.333333333])
    wi8=np.asarray([0.00001,0.127606346,0.872393654])
    Mi=np.asarray([32.042,92.142,90000.])
    volatile=np.asarray([True,True,False])
    wt=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,volatile)
    plt.plot(t,wt[:,0])
    plt.plot(t,wt[:,1])
    plt.plot(t,wt[:,2])


    from .PyCSAFT_nue import DlnaDlnx,vpure,SAFTSAC
    T=298.15
    p=1E5


    mi=np.asarray([1.5255, 2.8149, 2889.9])
    sigi=np.asarray([3.2300, 3.7169, 3.3972])
    ui=np.asarray([188.9, 285.69, 204.65])
    kapi=np.asarray([0.035176, 0., 0.])
    epsAiBi=np.asarray([2899.5, 0., 0.])
    N=np.asarray([1., 0., 1047.])
    Mw=np.asarray([32.042,  92.142, 90000.])
    vpures=vpure(p,T,mi,sigi,ui,epsAiBi,kapi,N)

    kij=D_Matrix(np.asarray([0.029,-0.05855362,0.027776682]),nc)

    for i in range(5):
        plt.plot(t,wt[:,0])
        plt.plot(t,wt[:,1])
        plt.plot(t,wt[:,2])
        Gammai=np.asarray([DlnaDlnx(T,vpures,np.ascontiguousarray(wt[i,:]),mi,sigi,ui,epsAiBi,kapi,N,Mw,kij) for i in range(nt)]).T
        Gammai2=np.asarray([DlnaDlnx(T,vpures,np.ascontiguousarray(wt[i,:]),mi,sigi,ui,epsAiBi,kapi,N,Mw,kij,idx=-1) for i in range(nt)]).T #leads to the same
        #Gammai=(Gammai*i/5+np.stack([np.eye(nc)]*nt).T*(5-i)/5)
        wt=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,volatile,Gammai=Gammai)
    # plt.plot(t,wt[:,0])
    # plt.plot(t,wt[:,1])
    # plt.plot(t,wt[:,2])
    plt.show()