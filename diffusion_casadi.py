import numpy as np
from scipy.integrate import solve_ivp,odeint
import matplotlib.pyplot as plt
import time
#from PyCSAFT_2 import mixture
from get_param import get_par
import os
from epcsaftpy import pcsaft,component
import casadi as cs

def dictkij2mat(dictkij,scomp):
    ncomp=len(scomp)
    kij=np.zeros((ncomp,ncomp))
    for i in range(ncomp):
        co1=scomp[i]
        for j in range(ncomp):
            if i!=j:
                try:
                    kij[i,j]=dictkij[co1+scomp[j]] if co1+scomp[j] in dictkij.keys() else dictkij[scomp[j]+co1]
                except:
                    kij[i,j]=0
    return kij
def get_mixture(components,T=298.15,path=os.path.join(os.getcwd(),"Nil","parameter.accdb")):
    df_par,dikij,_=get_par(components,T=T,path=path)
    a=[]
    for i,val in enumerate(df_par):
        a.append(component(name=val["name"],
        ms=val["mi"],
        Mw=val["Mi"]*1000,
        sigma=val["sigi"],
        eps=val["ui"],
        eAB=val["epsAiBi"],
        kappaAB=val["kapi"],
        sites=[0,val["N"],val["N"]]))
        if i>0:
            pars+=a[i] 
        else:
            pars=a[i]
    pars.KIJ0saft=dictkij2mat(dikij,components)
    eos1 = pcsaft(pars)
    return eos1
def get_lnphi(mix1,T,p,xi=None,wi=None,index=1):
    Mi=np.atleast_2d(mix1.Mw).T
    xorw=xi is None
    xi=wi/Mi/np.sum(wi/Mi,axis=0) if xorw else xi
    xtow=np.log(xi/wi) if xorw else 0
    lnphi=np.asarray([mix1.logfugef(xi[:,i],T,p,"L")[0][index] for i in range(len(wi[0,:]))])
    return wi[index,:],lnphi+xtow[index,:]


def SolveODEs(x,f,x0,opts,Time):
    dae={"x":x,"t":Time,"ode":f}#,"print_stats":True,"abstol":1E-12}#,max_num_steps":1000"collocation_scheme":"legendre","interpolation_order":2,"number_of_finite_elements":100,"rootfinder_options":{"error_on_fail": False}
    fint=cs.integrator("integrator","cvodes",dae,opts)
    F=fint(x0=x0)
    return F["xf"]  


def averaging(a):
    return (a[:-1]+a[1:])/2

def CDF(a):
    return (a[:-1]-a[1:])

def crank(t,d,l):
    pi = np.pi
    ns=30
    mt_minf = 1-sum([(8/pi**2)*(1/(2*n + 1)**2)*np.exp((-d*((n + 1/2)**2)*(pi**2)*t)/l**2) for n in range(ns)])
    return mt_minf

def Diffusion1D(t,L0,Ds,ws0,ws8,lnphs_fun=None):
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
    epsilon=Xsbar
    swelling=(1+epsilon)**-2
    Ds=Ds0*swelling/(1-wsbar) if ws8<0.95 else Ds0*swelling

    if lnphs_fun is None:
        lnphsint=np.asarray([0.,0.])
        wsint=np.asarray([0.,1.])
        lnphs_fun=cs.interpolant("lnphis_fun","linear",[wsint],lnphsint)
    
    lnfs= np.log(ws)+lnphs_fun(ws) 
    TDF= CDF(lnfs)
    j_bnd=np.asarray([0.])
    js= cs.vertcat(j_bnd,Xsbar*Ds*TDF)
    dXsdt_bnd=np.asarray([0.])
    dXsdt=cs.vertcat(CDF(js),dXsdt_bnd)
    opts = {"grid":t,"max_num_steps":10000,"regularity_check":True,"output_t0":True}
    start=time.time_ns()
    Xs_sol=SolveODEs(Xs,dXsdt,Xs_ini,opts,Time).full()
    end=time.time_ns()
    print((end-start)/1E9)
    Xst=np.sum(Xs_sol[:-1,:],axis=0)/(nz)
    wst=Xst/(1+Xst)
    return wst
def lnphi_int(index,DL,mix1):
    T=298.15
    p=1E5
    index=0
    ws=np.linspace(1E-26,1,1000)
    wp=(1-DL)*(1-ws)
    wa=(DL)*(1-ws)
    wi=np.stack((ws,wp,wa))
    wsint,lnphiint=get_lnphi(mix1,T,p,wi=wi,index=index)
    lnphi_fun=cs.interpolant("lnphis_fun","linear",[wsint],lnphiint)
    return lnphi_fun
def THFaktor(ws,lnphi_fun):
    x=cs.MX.sym("x")
    return cs.Function("TH",[x],[1+x*cs.gradient(lnphi_fun(x),x)])(ws)#ws*(lnphi_fun(ws+1E-8)-lnphi_fun(ws-1E-8))/2E-8+1

if __name__=="__main__":
    nt=50
    t=np.linspace(0,60,nt)*60
    L0=10E-6
    Ds0=1E-13
    ws0=0.1
    ws8=0.3
    components=["water","pvpva64","indomethacin"]
    DL=0.
    #mix1=get_mixture(components)
    mix1=get_mixture(components)
    lnphi_fun=lnphi_int(0,DL,mix1)

    wst1=Diffusion1D(t,L0,Ds0,ws0,ws8)

    wst1=Diffusion1D(t,L0,Ds0,ws0,ws8)

    
    wst2=Diffusion1D(t,L0,Ds0,ws0,ws8,lnphi_fun)

    plt.plot(t,wst1)
    plt.plot(t,wst2)
    plt.show()

