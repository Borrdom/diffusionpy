import numpy as np
from scipy.integrate import solve_ivp

def time_dep_surface_cryst(t,mobiles,immobiles,crystallize,wi0,wi8,rho0i,Mi,DAPI,sigma,kt,g,deltaHSL,TSL,cpSL,lngi_fun=None,wv_fun=None):
    """calculate the time dependent surface concentration during crystallization

    Args:
        t (array_like): vector of time
        mobiles (array_like): index array indicating the mobile components
        immobiles (array_like): index array indicating the immobile component
        crystallize (array_like): index array indicating the crystallizing components
        wi0 (array_like): initial mass fractions
        wi8 (array_like): mass fractions at time equals infinity
        rho0i (array_like): pure component densities
        Mi (array_like): molar mass of components
        DAPI (array_like): crystallizing components diffusion coefficient in the vector
        sigma (array_like): interfacial tension of crystal component and mixture
        kt (array_like): crystal growth rate constant
        g (array_like): crsystal growth exponent
        deltaHSL (array_like): melting enthalpy of crystallizing components
        TSL (array_like): melting temperature of crystallizing components
        cpSL (array_like): differfence in liquid/solid heat capacity of crystallizing components
        lngi_fun (array_like, optional): function of logarithmic activity coefficients
        wv_fun (array_like, optional): function or vector how the concentration of the volatile components changes with the 
        concentration of the cryystallizing components

    Returns:
        array_like: vector of mass fractions at the surface as a function of time
    """
    crystallizes=np.where(crystallize)[0]
    M=Mi[crystallizes]/1000.
    rho=rho0i[crystallizes]
    AR=100
    pre=rho*np.pi/(4*AR**2)
    R=8.31445
    NA=6.023E23
    kB=R/NA
    C0=rho/M*NA 
    temp=298.15
    lnai=lngi_fun(wi8)+np.log(wi8)
    lnaiSLE=-deltaHSL/(R*temp)*(1-temp/TSL)+cpSL/R*(TSL/temp-1-np.log(TSL/temp))
    dmu_sla0=lnai[crystallizes]-lnaiSLE
    r0=2*sigma/(C0*dmu_sla0*kB*temp)
    alpha0=pre*(r0)**3
    def dxdt(t,x):
        alpha,r=x[0],x[1]
        dalphadt,dalphadr=CNT(t,alpha,r,mobiles,immobiles,crystallize,wi0,wi8,rho0i,Mi,DAPI,sigma,kt,g,deltaHSL,TSL,cpSL,lngi_fun,wv_fun)
        return np.hstack((dalphadt,dalphadr))
    x0=np.hstack((alpha0,r0))
    xsol=solve_ivp(dxdt,(t[0],t[-1]),x0,method="Radau",t_eval=t)["y"]
    alpha,r=xsol[0],xsol[1]
    dl0=wi0[crystallizes]/np.sum(wi0[immobiles])
    dl_la = (dl0-alpha)/(1-alpha)
    wvB=wv_fun(dl_la)
    wiB=np.zeros((wi0.shape[0],t.shape[0]))
    wiB[mobiles,:]=wvB
    wiB[crystallizes,:]=(1-wvB)*dl_la
    wiB[immobiles[crystallizes!=immobiles],:]=(1-wvB)*(1-dl_la)
    return wiB,alpha,r

def CNT(t,alpha,r,mobiles,immobiles,crystallize,wi0,wi8,rho0i,Mi,DAPI,sigma,kt,g,deltaHSL,TSL,cpSL,lngi_fun=None,wv_fun=None):
    """Calculates the crystallization kinetics based on the classical nucleation theory for nucleation and a simple crystal growth model"""
    crystallizes=np.where(crystallize)[0]
    R=8.31445
    NA=6.023E23
    kB=R/NA
    temp=298.15
    AR=100
    M=Mi[crystallizes]/1000.
    rho=rho0i[crystallizes]
    pre=rho*np.pi/(4*AR**2)
    C0=rho/M*NA
    dl0=wi0[crystallizes]/np.sum(wi0[immobiles])
    X_la=dl0-alpha
    Xn_la=X_la/M
    wi=np.zeros_like(wi0)
    dl_la = (dl0-alpha)/(1-alpha)

    wv=wv_fun(dl_la) if callable(wv_fun) else wv_fun
    wi[mobiles]=wv
    wi[crystallizes]=(1-wv)*dl_la
    wi[immobiles[crystallizes!=immobiles]]=(1-wv)*(1-dl_la)
    lnaiSLE=-deltaHSL/(R*temp)*(1-temp/TSL)+cpSL/R*(TSL/temp-1-np.log(TSL/temp))
    lnai=lngi_fun(wi)+np.log(wi)
    dmu_sla=lnai[crystallizes]-lnaiSLE
    rstar = ((2*sigma)/(C0*kB*temp*dmu_sla))
    deltaG=sigma**3*(16*np.pi)/(3*(C0*kB*temp*dmu_sla)**2)
    wv=wi[mobiles]
    wv0=wi0[mobiles]
    wv8=wi8[mobiles]
    beta=np.fmin((wv-wv0)/(wv8-wv0),1)
    ze=(kB*temp/sigma)**(1.5)*C0/(8*np.pi)*dmu_sla**2
    f=4*np.pi*rstar*DAPI*Xn_la*NA*np.fmax(beta,1E-4)
    dNdt = ze*f*C0*np.exp(-deltaG/(kB*temp))#*cs.exp(-NA)/NA**0.5
    drdt = np.fmax(kt*np.fmax(beta,1E-4)*(dmu_sla)**g,0)
    dalphadr=3*alpha/r
    dalphadN=pre*(r)**3
    dalphadt=np.fmax((drdt*dalphadr+dNdt*dalphadN),0)
    return dalphadt,drdt

