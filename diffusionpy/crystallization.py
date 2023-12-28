import numpy as np
# from numba import njit
from scipy.integrate import solve_ivp
from numba import njit,config,prange,jit
import time


# config.DISABLE_JIT = True
# @njit(['Tuple((f8[:,::1], f8[:,::1]))(f8, f8[::1], f8[::1],  i8[::1], i8[::1],i8[::1], f8[::1], f8[::1],f8[::1],f8[::1], f8[::1], f8[::1], f8[::1], f8[::1],f8[::1],f8[::1],f8[::1],f8[:,::1],f8[:,::1])'],cache=True)
def CNT(t,alpha,r,mobiles,immobiles,crystallizes,wi0,wi8,rho0i,Mi,DAPI,sigma,kt,g,deltaHSL,TSL,cpSL,tnuc,temp,lngi,wv):
    """Calculates the crystallization kinetics based on the classical nucleation theory for nucleation and a simple crystal growth model"""
    
    R=8.31445
    NA=6.023E23
    kB=R/NA
    AR=100
    M=Mi[crystallizes]/1000.
    rho=rho0i[crystallizes]
    pre=rho*np.pi/(4*AR**2)
    C0=rho/M*NA
    X_la=1-alpha
    Xn_la=X_la/M
    nz_1=alpha.shape[0] if alpha.shape is not () else 1
    nTH=np.sum(mobiles)
    nc=len(wi0)
    wi=np.zeros((nc,nz_1))
    dl_la = (1-alpha)/(1/wi0[crystallizes]-alpha)
    # wv=wvfun(dl_la) if np.shape(wvfun)[0]==0 else wvfun
    wi[mobiles,...]=wv
    if immobiles.shape[0]>0.:
        wi[crystallizes,...]=(1-np.sum(wv,axis=0))*dl_la
        wi[immobiles[crystallizes!=immobiles],...]=(1-np.sum(wv,axis=0))*(1-dl_la)
        beta=np.fmin(np.fmax((wv[mobiles[0],:]-wi0[mobiles[0]])/(wi8[mobiles[0]]-wi0[mobiles[0]]),0),1)
    else:
        wi[crystallizes,...]=dl_la
        wi[mobiles,...]=wi[mobiles,...]/np.sum(wi[mobiles,...],axis=0)
        beta=np.ones_like(alpha)
    lnaiSLE=-deltaHSL/(R*temp)*(1-temp/TSL)+cpSL/R*(TSL/temp-1-np.log(TSL/temp))
    lnai=lngi.T+np.log(wi)
    dmu_sla=np.fmax(lnai[crystallizes,...]-lnaiSLE,1E-20)
    rstar = ((2*sigma)/(C0*kB*temp*dmu_sla))
    deltaG=sigma**3*(16*np.pi)/(3*(C0*kB*temp*dmu_sla)**2)
    ze=(kB*temp/sigma)**(1.5)*C0/(8*np.pi)*dmu_sla**2
    f=4*np.pi*rstar*DAPI*Xn_la*NA
    A=((2*sigma)/(C0*kB*temp))*4*np.pi*((2*sigma)/(C0*kB*temp))*DAPI*NA*(kB*temp/sigma)**(1.5)*C0/(8*np.pi)*C0
    B=(16*np.pi)/(3*(C0)**2)*(sigma/(kB))**3
    
    
    dNdt = np.fmax(beta*ze*f*C0*np.exp(-deltaG/(kB*temp)),0) #*cs.exp(-NA)/NA**0.5
    r=np.fmax(r,rstar)
    drdt = np.fmax(beta*kt*(dmu_sla)**g,0) 
    drdt=drdt if (t/60)>tnuc else np.zeros_like(drdt)
    # drdt = np.fmax(kt*temp*(1-1/np.exp(dmu_sla)),0)
    dalphadr=3*alpha/r
    dalphadN=pre*(r)**3
    dalphadt=drdt*dalphadr+dNdt*dalphadN 
    return dalphadt,drdt

# config.DISABLE_JIT = True
def crystallization_mode(wvinit,ode,mobiles,immobiles,crystallize,wi0,wi8,rho0i,Mi,deltaHSL,TSL,cpSL,tnuc,temp,DAPI,sigma,kt,g,lngi_tz):
    """alter the ode function in diffusionpy.Diffusion_MS, to also solve the crystallization

    Args:
        wvinit (array_like): vector of the mass fractions of the mobile components
        ode (array_like): ode fuinction which is modified by the function
        deltaHSL (array_like): Melting Enthalpy
        TSL (array_like): Melting temperature
        DAPI (array_like) Diffusioncoefficient of the crystallizing substance in its surroundings
        sigma (array_like): interfacial tension of crystal to its surrounding
        kt (array_like): rate of crystal growth kinetics
        g (array_like): order of crystal growth kinetics

    Returns:
        array_like: new modified ode function with the same format as the input ode function
    """

    _,nz_1=wvinit.shape
    crystallizes=np.where(crystallize)[0]
    M=Mi[crystallizes]/1000.
    rho=rho0i[crystallizes]
    dl0=wi0[crystallizes]/np.sum(wi0[immobiles])
    CNT1=njit(['Tuple((f8[:,::1], f8[:,::1]))(f8, f8[::1], f8[::1],  i8[::1], i8[::1],i8[::1], f8[::1], f8[::1],f8[::1],f8[::1], f8[::1], f8[::1], f8[::1], f8[::1],f8[::1],f8[::1],f8[::1],f8,f8,f8[:,::1],f8[:,::1])'],cache=True)(CNT)
    def crystallization_ode(t,x,tint,THFaktor,mobiles,immobiles,Mi,D,allflux,wi0,dmuext,wiB):
        """solves the genralized Maxwell model for relaxation"""
        nTH,nz_1=dmuext.shape
        wv=np.zeros((nTH,nz_1))
        for i in range(nTH):
            wvtemp=x[(nz_1)*(i):(nz_1)*(1+i)]
            wv[i,:]=wvtemp
        # wv=np.fmin(np.fmax(wv,1E-300),1.)
        alpha=np.fmax(x[(nz_1)*(nTH):(nz_1)*(nTH+1)],1E-30)
        r=x[(nz_1)*(nTH+1):(nz_1)*(nTH+2)]
        # if immobiles.shape[0]==0.:
        #     dl_la = (1-alpha)/(1/wi0[crystallizes]-alpha) 
        #     wv[crystallizes,...]=(1-np.sum(wv,axis=0))*dl_la  if immobiles.shape[0]>0. else dl_la
        #     wv[mobiles,...]=wv[mobiles,...]/np.sum(wv[mobiles,...],axis=0)
        wv=np.ascontiguousarray(wv)
        for i in range(nTH):
            wv[i,-1]=np.interp(t,tint,wiB[:,mobiles[i]])
            # wv[i,:]=np.interp(t,tint,wiB[:,mobiles[i]])*np.ones_like(wv[i,:])
        # omega=(1-np.sum(wv,axis=0))/np.sum(wi0[immobiles])
        # porosity=(1-alpha*omega)[None,:,None,None]
        # eta=1.5
        # dwvdt=ode(t,np.ascontiguousarray(wv.flatten()),tint,THFaktor*porosity**eta,mobiles,immobiles,Mi,D,allflux,wi0,dmuext,wiB)
        dwvdt=ode(t,np.ascontiguousarray(wv.flatten()),tint,THFaktor,mobiles,immobiles,Mi,D,allflux,wi0,dmuext,wiB)
        
        dalphadt,drdt=CNT1(t,np.ascontiguousarray(alpha),np.ascontiguousarray(r),mobiles,immobiles,crystallizes,wi0,wi8,rho0i,Mi,DAPI,sigma,kt,g,deltaHSL,TSL,cpSL,tnuc,temp,lngi_tz(t),wv)
        fvec=np.hstack((dwvdt.flatten(),dalphadt.flatten(),drdt.flatten()))
        return fvec

    AR=100
    pre=rho*np.pi/(4*AR**2)
    R=8.31445
    NA=6.023E23
    kB=R/NA
    C0=rho/M*NA 
    logwi0tz=np.ones(nz_1)*np.log(wi0[crystallizes])
    logwi0tz[-1]=np.log(wi8[crystallizes])
    lnai=lngi_tz(0)[:,crystallizes].flatten()+logwi0tz
    lnaiSLE=-deltaHSL/(R*temp)*(1-temp/TSL)+cpSL/R*(TSL/temp-1-np.log(TSL/temp))
    dmu_sla0=lnai-lnaiSLE 
    r0=2*sigma/(C0*dmu_sla0*kB*temp)*np.ones(nz_1)
    alpha0=pre*(r0)**3
    xinit=np.hstack((wvinit.flatten(),alpha0.flatten(),r0.flatten()))
    return xinit,crystallization_ode

def time_dep_surface_cryst(t,mobiles,immobiles,crystallize,wi0,wi8,rho0i,Mi,DAPI,sigma,kt,g,deltaHSL,TSL,cpSL,tnuc=0.,temp=298.15,lngi=None,wv_fun=None):
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
    logwi0tz=np.log(wi0[crystallizes])
    lnai=lngi(0)[crystallizes].flatten()+logwi0tz
    lnaiSLE=-deltaHSL/(R*temp)*(1-temp/TSL)+cpSL/R*(TSL/temp-1-np.log(TSL/temp))
    dmu_sla0=lnai-lnaiSLE
    r0=2*sigma/(C0*dmu_sla0*kB*temp)
    alpha0=pre*(r0)**3
    def dxdt(t,x):
        alpha,r=x[0],x[1]
        dl_la = (1-alpha)/(1/wi0[crystallizes]-alpha)
        wv=wv_fun(dl_la)
        dalphadt,dalphadr=CNT(t,alpha,r,mobiles,immobiles,crystallizes,wi0,wi8,rho0i,Mi,DAPI,sigma,kt,g,deltaHSL,TSL,cpSL,tnuc,temp,lngi(t)[None,:],wv)
        return np.hstack((dalphadt[0],dalphadr[0]))
    x0=np.hstack((alpha0,r0))
    xsol=solve_ivp(dxdt,(t[0],t[-1]),x0,method="Radau",t_eval=t)["y"]
    alpha,r=xsol[0],xsol[1]
    dl0=wi0[crystallizes]/np.sum(wi0[immobiles])
    dl_la = (1-alpha)/(1/wi0[crystallizes]-alpha)
    wvB=wv_fun(dl_la)
    wiB=np.zeros((wi0.shape[0],t.shape[0]))
    wiB[mobiles,:]=wvB
    wiB[crystallizes,:]=(1-wvB)*dl_la
    wiB[immobiles[crystallizes!=immobiles],:]=(1-wvB)*(1-dl_la)
    return wiB,alpha,r

