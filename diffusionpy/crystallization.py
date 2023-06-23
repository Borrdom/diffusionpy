import numpy as np
from numba import njit
from scipy.interpolate import interp1d
from.PyCSAFT_nue import lngi,vpure,dlnai_dlnxi
from.Stefan_Maxwell_segmental import D_Matrix,Diffusion_MS_iter,Diffusion_MS
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def Bound(t,mobiles,immobiles,crystallize,wi0,wi8,rho0i,Mi,DAPI,sigma,kt,g,deltaHSL,TSL,cpSL,lngi_fun=None,wv_fun=None):
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
        dalphadt,dalphadr=Kris(t,alpha,r,mobiles,immobiles,crystallize,wi0,wi8,rho0i,Mi,DAPI,sigma,kt,g,deltaHSL,TSL,cpSL,lngi_fun,wv_fun)
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

def Kris(t,alpha,r,mobiles,immobiles,crystallize,wi0,wi8,rho0i,Mi,DAPI,sigma,kt,g,deltaHSL,TSL,cpSL,lngi_fun=None,wv_fun=None):
    """function calculating the crystallization"""
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






crystallize=np.asarray([False,False,True])
mobile=np.asarray([True,False,False])
mobiles=np.where(mobile)[0]
immobiles=np.where(~mobile)[0]
deltaHSL=np.asarray([31500.])
TSL=np.asarray([429.47])
cpSL=np.asarray([87.44])
ww=np.asarray([0.27087,0.22302, 0.13792, 0.09208, 0.06118])
dl=np.asarray([0,0.1 , 0.3 , 0.5 , 0.68])
DAPI=np.asarray([6.6E-17])
sigma=np.asarray([1.98E-02])
kt=np.asarray([5.1E-12])
g=np.asarray([3.2])
Mi=np.asarray([18.015,65000.,230.26])
rho0i=np.asarray([997.,1180.,1320.])

wv_fun=interp1d(dl,ww,bounds_error=False,fill_value=(0.27087,0.06118))

nc=3
wv0=0.0001
dl0=0.68
wi0=np.asarray([wv0,(1-wv0)*(1-dl0),(1-wv0)*dl0])
wv8=wv_fun(dl0)
wi8=np.asarray([wv8,(1-wv8)*(1-dl0),(1-wv8)*dl0])
T=298.15
p=1E5

kij=D_Matrix(np.asarray([-0.128,0.00648,-0.0574]),nc)
par={"mi":np.asarray([1.20469,2420.99, 8.105152]),
"si": np.asarray([2.797059952,2.947, 2.939]),
"ui" :np.asarray([353.95,205.27, 229.45]),
"eAi" :np.asarray([2425.67,0., 934.2]),
"kAi":np.asarray([0.04509,0.02, 0.02]),
"NAi":np.asarray([1.,653., 2.]),
"Mi": Mi,
"kij":kij}

vpures=vpure(p,T,**par)

par["vpure"]=vpures
lngi_fun=lambda wi: lngi(T,wi,**par)


nt=300
t=np.linspace(0,1000,nt)*60
witB,alpha,r=Bound(t,mobiles,immobiles,crystallize,wi0,wi8,rho0i,Mi,DAPI,sigma,kt,g,deltaHSL,TSL,cpSL,lngi_fun,wv_fun)
tmin=t/60
# fig1,ax1=plt.subplots()
# fig2,ax2=plt.subplots()
# ax1.plot(tmin,alpha)
# rmu=r/1E-6
# ax2.plot(tmin,rmu)
# import pandas as pd
# pd.DataFrame((tmin,alpha,rmu)).to_clipboard()
# plt.plot(t/60,witB[0,:])
# plt.plot(t/60,witB[1,:])
# plt.plot(t/60,witB[2,:])
Dvec=np.asarray([2E-13,2E-13,2E-13])
L=2.5E-5
dlnai_dlnwi_fun=lambda wi: dlnai_dlnxi(T,wi,**par)
wt=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,swelling=True,witB=witB.T)
# wt=Diffusion_MS_iter(t,L,Dvec,wi0,wi8,Mi,mobile,swelling=True,witB=witB.T,dlnai_dlnwi_fun=dlnai_dlnwi_fun)
XwL=wt[:,0]/(1-wt[:,0])
Xw=XwL*(1-alpha)
# plt.plot(t/60,witB.T[:,0],"bo")
plt.plot(t/60,Xw,'b',label="Non-ideal and wiB(t)")
plt.legend("")
plt.show()
import pandas as pd
pd.DataFrame((tmin,Xw)).T.to_clipboard()