import numpy as np
from numba import njit
from scipy import interp1d
from.PyCSAFT_nue import lngi


crystallize=np.asarray([False,False,True])
deltaHSL=np.asarray([39300.])
TSL=np.asarray([433.25])
cpSL=np.asarray([116.95])
ww=np.asarray([0.22302, 0.13792, 0.09208, 0.06118])
dl=np.asarray([0.1 , 0.3 , 0.5 , 0.68])
DAPI=np.asarray([0.59010330e-19])
sigma=np.asarray([2.97730286e-02])
kt=np.asarray([8.0778700e-13])
g=np.asarray([3.92])
Mi=np.asarray([18.015,25700.,300.])
rho0i=np.asarray([997.,1180.,1320.])

ww_fun=np.interp1(dl,ww)

def Bound(wv,alpha,r,mobiles,immobiles,crystallize,wi0,wi8,rho0i,Mi,DAPI,sigma,kt,g,deltaHSL,TSL,cpSL,lngi_fun=None,wv_fun=None):
    crystallizes=np.where(crystallize)
    dalphadt,dalphadr=Kris(wv,alpha,r,mobiles,immobiles,crystallize,wi0,wi8,rho0i,Mi,DAPI,sigma,kt,g,deltaHSL,TSL,cpSL,lngi_fun)
    def solve(): 1,2
    alpha0=0
    r0=0
    alpha,r=solve(dalphadt,dalphadr,alpha0,r0)
    return wv

def Kris(alpha,r,mobiles,immobiles,crystallize,wi0,wi8,rho0i,Mi,DAPI,sigma,kt,g,deltaHSL,TSL,cpSL,lngi_fun=None,wv_fun=None):
    crystallizes=np.where(crystallize)
    R=8.31445
    NA=6.023E23
    kB=R/NA
    temp=298.15
    AR=100
    M=Mi[crystallizes]
    rho=rho0i[crystallizes]
    pre=rho*np.pi/(4*AR**2)
    scale_r=1E-6
    C0=rho/M*NA
    X_la=wi0[crystallizes]-alpha
    Xn_la=X_la/M
    wi=np.zeros_like(wi0)
    dl_la = (wi0[crystallizes]-alpha)/(1-alpha)
    wv=wv_fun(dl_la) if callable(wv_fun) else wv_fun
    wi[mobiles]=wv
    wi[immobiles]=(1-np.sum(wv))*wi0[immobiles]/np.sum(wi0[immobiles])
    lnaSLE=-deltaHSL/(R*temp)*(1-temp/TSL)+cpSL/R*(TSL/temp-1-np.log(TSL/temp))
    lnai=lngi_fun(wi)-np.log(wi)
    lnS_fun=lnai[crystallizes]-lnaSLE
    dmu_sla=lnS_fun(wi)
    rstar = ((2*sigma)/(C0*kB*temp*dmu_sla))
    deltaG=sigma**3*(16*np.pi)/(3*(C0*kB*temp*dmu_sla)**2)
    wv=wi[mobiles]
    wv0=wi0[mobiles]
    wv8=wi8[mobiles]
    beta=np.fmin((wv-wv0)/(wv8-wv0),1)
    ze=(kB*temp/sigma)**(1.5)*C0/(8*np.pi)*dmu_sla**2
    f=4*np.pi*rstar*DAPI*Xn_la*NA*np.fmax(beta,1E-4)
    dNdt = ze*f*C0*np.exp(-deltaG/(kB*temp))#*cs.exp(-NA)/NA**0.5
    drdt = np.fmax(1/scale_r*kt*np.fmax(beta,1E-4)*(dmu_sla)**g,0)
    dalphadr=3*alpha/r
    dalphadN=pre*(r*scale_r)**3
    dalphadt=np.fmax((drdt*dalphadr+dNdt*dalphadN),0)
    return dalphadt,dalphadr

def Krisnewerbut(wi,alpha,r,ww_fun,mobiles,crystallize,wi0,wi8,rho0i,Mi,DAPI,sigma,kt,g,deltaHSL,TSL,cpSL,lngi_fun=None):
    crystallizes=np.where(crystallize)
    R=8.31445
    NA=6.023E23
    kB=R/NA
    temp=298.15
    AR=100
    M=Mi[crystallizes]
    rho=rho0i[crystallizes]
    pre=rho*np.pi/(4*AR**2)
    scale_r=1E-6
    C0=rho/M*NA
    X_la=wi0[crystallizes]-alpha
    Xn_la=X_la/M
    dl_la = (wi0[crystallizes]-alpha)/(1-alpha)
    lnaSLE=-deltaHSL/(R*temp)*(1-temp/TSL)+cpSL/R*(TSL/temp-1-np.log(TSL/temp))
    lnai=lngi_fun(wi)-np.log(wi)
    lnS_fun=lnai[crystallizes]-lnaSLE
    dmu_sla=lnS_fun(wi)
    rstar = ((2*sigma)/(C0*kB*temp*dmu_sla))
    deltaG=sigma**3*(16*np.pi)/(3*(C0*kB*temp*dmu_sla)**2)
    wv=wi[mobiles]
    wv0=wi0[mobiles]
    wv8=wi8[mobiles]
    beta=np.fmin((wv-wv0)/(wv8-wv0),1)
    ze=(kB*temp/sigma)**(1.5)*C0/(8*np.pi)*dmu_sla**2
    f=4*np.pi*rstar*DAPI*Xn_la*NA*np.fmax(beta,1E-4)
    dNdt = ze*f*C0*np.exp(-deltaG/(kB*temp))#*cs.exp(-NA)/NA**0.5
    drdt = np.fmax(1/scale_r*kt*np.fmax(beta,1E-4)*(dmu_sla)**g,0)
    dalphadr=3*alpha/r
    dalphadN=pre*(r*scale_r)**3
    dalphadt=np.fmax((drdt*dalphadr+dNdt*dalphadN),0)
    rhovB=ww_fun(dl_la)
    return

def Krisold(rhov,alpha,r,ww_fun,mobiles,immobiles,crystallizes,wi0,rho,M,DAPI,sigma,kt,g,deltaHSL,TSL,cpSL,mnull,lngi_fun=None):
    _,nz_1=rhov.shape
    #vec=np.linspace(0,1-self.w0_imobiles,100)
    # CNT DGL
    R=8.31445
    NA=6.023E23
    kB=R/NA
    temp=298.15
    AR=100
    pre=rho*np.pi/(4*AR**2)
    scale_r=1E-6
    C0=rho/M*NA
    X_la=wi0[crystallizes]-alpha
    Xn_la=X_la/M
    dl_la = (wi0[crystallizes]-alpha)/(1-alpha)
    X2II=rhov/rho
    w2II=X2II/(X2II+1.)
    c_la=dl_la*(1-w2II)
    # c_la_0=wi0[mobiles]*(1-rhov[0,0]/rho)
    # c_GGW=0.
    # ideal=0
    lnaSLE=deltaHSL/R/temp*(1/TSL-1/temp)+cpSL*np.log(1-temp/TSL)
    lnai=lngi_fun(w2II)-np.log(w2II)
    lnS_fun=lnai[crystallizes]-lnaSLE

    # dmu_sla0=np.log(c_la_0/c_GGW) if ideal else lnS_fun(w2II)
    dmu_sla0=lnS_fun(wi0)
    
    dmu_sla=np.ones(nz_1)
    for i in range(nz_1):
        dmu_sla[i]=lnS_fun(w2II)
        # dmu_sla[i]=np.log(c_la[i]/c_GGW) if ideal else lnS_fun(w2II)
    rstar = ((2*sigma)/(C0*kB*temp*dmu_sla))
    
    #r0=np.ones(self.nz+1)/1E12
    r0=np.ones(nz_1)*2*sigma/(C0*dmu_sla0*kB*temp)*scale_r #taking r_star might solve the numerical problems in the volume element before the surface. Critical radius
    deltaG=sigma**3*(16*np.pi)/(3*(C0*kB*temp*dmu_sla)**2)
    beta=np.fmin((rhov-rhov[0,0])/(rhov[-1,0]-rhov[0,0]),1)
    #tau=(self.L0/2)**2/self.D12
    #beta=1-cs.exp(-Time/tau)
    ze=(kB*temp/sigma)**(1.5)*C0/(8*np.pi)*dmu_sla**2
    f=4*np.pi*rstar*DAPI*Xn_la*NA*np.fmax(beta,1E-4)
    
    dNdt = ze*f*C0*np.exp(-deltaG/(kB*temp))#*cs.exp(-NA)/NA**0.5
    #dNdt =NA*C0*cs.fmax(beta,1E-4)*DAPI*(kB*temp/sigma)**0.5*dmu_sla*Xn_la*cs.exp(-16/3*np.pi/dmu_sla**2/C0**2*(sigma/kB/temp)**3)
    drdt = np.fmax(1/scale_r*kt*np.fmax(beta,1E-4)*(dmu_sla)**g,0)
    dalphadr=3*alpha/r
    dalphadN=pre*(r*scale_r)**3
    #dalphadN=rho*np.pi/(4*AR**2)*(r*scale_r)**3
    dalphadt=np.fmax((drdt*dalphadr+dNdt*dalphadN),0)


    rhovB=ww_fun(dl_la)

    #paper=np.pi/(4*AR**2)*(r*scale_r)**3*C0**2*cs.fmax(beta,1E-4)*DAPI*(kB*temp/sigma)**0.5*dmu_sla*X_la*cs.exp(-16/3*np.pi/dmu_sla**2/C0**2*(sigma/kB/temp)**3)

    #dalphadt=cs.fmax((drdt*dalphadr+paper),0)
    return 