import numpy as np
from numba import njit
from scipy.interpolate import interp1d
from.PyCSAFT_nue import lngi,vpure
from.Stefan_Maxwell_segmental import D_Matrix
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


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

kij=D_Matrix(np.asarray([-0.156,0.00648,-0.0574]),nc)
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
nt=300
t=np.linspace(0,3000,nt)*60

lngi_fun=lambda wi: lngi(T,wi,**par)




def Bound(t,mobiles,immobiles,crystallize,wi0,wi8,rho0i,Mi,DAPI,sigma,kt,g,deltaHSL,TSL,cpSL,lngi_fun=None,wv_fun=None):
    crystallizes=np.where(crystallize)[0]
    alpha0=0.01
    r0=0.1
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
    crystallizes=np.where(crystallize)[0]
    R=8.31445
    NA=6.023E23
    kB=R/NA
    temp=298.15
    AR=100
    M=Mi[crystallizes]/1000.
    rho=rho0i[crystallizes]
    pre=rho*np.pi/(4*AR**2)
    scale_r=1E-6
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
    # plt.plot(t,np.exp(dmu_sla),"kx")
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

wiB,alpha,r=Bound(t,mobiles,immobiles,crystallize,wi0,wi8,rho0i,Mi,DAPI,sigma,kt,g,deltaHSL,TSL,cpSL,lngi_fun,wv_fun)
import matplotlib.pyplot as plt


# plt.plot(t/60,alpha)
# plt.plot(t/60,r)
# plt.plot(t/60,wiB[0,:])
# plt.plot(t/60,wiB[1,:])
# plt.plot(t/60,wiB[2,:])
plt.show()

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