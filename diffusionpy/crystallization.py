import numpy as np
from numba import njit
from scipy import InterpolatedUnivariateSpline

def Kris(rho2II,alpha,r,wwSIMraw,wAPISIMraw,rhoSIMraw,w0_imobiles,rho):
    # WasserFeed Information
    XwSIMraw=wwSIMraw/(1-wwSIMraw)
    rho2GGW_fun1=InterpolatedUnivariateSpline(wAPISIMraw[::-1], (rhoSIMraw*XwSIMraw)[::-1],k=1)
    rho2II_his[-1,0]=rho2GGW_fun1(1-w0_imobiles)
    #vec=np.linspace(0,1-self.w0_imobiles,100)
    # CNT DGL
    R=self.R
    NA=6.023E23
    kB=R/NA
    temp=self.T
    rho=self.rho03
    M=self.M3
    DAPI=self.DAPI
    g=self.g
    sigma=self.sigma
    kt=self.kt
    AR=100
    pre=rho*np.pi/(4*AR**2)
    scale_r=1E-6
    C0=rho/M*NA
    X_la=(1-self.w0_imobiles)-alpha
    
    Xn_la=X_la/M
    dl_la = (1-self.w0_imobiles-alpha)/(1-alpha)
    
    c_la=dl_la*(1-w2II)
    c_la_0=(1-self.w0_imobiles)*(1-rho2II_his[0,0]/rho)
    c_GGW=solubility if self.ideal else 0
    dmu_sla0=np.log(c_la_0/c_GGW) if self.ideal else deltamu_fun(rho2II_his[0,0]/rho,0)
    

    dmu_sla=cs.MX.ones(self.nz+1)
    for i in range(self.nz+1):
        dmu_sla[i]=np.log(c_la[i]/c_GGW) if self.ideal else deltamu_fun(w2II[i],alpha[i])
    rstar = ((2*sigma)/(C0*kB*temp*dmu_sla))
    
    #r0=np.ones(self.nz+1)/1E12
    r0=np.ones(self.nz+1)*2*sigma/(C0*dmu_sla0*kB*temp)*scale_r #taking r_star might solve the numerical problems in the volume element before the surface. Critical radius
    alpha0=pre*(r0/scale_r)**3/self.mnull #1 Nuceli per kg of dry mass
    deltaG=sigma**3*(16*np.pi)/(3*(C0*kB*temp*dmu_sla)**2)
    beta=cs.fmin((rho2II-rho2II_his[0,0])/(rho2II_his[-1,0]-rho2II_his[0,0]),1)
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
    #paper=np.pi/(4*AR**2)*(r*scale_r)**3*C0**2*cs.fmax(beta,1E-4)*DAPI*(kB*temp/sigma)**0.5*dmu_sla*X_la*cs.exp(-16/3*np.pi/dmu_sla**2/C0**2*(sigma/kB/temp)**3)

    #dalphadt=cs.fmax((drdt*dalphadr+paper),0)


    #BOUNDARY
    rho2GGW=rho2GGW_fun1(alphasym)
    return 