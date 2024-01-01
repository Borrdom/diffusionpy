import casadi as cs
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
R=8.3145
T=298.15
rho02=997
rho01=1180
M2=0.01815
RV=R*T*(rho02/M2)
des=False
activity2=0.1
activity1=0.001
wGGW1=0.001
wGGW2=0.1
XGGW1=wGGW1/(1-wGGW1)
XGGW2=wGGW2/(1-wGGW2)
rho1GGW=XGGW1*rho01
rho2GGW=XGGW2*rho02
global ax,ax1,ax2,ax3
fig,ax=plt.subplots()
fig,ax1=plt.subplots()
fig,ax2=plt.subplots()
fig,ax3=plt.subplots()


def TestDamper(activity1,activity2,rho1GGW,rho2GGW,RV, des):
    a=activity2 if not des else activity1
    b=activity1 if not des else activity2
    c=rho1GGW if not des else rho2GGW
    d=rho2GGW if not des else rho1GGW
    A=RV
    B=1E12
    tau=B
    p=1
    m=(a-b)/(d-c)
    act=lambda rho: m*(rho-c)+b
    m=(a-b)/(d**p-c**p)
    ord=b-m*c**p
    act=lambda rho: m*rho**p+ord
    rhospan=np.linspace(rho1GGW,rho2GGW,100)


    ax2.plot(rhospan,act(rhospan))
    ax2.set_ylabel(r"$a_s/-$")
    ax2.set_xlabel(r"$\rho_s/kg/m^3$")
    #print(m)
    #act=lambda rho: b*(rho/d)**m
    #print(act(c))
    #print(act(d))
    #act=lambda rho: np.log(rho*c2)*c1
    t=np.linspace(0,10000**(1/2),100)**2
    rhs=lambda rho,t:rho02*A/B*np.log(a/act(rho))

    rhow=odeint(rhs, c, t)
    Q=-(rhow-rhow[0])/(rhow[0]-rhow[-1])
    #rhow=rhow if not des else -rhow+rhow[0]+rhow[-1]
    ax1.plot(np.sqrt(t/60),rhow)
    ax1.set_xlabel(r"$t^{0.5}/s^{0.5}$")
    ax1.set_ylabel(r"$\rho_s/kg/m^3$")
    #ax1.plot(np.sqrt(t/60),Q)

    ax.plot(rhospan,rhs(rhospan,t))
    ax.set_xlabel(r"$\rho_s/kg/m^3$")
    ax.set_ylabel(r"$\frac{d\rho_s}{dt}/kg/m^3/s$")

    drhsdrho=lambda rho:-rho02*A/B*a/act(rho)*m
    ax3.plot(rhospan,drhsdrho(rhospan))
    ax3.set_ylabel(r"$\frac{d\rho_s}{dtd\rho_w}/s$")
    ax3.set_xlabel(r"$\rho_s/kg/m^3$")
TestDamper(activity1,activity2,rho1GGW,rho2GGW,RV, True)
TestDamper(activity1,activity2,rho1GGW,rho2GGW,RV, False)
#etavec=10**np.linspace(3,13)
#tauvecDamper=etavec/RV
#E=80775.63056
#E=7.54134E13
#D=1E-13
#L=1E-5
#tauDiff=(2*L)**2/D
#tauDiff=360
#De=10
#taudeltat=tauDiff/De
#tauvec=etavec/E
#plt.plot(np.log10(etavec),np.log10(tauvecDamper),'k-')
#plt.plot(np.log10(etavec),np.log10(tauvec),'k--')
#plt.plot(np.log10(etavec),np.ones_like(etavec)*np.log10(tauDiff),'b-')
#plt.plot(np.log10(etavec),np.ones_like(etavec)*np.log10(taudeltat),'r-')
plt.show()
