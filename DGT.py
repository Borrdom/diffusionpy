import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

Mi=np.asarray([18.015,92.14])
ri=Mi/np.min(Mi)

def deltaG(phi,ri,betaij,gammaij):
    deltaS=np.sum(phi/ri*np.log(phi))
    nc=len(ri)
    deltaH=0
    for i in range(nc):
        for j in range(nc):
            if j>i: deltaH+=phi[i]*phi[j]/(1.-gammaij*phi[j])*betaij
    return deltaS+deltaH
betaij=2.6
gammaij=0.59

L=50/1000000
nz=20
nt=40
dz=L/nz
L12=1.4E-16/dz**2
kbin=-8.8E-7/dz**2
def mu(phi,ri,betaij,gammaij):
    h=1E-21
    dphi=1j*h
    phi_=np.zeros_like(phi).astype(complex)
    mu_=np.zeros_like(phi)
    nc=len(phi)
    for i in range(nc):
        phi_=np.zeros_like(phi).astype(complex)
        phi_+=phi
        phi_[i]+=dphi
        phi_[i!=np.arange(0,nc)]-=dphi
        mu_[i]=deltaG(phi_,ri,betaij,gammaij).imag/h
    return mu_

phi1=np.linspace(0,1,1000)
mu1=np.zeros_like(phi1)
mu2=np.zeros_like(phi1)
for i in range(1000):
    phi_=np.hstack((phi1[i],1-phi1[i]))
    mui=mu(phi_,ri,betaij,gammaij)
    mu1[i]=mui[0]
    mu2[i]=mui[1]
# plt.plot(phi1,mu1)
# plt.plot(phi1,mu2)
delg=phi1*mu1/ri[0]+(1-phi1)*mu2/ri[1]#-phi1/ri[0]*np.log(phi1)-(1-phi1)/ri[1]*np.log(1-phi1)
# plt.plot(phi1,delg)
# plt.show()

phi=np.asarray([0.5,0.5])
mu1=mu(phi,ri,betaij,gammaij)
print(mu1)

def dphidt(t,phi1):
    deltaphi1=np.diff(phi1)
    deltaphi1=np.hstack((0,deltaphi1))
    delta2phi1=np.diff(deltaphi1)
    delta2phi1=np.hstack((delta2phi1,0))
    mu1=np.zeros_like(phi1)
    mu2=np.zeros_like(phi1)
    r1=ri[0]
    r2=ri[1]
    for i in range(nz):
        phi_=np.hstack((phi1[i],1-phi1[i]))
        mui=mu(phi_,ri,betaij,gammaij)
        mu1[i]=mui[0]
        mu2[i]=mui[1]
    deltamu=np.diff((mu1/r1-mu2/r2+kbin*delta2phi1))
    phi1bar=(phi1[1:]+phi1[:-1])/2
    dphi1dt=L12*phi1bar*(1-phi1bar)*deltamu
    dphi1dt=np.hstack((dphi1dt,0))
    dphi1dt[phi1<0.01]=0
    return dphi1dt
phi10=np.ones((nz))/10
phi10[-1]=0.9
t=np.linspace(0,100,nt)
z=np.linspace(0,L,nz)
sol=solve_ivp(dphidt,(t[0],t[-1]),phi10,method="Radau",t_eval=t)

print(sol)


plt.plot(z,sol["y"])
plt.show()