
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
def langmuir(phiw,phiwinfty,k):
    return np.fmin(1/k*(np.fmin(phiw,phiwinfty)/(phiwinfty-np.fmin(phiw,phiwinfty))),1)
def flory(phiw,xhi):
    return np.log(phiw)+(1-phiw)+xhi*(1-phiw)**2
    #return np.log(phiw)+1-phiw/xi+xhi*(1-phiw)**2
    #return np.log(phiw)+(1-phiw)*(1-(0.0182/997)/(25.7/1216))+xhi*(1-phiw)**2
    #xi=phiw*997/18.2/(phiw*997/18.2+(1-phiw)*1216/25700)
    #return np.log(phiw)+1-phiw/xi+xhi*(1-phiw)**2
def henry(phiw,H):
    return phiw/H
def henry2(aw,H):
    return aw*H
def langmuir2(aw,phiwinfty,k):
    return phiwinfty*k*aw/(1+k*aw)
def flory2(aw,xhi):
    aw=np.atleast_1d(aw)
    floryaw=lambda aw,phi,xhi: np.log(aw)-flory(phi,xhi)
    phiwf=np.asarray([least_squares((lambda phi: floryaw(val,phi,xhi)),val)["x"][0] for i,val in enumerate(aw)])
    return phiwf

if __name__=="__main__":
    phiw=np.linspace(0,1,100)
    
    H=1
    k=10
    phiwinfty=0.3
    xhi=0.2
    awl=langmuir(phiw,phiwinfty,k)
    awh=henry(phiw,H)
    awf=np.exp(flory(phiw,xhi))
    fig,ax=plt.subplots()
    fig1,ax1=plt.subplots()
    aw=np.linspace(1E-10,1,1000)
    ax.plot(awl,phiw)
    ax.plot(awh,phiw)
    ax.plot(awf,phiw)
    #ax.plot(aw,henry2(aw,H)+langmuir2(aw,phiwinfty,k))
    ax.set_xlabel("RH/-")
    ax.set_ylabel("$\phi_w$/-")
    #log(aw)=log(phiw)+(1-phiw)+xhi*(1-phiw)**2
    
    phiwf=flory2(aw,xhi)
    phiwh=henry(aw,H)
    phiwl=langmuir2(aw,phiwinfty,k)
    ax.plot(aw,phiwf+phiwl,'k-')
    ax1.plot(flory(phiw,xhi),phiw)
