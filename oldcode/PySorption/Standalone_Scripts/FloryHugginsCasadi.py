import casadi as cs
import numpy as np
from PCSAFTPolynomial import k,Nav
import matplotlib.pyplot as plt
def Flory(w1,chi12):
    nc=2
    xi=cs.SX.sym("xi",nc)
    R=k*Nav
    
    T=298.15
    Mi=cs.SX.ones(nc)
    Mi[0]=18.02
    Mi[1]=357.79
    rho0i=cs.SX.ones(nc)
    rho0i[0]=997
    rho0i[1]=1320
    phii=xi*Mi/rho0i/cs.sum1(xi*Mi/rho0i)
    wi=xi*Mi/cs.sum1(xi*Mi)
    deltagE=R*T*(xi[0]*cs.log(phii[1])+xi[1]*cs.log(phii[0])+chi12*xi[0]*phii[1])
    gammai=cs.exp(cs.gradient(deltagE,xi))
    gammai=cs.exp(cs.log(phii/xi)+1-phii/xi+chi12*(1-phii)**2)
    gammai_fun=cs.Function("gammai_fun",[xi],[gammai])
    wi_fun=cs.Function("gammai_fun",[xi],[wi])
    n=1000
    
    x1vec=cs.linspace(0,1,n)
    x2vec=1-x1vec
    
    gammaivec=np.asarray([gammai_fun(cs.vertcat(x1vec[i],x2vec[i])).full().T[0] for i in range(n)])
    wivec=np.asarray([wi_fun(cs.vertcat(x1vec[i],x2vec[i])).full().T[0] for i in range(n)])
    
    gamma1vec=gammaivec[:,0]
    gamma2vec=gammaivec[:,1]
    w1vec=wivec[:,0]
    w2vec=wivec[:,1]
    RH=np.nan_to_num(gamma1vec*x1vec,0)
    from scipy.interpolate import InterpolatedUnivariateSpline
    
    RH_fun=InterpolatedUnivariateSpline(w1vec,RH)
    return RH_fun(w1)
from scipy.optimize import curve_fit
chi120=2.5
w10=np.asarray([0.01,0.015,0.02])
RH=Flory(w10,chi120)
fig,ax=plt.subplots()
#ax.plot(RH,w1vec)
#ax.plot(x1vec,gamma2vec)
import pandas as pd
from os import getcwd
from os.path import join
import numpy as np

cwd=getcwd()
filename=join(cwd,"LeastFit","indometacin.csv")
Lit=pd.read_csv(filename).dropna() 
RHbar=Lit["RHbar[-]"].values
wH2Obar=Lit["wH20bar[-]"].values
ax.plot(RHbar,wH2Obar,'kx')
plt.xlim([0,1])
plt.ylim([0,0.03])
n=100
w1vec=cs.linspace(0,0.03,n)
popt, pcov = curve_fit(Flory, wH2Obar, RHbar,p0=chi120)
RHopt=Flory(w1vec,popt)
plt.plot(RHopt,w1vec)
