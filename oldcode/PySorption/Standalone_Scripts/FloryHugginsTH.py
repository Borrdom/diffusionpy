import numpy as np
from PCSAFTPolynomial import k,Nav
import matplotlib.pyplot as plt
def Flory(w1,chi12):
    R=k*Nav
    T=298.15
    from PureDataBase import indometacin,water
    solvent=water
    Matrix=indometacin
    M1=solvent["Mi"]
    M2=Matrix["Mi"]
    Mi=np.asarray([M1,M2])
    rho01=solvent["rho0"]
    rho02=Matrix["rho0"]
    rho0i=np.asarray([rho01,rho02])
    phii=lambda xi: xi*Mi/rho0i/np.sum(xi*Mi/rho0i)
    wi=lambda xi: xi*Mi/np.sum(xi*Mi)
    gammai= lambda xi: np.exp(np.log(phii(xi)/xi)+1-phii(xi)/xi+chi12*(1-phii(xi))**2)
    n=1000 
    x1vec=np.linspace(0,1,n)
    x2vec=1-x1vec
    xivec=np.vstack((x1vec,x2vec))
    gammaivec=np.asarray([gammai(xivec[:,i]) for i,val in enumerate(xivec[0,:])])
    wivec=np.asarray([wi(xivec[:,i]) for i,val in enumerate(xivec[0,:])])
    gamma1vec=gammaivec[:,0]
    gamma2vec=gammaivec[:,1]
    w1vec=wivec[:,0]
    w2vec=wivec[:,1]
    RH=np.nan_to_num(gamma1vec*x1vec,0)
    from scipy.interpolate import InterpolatedUnivariateSpline
    fig,ax=plt.subplots() 
    ax.plot(w1vec,1+x1vec*np.gradient(np.log(gamma1vec),x1vec))
    RH_fun=InterpolatedUnivariateSpline(w1vec,RH)
    return RH_fun(w1)
def GetIsothermData(molecule):
    import pandas as pd
    from os import getcwd
    from os.path import join
    import numpy as np
    
    cwd=getcwd()
    filename=join(cwd,"LeastFit",molecule+".csv")
    Lit=pd.read_csv(filename).dropna() 
    RHbar=Lit["RHbar[-]"].values
    wH2Obar=Lit["wH20bar[-]"].values
    return RHbar,wH2Obar

def FitToCurve():
    from scipy.optimize import curve_fit
    RHbar,wH2Obar=GetIsothermData("indometacin")
    chi120=2.5
    w10=np.asarray([0.01,0.015,0.02])
    RH=Flory(w10,chi120)
    fig,ax=plt.subplots()
    #ax.plot(RH,w1vec)
    ax.plot(x1vec,gamma2vec)
    ax.plot(RHbar,wH2Obar,'kx')
    plt.xlim([0,1])
    plt.ylim([0,0.03])
    n=100
    w1vec=np.linspace(0,0.03,n)
    popt, pcov = curve_fit(Flory, wH2Obar, RHbar,p0=chi120)
    RHopt=Flory(w1vec,popt)
    plt.plot(RHopt,w1vec)

if __name__=="__main__":
    FitToCurve()