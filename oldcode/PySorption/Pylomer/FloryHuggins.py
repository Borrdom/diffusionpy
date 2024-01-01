import casadi as cs
import numpy as np
from scipy.optimize import curve_fit
import pandas as pd
from os import getcwd
from os.path import join

import matplotlib.pyplot as plt
# from playsound import playsound
# playsound('_CLICK_Nice.mp3')
from scipy.interpolate import InterpolatedUnivariateSpline
from PyCSAFT import get_pcpar

R,Nav=8.31448,6.022E23
k=R/Nav
class FloryFilm:

    def __init__(self,matrix,solvent,wmax=1):



        self.Matrix,self.solvent=matrix,solvent
        self.wmax=wmax
        self.FitToCurve()


    def FloryRH(self,w2,chi12):
        nc=2
        xi=cs.SX.sym("xi",nc)
        R=k*Nav
        T=298.15

        Mi=cs.SX.ones(nc)
        Voi=cs.SX.ones(nc)
        Mi[1]=self.solvent["Mi"]
        Mi[0]=self.Matrix["Mi"]
        rho0i=cs.SX.ones(nc)
        rho0i[1]=self.solvent["rho0"]
        rho0i[0]=self.Matrix["rho0"]
        Voi[1]=self.solvent["Mi"]/self.solvent["rho0"]
        Voi[0]=self.Matrix["Mi"]/self.Matrix["rho0"]

        phii=xi*Mi/rho0i/cs.sum1(xi*Mi/rho0i)
        wi=xi*Mi/cs.sum1(xi*Mi)
        gammai=cs.exp(cs.log(phii/xi)+1-phii/xi+chi12*(1-phii)**2) #might be wrong 05.01.2021
        gammai=cs.exp(cs.log(phii/xi)+(1-phii)*(1-Voi[1]/+Voi[0])+chi12*(1-phii)**2) # makes "almost" entirely no difference even in the API case. Its identical
        gammai_fun=cs.Function("TH_fun",[xi],[gammai])
        gammai=gammai_fun(cs.vertcat(1-xi[1],xi[1]))

        #Added 22.05.20 In order caluclate critical Point of LLE
        gammai2=cs.exp(cs.log(phii/xi)+1-phii/xi+chi12*(1-phii)**2)
        rhoi=cs.SX.sym("rhoi")
        rho=cs.SX.sym("rho")
        rho0imol=rho0i/(Mi)
        rhoeq=cs.sum1(xi/rho0imol)**-1
        muid=R*T*cs.log(xi*rho*R*T)+R*T
        #muid=0
        mui=muid+R*T*cs.log(xi*gammai2)
        mui2=R*T*cs.log(xi*rhoeq*R*T)+R*T+R*T*cs.log(xi*gammai2)
        mui_fun=cs.Function("mui_fun",[xi],[mui2])
        G=cs.sum1(mui*xi)

        xieq=rhoi/rho

        G_fun=cs.Function("G_fun",[xi],[G])
        G=G_fun(xieq)
        dGdrhoi=cs.gradient(G*rho,rhoi)
        rhoieq=xi*rho
        dGdrhoi_fun=cs.Function("dGdrhoi_fun",[rhoi],[dGdrhoi])
        dGdrhoi=dGdrhoi_fun(rhoieq)
        dGdrhoi_fun=cs.Function("dGdrhoi_fun",[rho],[dGdrhoi])
        dGdrhoi=dGdrhoi_fun(rhoeq)
        dGdrhoi_fun=cs.Function("dGdrhoi_fun",[xi],[dGdrhoi])



        THFaktor=1+xi[1]*cs.gradient(cs.log(gammai)[1],xi[1])
        THFaktor_fun=cs.Function("TH_fun",[xi],[THFaktor])
        gammai_fun=cs.Function("gammai_fun",[xi],[gammai])
        wi_fun=cs.Function("gammai_fun",[xi],[wi])
        n=1000
        x2vec=cs.linspace(0,1,n)
        x1vec=1-x2vec
        gammaivec=np.asarray([gammai_fun(cs.vertcat(x1vec[i],x2vec[i])).full().T[0] for i in range(n)])
        wivec=np.asarray([wi_fun(cs.vertcat(x1vec[i],x2vec[i])).full().T[0] for i in range(n)])
        THFaktorvec=np.asarray([THFaktor_fun(cs.vertcat(x1vec[i],x2vec[i])).full().T[0] for i in range(n)])
        dGdrhoivec=np.asarray([dGdrhoi_fun(cs.vertcat(x1vec[i],x2vec[i])).full().T[0] for i in range(n)])
        muivec=np.asarray([mui_fun(cs.vertcat(x1vec[i],x2vec[i])).full().T[0] for i in range(n)])

        #fig,ax=plt.subplots()

        gamma1vec=gammaivec[:,0]
        gamma2vec=gammaivec[:,1]
        w1vec=wivec[:,0]
        w2vec=wivec[:,1]
        self.THFaktor_fun=InterpolatedUnivariateSpline(w2vec,np.nan_to_num(THFaktorvec,0))
        self.gamma2_fun=InterpolatedUnivariateSpline(w2vec,np.nan_to_num(gamma2vec,0))
        RH=np.nan_to_num(gamma2vec*x2vec,0)
        #x2bar=wH2Obar/FilmFlory.solvent["Mi"]/(wH2Obar/FilmFlory.solvent["Mi"]+(1-wH2Obar)/FilmFlory.Matrix["Mi"])
        #ax.plot(w2vec,dGdrhoivec)
        #ax.plot(w2vec,muivec)
        self.w2vec=w2vec
        self.gamma2vec=gamma2vec
        RH_fun=InterpolatedUnivariateSpline(w2vec,RH)
        return RH_fun(w2)

    def GetIsothermData(self,molecule):


        cwd=getcwd()
        filename=join(cwd,"LeastFit",molecule+".csv")
        Lit=pd.read_csv(filename)
        RHbar=Lit["RHbar[-]"].dropna().values
        wH2Obar=Lit["wH20bar[-]"].dropna().values
        self.RHbar=RHbar
        self.wH2Obar=wH2Obar
        return RHbar,wH2Obar

    def FitToCurve(self):

        RHbar,wH2Obar=self.GetIsothermData(self.Matrix["name"])
        chi120=2.5
        w20=np.asarray([0.01,0.015,0.02])
        RH=self.FloryRH(w20,chi120)

        #ax.plot(RH,w1vec)
        #ax.plot(x1vec,gamma2vec)

        n=1000
        w2vec=np.linspace(0,self.wmax,n)

        popt, pcov = curve_fit(self.FloryRH, wH2Obar, RHbar,p0=chi120)
        print(popt)

        RHopt=self.FloryRH(w2vec,popt)
        self.THFaktorvec=self.THFaktor_fun(w2vec)
        self.gamma2vec=self.gamma2_fun(w2vec)
        self.RHopt,self.w2,self.chi12=RHopt,w2vec,popt



if __name__=="__main__":
    T=298.15
    pure, kij=get_pcpar.get_par(["indomethacin","water","indomethacin"],T=T)
    # pure, kij=get_pcpar.get_par(["pvpva64","water","nifedipin"],T=T)
    pol,water,api=pure
    FilmFlory=FloryFilm(pol,water,wmax=1)
    fig,ax=plt.subplots()
    wH2Obar=FilmFlory.wH2Obar
    ax.plot(FilmFlory.RHbar,FilmFlory.wH2Obar,'kx')
    x2bar=wH2Obar/FilmFlory.solvent["Mi"]/(wH2Obar/FilmFlory.solvent["Mi"]+(1-wH2Obar)/FilmFlory.Matrix["Mi"])
    plt.xlim([0,1])
    plt.ylim([0,0.03])
    ax.plot(FilmFlory.RHopt,FilmFlory.w2)
    fig1,ax1=plt.subplots()
    ax1.plot(FilmFlory.w2,FilmFlory.THFaktorvec)
    fig2,ax2=plt.subplots()
    ax2.plot(FilmFlory.RHbar,x2bar,'kx')
