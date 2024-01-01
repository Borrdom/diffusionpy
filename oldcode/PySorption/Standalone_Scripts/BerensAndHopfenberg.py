# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 14:26:56 2021

@author: domin
"""
import numpy as np
from scipy.optimize import least_squares
from scipy.ndimage import median_filter
def Diffusion(t,kf,mfinfty):
    ninfty=100
    mtdiff=mfinfty*(1-sum([8/np.pi**2*1/(2*n+1)**2*np.exp(-(2*n+1)**2*kf*t) for n in range(ninfty)]))
    return mtdiff

def Relaxation(t,kr,nr,mrinfty):
    return sum([(1-np.exp(-kr[n]*t))*mrinfty[n] for n in range(nr)])
def BHModel(t,kf,kr,nr,mfinfty,mrinfty):
    return Diffusion(t,kf,mfinfty)+Relaxation(t,kr,nr,mrinfty)
def BHX(t,kf,kr,nr,X0,Xinfty,mfinfty,mrinfty):
    M_M=BHModel(t,kf,kr,nr,mfinfty,mrinfty)
    return (1-M_M)*X0+M_M*Xinfty
def XtoW(X):
    return X/(1+X)
def WtoX(W):
    return W/(1-W)

def InterpretExperiment(texp,wexp,nr,wbool=True):
    X=WtoX(wexp) if wbool else wexp
    X0=X[0]
    Xinfty=X[-1]

    X_filt=median_filter(X,10,mode="nearest")
    dXdt=np.gradient(X_filt,texp**(1/2))
    plt.plot(dXdt)
    plt.plot()
    Xfinfty=X[np.where(dXdt<1E-4)[0][0]]
    mfinfty=(Xfinfty-X0)/(Xinfty-X0)
    mrinfty=np.ones(nr)*(1-mfinfty)/nr 
    return X0,Xinfty,mfinfty,mrinfty

    
    
    
if __name__=="__main__":
    import matplotlib.pyplot as plt
    from tkinter.filedialog import askopenfilename
    import pandas as pd
    filename='C:/Users/domin/OneDrive/Promotion/ProgrammCode/Diffusionsmodell/Python_Boettcher/BHModel.xlsx'
    filename=askopenfilename() if not ('filename' in locals()) else filename
    dF=pd.read_excel(filename)
    texp1,wexp1=dF["t1"].values*60,dF["w1"].values
    texp2,wexp2=dF["t2"].values*60,dF["w2"].values
    nr=2
    X01,Xinfty1,mfinfty1,mrinfty1=InterpretExperiment(texp1,wexp1,nr,wbool=True)
    X02,Xinfty2,mfinfty2,mrinfty2=InterpretExperiment(texp2,wexp2,nr,wbool=True)
    fig1,ax1=plt.subplots()
    fig2,ax2=plt.subplots()

    ax1.set_xlabel("$t^{0.5}$/min")
    ax1.set_ylabel("$X_w$/-")
    ax2.set_xlabel("$t^{0.5}$/min")
    ax2.set_ylabel("$X_w$/-")
    Xexp1,Xexp2=WtoX(wexp1),WtoX(wexp2)
    ax1.plot(((texp1)**(1/2))/60,Xexp1,"kx")
    ax2.plot(((texp2)**(1/2))/60,Xexp2,"kx")
    D=1E-12
    l=7E-6
    kf=np.pi**2*D/l**2
    kr=np.ones(nr)*2E-4
    kr[0]=kr[0]/5
    X1=BHX(texp1,kf,kr,nr,X01,Xinfty1,mfinfty1,mrinfty1)
    X2=BHX(texp2,kf,kr,nr,X02,Xinfty2,mfinfty1,mrinfty1)
    ax1.plot(((texp1)**(1/2))/60,X1,"k-")
    ax2.plot(((texp2)**(1/2))/60,X2,"k-")
    w1,w2=XtoW(X1),XtoW(X2)
    #ax2.plot((t2**(1/2))/60,mt2)
    # wt1=mt1/(1+mt1)
    # wt2=mt2/(1+mt2)
    # nr=1
    # kr=np.ones(nr)*9E-5
    # #mrinfty=np.ones(nr)*0.8
    # #mfinfty=1

    # mrinfty=np.ones(nr)*0.5807
    # mrinfty=mrinfty/((1-mrinfty))
    # mfinfty=0.50271
    # mfinfty=mfinfty/((1-mfinfty))
    # m0=0.3131
    # m0=m0/((1-m0))
    # percentage=(mfinfty-m0)/(np.max(mrinfty)-m0)
    # D=4E-11
    # l=7E-6
    # t0=0
    # tend=1400*60
    # nt=300

    # t1=(np.linspace(t0,tend**(1/2),nt))**2
    # tmin1=t1/60
    # kf=np.pi**2*D/l**2
    # mt1=BHModel(t1,kf,kr,nr,mfinfty,mrinfty,m0)
    
    # m0=0.15269
    # m0=m0/((1-m0))
    # tend=160*60
    # t2=(np.linspace(t0,tend**(1/2),nt))**2
    # tmin2=t2/60
    # mrinfty=np.ones(nr)*0.28536
    # mrinfty=mrinfty/((1-mrinfty))
    # mfinfty=percentage*(np.max(mrinfty)-m0)+m0
    # mt2=BHModel(t2,kf,kr,nr,mfinfty,mrinfty,m0)
    
    
    # fig1,ax1=plt.subplots()
    # fig2,ax2=plt.subplots()
    # ax1.plot((t1**(1/2))/60,mt1)
    # ax1.set_xlabel("$t^{0.5}$/min")
    # ax1.set_ylabel("$X_w$/-")
    # ax2.set_xlabel("$t^{0.5}$/min")
    # ax2.set_ylabel("$X_w$/-")
    # import pandas as pd
    # from tkinter.filedialog import askopenfilename
    # filename='C:/Users/domin/OneDrive/Promotion/ProgrammCode/Diffusionsmodell/Python_Boettcher/BHModel.xlsx'
    # filename=askopenfilename() if not ('filename' in locals()) else filename
    # dF=pd.read_excel(filename)
    
    # texp1,wexp1=dF["t1"].values,dF["w1"].values
    # texp2,wexp2=dF["t2"].values,dF["w2"].values
    # mexp1=wexp1/(1-wexp1)
    # mexp2=wexp2/(1-wexp2)
    # ax1.plot(((texp1*60)**(1/2))/60,mexp1,"kx")
    # ax2.plot(((texp2*60)**(1/2))/60,mexp2,"kx")
    # ax2.plot((t2**(1/2))/60,mt2)
    # wt1=mt1/(1+mt1)
    # wt2=mt2/(1+mt2)
