import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def ReadFittings(sheetnumber,component):
    Jumps=[0,1,2,3,4,5]

    ending="_fitting.xlsx"
    values=[]
    eta0=[]
    E0=[]
    tau0=[]
    D12=[]
    Tginit=[]
    TgGGW=[]
    Tgaverage=[]
    wGGW=[]
    winit=[]
    
    for i,val in enumerate(Jumps):
        current=str(val)+component+ending
        data=(pd.read_excel(current,sheet_name=sheetnumber))
        values.append(dict(zip(data["name"],data["value"])))
    
        tau0.append(values[i]["tau00"])
        E0.append(values[i]["E00"])
        D12.append(values[i]["D12"])
        eta0.append(values[i]["eta00"])
        Tginit.append(values[i]["Tginit"])
        TgGGW.append(values[i]["TgGGW"])
        Tgaverage.append(values[i]["Tgaverage"])
        wGGW.append(values[i]["wGGW"])
        winit.append(values[i]["winit"])
    fig,ax=plt.subplots()
    fig1,ax1=plt.subplots()
    etawater=0.001#
    Tginit=np.asarray(Tginit)
    Tgaverage=np.asarray(Tgaverage)
    TgGGW=np.asarray(TgGGW)
    eta0=np.asarray(eta0)
    #Terr=np.vstack((-Tgaverage+Tginit,Tgaverage-TgGGW))
    ax.plot(Tginit,np.log10(np.asarray(eta0)/3),'bx')
    ax.plot(Tgaverage,np.log10(np.asarray(eta0)/3),'gx')
    ax.plot(TgGGW,np.log10(np.asarray(eta0)/3),'rx')
    #ax.errorbar(xerr=Tgaverage-Tginit,x=Tgaverage,y=np.log10(np.asarray(eta0)*(1-np.asarray(wGGW))**-2),fmt='bx')
    eta0=np.asarray(eta0)
    wGGW=np.asarray(wGGW)
    fig2,ax2=plt.subplots()
    #ax2.plot(Tgaverage,np.log10((eta0/etawater**(wGGW))**(1/(1-wGGW))),'bx') 
    ax2.plot(Tgaverage,np.log10(D12),'kx')
    ax1.plot(winit,np.log10(np.asarray(eta0)),'bx')
    ax1.plot(wGGW,np.log10(np.asarray(eta0)),'rx')
    poly=np.polyfit((winit[3:]+0*wGGW[3:]),np.log(eta0)[3:],1)
    polyval=np.polyval(poly,wGGW)*np.log10(np.e)
    print(poly)
    ax1.plot(wGGW,polyval)
    ax.set_ylim([0,20])
    
    # -*- coding: utf-8 -*-
    """
    Created on Sun Oct 25 18:43:15 2020
    
    @author: domin
    """
    from TgGordnonTaylor import TG2
    
    
    Tg0=np.asarray([383.9,136,323.15]) # PVPVARitonavir
    #Tg0=np.asarray([383.9,136,317.6])#PVPVA
    #rho0=np.asarray([1210,997,1320])
    rho0=np.asarray([1180,997,1320])#PVPVA64 indo
    C1=4.68
    C2=25.2
    Tref=110+273.15
    etaref=8.42E8
    Tvec=np.linspace(273.15,440,100)
    def WLF(T,Tref,C1,C2):
        return 10**(-C1*(T-Tref)/(C2+T-Tref))*etaref
    
    w2=np.linspace(0.1,0.2,100)
    
    Tg=TG2(w2,[0],Tg0,rho0)
    T=298.15
    Tapp=T+Tg[0]-Tg
    
    Tvec=np.linspace(370,440,100)
    Shif=WLF(Tvec,Tref,C1,C2)

    ax.plot(Tvec,np.log10(Shif))
    #ax.plot([Tref,Tref],[0,20])

    
ReadFittings(0,"pvpva64")
ReadFittings(1,"pvpva64")
ReadFittings(2,"pvpva64")
ReadFittings(3,"pvpva64")
#ReadFittings(0,"pvp")
#ReadFittings(1,"pvp")
#ReadFittings(2,"pvp")