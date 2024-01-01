import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline

Tg0=np.asarray([441.51,136,317.6])#PVP #303.6 Indometacin to hit right transition

#Tg0=np.asarray([383.9,136,323.15]) # PVPVARitonavir
#Tg0=np.asarray([383.9,136,317.6])#PVPVA
rho0=np.asarray([1210,997,1320])
#rho0=np.asarray([1180,997,1320])#PVPVA64 indo
#rho0=np.asarray([1270,997,1150]) # Changed PVPVA density to fit Birtes data point
#rho0=np.asarray([1180,997,1150])
#rho0=np.asarray([1380,997,1150]) #Changed to fit data points from paper
#TG(wF3,Tg0,rho0)



def TG(wF3,Tg0=Tg0,rho0=rho0):
    #Tg0=np.asarray([441.51,136,317.6])
    #Tg0=np.asarray([383.9,136,317.6])
    #rho0=np.asarray([1210,997,1320])
    #rho0=np.asarray([1190,997,1320])
    k1=Tg0[1]*rho0[1]/(Tg0[2]*rho0[2])
    k2=Tg0[1]*rho0[1]/(Tg0[0]*rho0[0])
    w2=np.linspace(0,1,100)
    w3=lambda wF3:(1-w2)*(wF3)
    w1=lambda wF3:1-w2-w3(wF3)
    Tg=lambda wF3:(w2*Tg0[1]+k1*w3(wF3)*Tg0[2]+k2*w1(wF3)*Tg0[0])/(w2+k1*w3(wF3)+k2*w1(wF3))
    T=273.15+25
    #Gordon Taylor Gleichung nach Wassermasse umgestellt
    wg=(k1*wF3*(Tg0[2]-T)+k2*(1-wF3)*(Tg0[0]-T))/((T-Tg0[1])+wF3*k1*(Tg0[2]-T)+(1-wF3)*k2*(Tg0[0]-T))
    #
    #fig,ax=plt.subplots()
    Tgs=[]
    wgs=[]
    for i,val in enumerate(wF3):
        #ax.plot(w2,Tg(val)-273.15)
        Tgs.append(InterpolatedUnivariateSpline(w2,Tg(val)-273.15-25))
        #ax.plot(w2,Tgs[i](w2)+25,'--')
        wgs.append(Tgs[i].roots())
        #print(Tg(val)-273.15)
    wgs=np.asarray(wgs)
    #ax.plot([0,1],[25,25])
    #ax.plot(wgs,np.ones_like(wgs)*25,'kx')
    return wgs,wg
def TG2(w2,wF3,Tg0=Tg0,rho0=rho0):

    k1=Tg0[1]*rho0[1]/(Tg0[2]*rho0[2])
    k2=Tg0[1]*rho0[1]/(Tg0[0]*rho0[0])
    #w2=np.linspace(0,1,100)
    w3=lambda wF3:(1-w2)*(wF3)
    w1=lambda wF3:1-w2-w3(wF3)
    Tg=lambda wF3:(w2*Tg0[1]+k1*w3(wF3)*Tg0[2]+k2*w1(wF3)*Tg0[0])/(w2+k1*w3(wF3)+k2*w1(wF3))
    return Tg(wF3)
w3F=np.linspace(0,1,11)

import pandas as pd
wg=TG(w3F,Tg0,rho0)[1]
data=np.vstack(((1-wg)*(w3F),wg,(1-wg)*(1-w3F))).T
pd.DataFrame(data).to_excel("Tgs.xlsx",header=["wa[-]","ww[-]","wp[-]"])
#wg,wgs=TG(w3F)
#w2=np.linspace(0,.4,100)
#Tgs=np.asarray([(TG2(val,w3F)-273.15)[0] for i,val in enumerate(w2)])
#print(TG(w3F))