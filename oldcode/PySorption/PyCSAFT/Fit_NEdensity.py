from asyncio import WindowsProactorEventLoopPolicy
from PyCSAFT.PyCSAFTNil import Mixture
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import least_squares
from PyCSAFT import get_pcpar
import numpy as np
T=303.15
    #pure, kij=get_pcpar.get_par(["pvpva64","water","naproxen"],T=T)
pure, kij,ksw=get_pcpar.get_par(["pvpva64","water","indomethacin","ethanol"],T=T)

frame=pd.read_excel("data.xlsx")

yexp1=frame["w0"]
yexp2=frame["w1"]
yexp3=frame["w2"]
yexp4=frame["w3"]

yexp=np.asarray([yexp1,yexp2,yexp3,yexp4])


pol,water,api,solvent=pure
lij={}
ksw["waterpvpva64"]=0
ksw["pvpva64water"]=0
ksw["ethanolpvpva64"]=0
ksw["pvpva64ethanol"]=0
kij["ethanolpvpva64"]=-0.04
kij["pvpva64ethanol"]=-0.04
lij["ethanolpvpva64"]=0
lij["pvpva64ethanol"]=0



fig,ax=plt.subplots()


def Objective(x,yexp,p):
    Objective.counter += 1
    pol["rho0Poly0"]=p[0]*1000
    Film=Mixture(pol,water,api,solvent,dikij=kij,diksw=ksw,dilij=lij)
    Film.idx=np.asarray([1,3])
    sol,solnet=Film.Isotherm(wPolyASD=1,T=303.15,psys=1E5,NET=True,RSvec=x)
    resvec=[]
    wnet=np.asarray([solnet[i]["wi"] for i,val in enumerate(solnet)])
    for j,val in enumerate (yexp):
        resvec.append(((yexp[j]-wnet[:,j])/yexp[j])**2) if j!=2 else None #API not considered here
    dis=sum(resvec)**0.5
    mindis=np.min(dis)
    #ax.plot(Objective.counter,mindis,'kx')

    #plt.show(block=False)
    #plt.pause(0.001)
    return mindis
Objective.counter = 0

RHlin=np.asarray([0.3])
RHlin2=np.linspace(0.01,0.3,40)
RS1,RS2=np.meshgrid(RHlin,RHlin2)
RHRS=np.vstack((RS1.flatten(),RS2.flatten()))
RSvec=[]
rhooptvec=[]

for i,val in enumerate (yexp[0,:]):
    opt=least_squares(lambda p :Objective(RHRS,yexp[:,i],p),x0=1.200,bounds=(0.400,1.500))    
    rhoopt=opt["x"][0]  #1390
    pol["rho0Poly0"]=rhoopt*1000
    Film=Mixture(pol,water,api,solvent,dikij=kij,diksw=ksw)
    sol,solnet=Film.Isotherm(wPolyASD=1,T=303.15,psys=1E5,NET=True,RSvec=RHRS)
    resvec=[]
    wnet=np.asarray([solnet[i]["wi"] for i,val in enumerate(solnet)])
    for j,val in enumerate (yexp[0,:]):
        resvec.append(((yexp[:,j]-wnet[:,j])/yexp[:,j])**2) if j!=2 else None #API not considered here
    dis=sum(resvec)**0.5
    mindiidx=np.argmin(dis)
    RSvec.append(RHlin2[mindiidx])
    #yopt=yexp[:,i]*(1+opt[0])
    #yvec.append(yopt)
    rhooptvec.append(rhoopt)
rhooptvec=np.asarray(rhooptvec)
RSvec=np.asarray(RSvec)

pd.DataFrame(columns=["rho","RS"],data=np.asarray([rhooptvec,RSvec]).T).to_excel("results.xlsx")