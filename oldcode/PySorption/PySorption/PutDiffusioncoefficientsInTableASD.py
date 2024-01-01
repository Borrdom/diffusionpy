import pandas as pd
import casadi as cs
from os import getcwd
from os.path import join
import xlsxwriter
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from PyCSAFT.get_pcpar import get_par
from matplotlib import cm
from scipy.interpolate import InterpolatedUnivariateSpline
from Pylomer.PyGT import Tggt
def Tg(wi,rhoi,Tgi):
    return np.sum(wi*1/rhoi,2)/np.sum(wi*1/rhoi/Tgi,2)
# Create an new Excel file and add a worksheet.
workbook = xlsxwriter.Workbook('TableDiffusioncoeff.xlsx')
worksheet = workbook.add_worksheet()

def roundFirst(x):
    res=[]
    for i,val in enumerate(x):
        val=np.fmax(val,1E-17)
        mag=10**np.floor(np.log10(np.abs(val)))
        res.append(np.round(np.round(val/mag,2)*mag,2+int(-(np.log10(mag)))))
    return np.asarray(res)

class Film():
    def __init__(self):
        self.molecule=1
        self.molecule2=1
self=Film()
self.molecule2="indomethacin"
self.molecule="pvpva64"
ASD=self.molecule+"_"+self.molecule2
dl1=0.8
dl2=0.5
ASD1=ASD+"_"+str(dl1)
ASD2=ASD+"_"+str(dl2)
#self.molecule="indomethacin"
#self.molecule2="ritonavir"
cwd=getcwd()
pathdiffcoeff1=join(cwd,"LeastFit","Fit",ASD1+".xlsx")
pathdiffcoeff2=join(cwd,"LeastFit","Fit",ASD2+".xlsx")




pure,kij=get_par([self.molecule,"water",self.molecule2],T=298.15)
Tgi=np.asarray([val["Tg"] for i,val in enumerate(pure)])
rhoi=np.asarray([val["rho0"] for i,val in enumerate(pure)])
if self.molecule2=="indomethacin":
    kfitindo=0.11
    rhofit=(rhoi[1]*Tgi[1])/(kfitindo*Tgi[2])
    rhoi[2]=rhofit
    print(rhoi[0])
if self.molecule2=="ritonavir":
    kfitindo=0.244
    rhofit=(rhoi[1]*Tgi[1])/(kfitindo*Tgi[2])
    rhoi[2]=rhofit
    print(rhoi[2])
w2=np.linspace(0,1,1000)
w101=np.ones_like(w2)*dl1
w102=np.ones_like(w2)*dl2
w301=1-w101
w302=1-w102
wnull=np.zeros_like(w2)

w11=(1-w2)*w101
w31=(1-w2)*w301
w12=(1-w2)*w102
w32=(1-w2)*w302
wwet1=np.asarray([w11,w2,w31])
wwet2=np.asarray([w12,w2,w32])
wdry1=np.asarray([w101,wnull,w301])
wdry2=np.asarray([w102,wnull,w302])
Tg1=Tggt(wwet1.T,rhoi,Tgi,q=0)
Tg2=Tggt(wwet2.T,rhoi,Tgi,q=0)
Tgd1=Tggt(wdry1.T,rhoi,Tgi,q=0)
Tgd2=Tggt(wdry2.T,rhoi,Tgi,q=0)
T=298.15#Tgmix[-1,-1]
Psi1=(Tgd1-Tg1)/(Tgd1-T)

Psi2=(Tgd2-Tg2)/(Tgd2-T)

Psi1_fun=InterpolatedUnivariateSpline(w2,Psi1)
Psi2_fun=InterpolatedUnivariateSpline(w2,Psi2)
#Tgmix[-1,-1] wasser
#Tgmix[0,0] molecule 1
#Tgmix[-1,0] water
#Tgmix[0,-1] molecule 2


def read_diffcoeff(filename):
    Lit=pd.read_excel(filename).dropna()
    DF=roundFirst(Lit["DCrancave/m/s^2 0"].values*10**15)
    DM=roundFirst(Lit["DStefanave/m/s^2 0"].values*10**15)
    DFs=roundFirst(Lit["DCrancstd/m/s^2 0"].values*10**15)
    DMs=roundFirst(Lit["DFStefanstd/m/s^2 0"].values*10**15)
    w=roundFirst(Lit["w2ave/- 0"].values*100)
    ws=roundFirst(Lit["w2std/- 0"].values*100)
    RH=roundFirst(Lit["RHave/- 0"].values*100)
    RHs=roundFirst(Lit["RHstd/- 0"].values*100)
    wm=roundFirst(Lit["w2barave/- 0"].values*100)
    #Dfun=cs.interpolant('LUT','linear',[w],DM)

    return DF,DM,w,RH,DFs,DMs,ws,RHs,wm

DF12,DM12,w12,RH12,DFs12,DMs12,ws12,RHs12,wm12= read_diffcoeff(pathdiffcoeff1)
DF23,DM23,w23,RH23,DFs23,DMs23,ws23,RHs23,wm23= read_diffcoeff(pathdiffcoeff2)
# Write some numbers, with row/column notation.
Psi12=roundFirst(Psi1_fun(wm12/100))*100
Psi23=roundFirst(Psi2_fun(wm23/100))*100
Table=[RH12,w12,w23,DF12,DF23,DM12,DM23,Psi12,Psi23]
Tablestd=[RHs12,ws12,ws23,DFs12,DFs23,DMs12,DMs23,np.zeros_like(Psi12),np.zeros_like(Psi23)]
for i, vali in enumerate(Table):
    for j, valj in enumerate(vali):
        worksheet.write(j, i, str(valj)+"\u00B1"+str(Tablestd[i][j]))
workbook.close()
