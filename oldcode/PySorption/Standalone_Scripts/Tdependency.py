import pandas as pd
import os
from os.path import join
import numpy as np
import lmfit
from lmfit import minimize,Parameters
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
cwd=os.getcwd()
filename=join(cwd,"EGModul_PVPVA.xlsx")
filenameT="DMA_PVPVA_Tscan3.xlsx"
def readcolumforstr(book,valj):
    lisidx=[np.flatnonzero(book[val]==valj) for i,val in enumerate(book) if (book[val]==valj).any()]
    lisbool=[(book[val]==valj).any() for i,val in enumerate(book)]
    lishead=[book.columns[i]  for i,val in enumerate(book) if (book[val]==valj).any()]
    return book[lishead[0]].values[3::].astype(float)
def readEGModulxlsx(filename,sheet):
    book=pd.read_excel(filename, encoding_override = "utf-8",sheet_name=sheet)
    tagfreq="Frequenz"
    tagstore="E'"
    tagloss="E''"
    Data=[]
    for j,valj in enumerate([tagfreq,tagstore,tagloss]):
        Data.append(readcolumforstr(book,valj))
    Data[0]=np.log10(Data[0])
    Data=tuple(Data)
    return  Data
def readEGModulxlsxT(filename):
    book=pd.read_excel(filename, encoding_override = "utf-8")
    tagTl="TL[°C]"
    tagTs="TS[°C]"
    tagstore="S[MPa]"
    tagloss="L[MPa]"
    
    Data=[]
    for j,valj in enumerate([tagTl,tagTs,tagstore,tagloss]):
        Data.append(readcolumforstr(book,valj))
    Data=tuple(Data)
    return  Data
sheet="120°C"
logomega,storage,loss=readEGModulxlsx(filename,sheet)
TL,TS,S,L=readEGModulxlsxT(filenameT)
S,L=S*1E6,L*1E6 #MPa in Pascal
omegaT0=0
storage_fun=InterpolatedUnivariateSpline(storage,logomega)
loss_fun=InterpolatedUnivariateSpline(storage,logomega)

tTS=storage_fun(S)
tTL=loss_fun(L)
lb=logomega[0]
ub=logomega[-1]


# remove storage bound
ETS=S[tTS>lb]
TTS=TS[tTS>lb]
tTS=tTS[tTS>lb]

ETS=ETS[tTS<ub]
TTS=TTS[tTS<ub]
tTS=tTS[tTS<ub]

#lremove loss bound
ETL=L[tTL>lb]
TTL=TL[tTL>lb]
tTL=tTL[tTL>lb]

ETL=ETL[tTL<ub]
TTL=TTL[tTL<ub]
tTL=tTL[tTL<ub]



ShiftingS,ShiftingL=tTS-omegaT0,tTL-omegaT0
T0=105
delTTS=TTS-T0
delTTL=TTL-T0

from scipy.optimize import curve_fit
C1=-8.86
C2=101.6
def WLF(delT,C1,C2):
    return C1*delT/(C2+delT)
Shifting_sim=WLF(delTTS,C1,C2)

popt, pcov =curve_fit(WLF,delTTS,ShiftingS,p0=[C1,C2])
C1opt,C2opt=tuple(popt)
Shifting_opt=WLF(delTTS,C1opt,C2opt)
plt.plot(delTTS,Shifting_sim)
plt.plot(delTTS,Shifting_opt)
plt.plot(delTTS,ShiftingS,'rx')
plt.plot(delTTL,ShiftingL,'bx')