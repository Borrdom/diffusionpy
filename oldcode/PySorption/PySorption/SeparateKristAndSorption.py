### read file
from scipy.ndimage import median_filter
from scipy.interpolate import InterpolatedUnivariateSpline
import pandas as pd
from os.path import join
from os import getcwd
import tkinter
from tkinter import *
from tkinter.filedialog import askopenfilenames
import matplotlib.pyplot as plt
import numpy as np

def Extract(original):
    Data=pd.read_excel(original)
    Datanumberheader=Data[Data['DVS-Advantage-Plus-Data-File'].str.contains("minutes").fillna(False)]
    idx=Datanumberheader.index[0]
    Dataupper=Data.loc[:idx-1]
    Datared=Data.loc[idx:]
    Datared2=Data.loc[idx+1:]
    Datared2.columns=Datared.iloc[0]
    
    Time=Datared2["Time [minutes]"].values
    Mass=Datared2["m_Korr[mg]"].values
    Humidity=Datared2["Target Partial Pressure (Solvent A) [%]"].values
    indexofchange=np.where(np.diff(Humidity))[0]
    TDry=Time[:indexofchange[0]]
    MDry=Mass[:indexofchange[0]]
    MDry2=Mass[indexofchange[-1]:]
    TDry2=Time[indexofchange[-1]:]
    M0=np.fmin(np.min(MDry),np.min(MDry2))
    
    indexofchange=indexofchange if len(indexofchange)>1 else np.hstack((indexofchange,-1))
    T=Time[indexofchange[0]+1:indexofchange[1]].astype(float)
    M=Mass[indexofchange[0]+1:indexofchange[1]].astype(float)
    
    from scipy.signal import argrelextrema
    #from scipy.interpolate import Rbf
    
    W=(M-M0)/M
    spacing=10
    
    from interparc import interparc
    T,W=interparc(T,W,5000)
    
    Tred,Wred=T[::spacing],W[::spacing]
    
    dWdT=np.gradient(Wred)
    dWdT_smooth=median_filter(dWdT,50,mode="nearest")
    ddWddT_smooth=median_filter(np.gradient(dWdT_smooth),50,mode="nearest")
    
    #dWdT_smooth=Rbf(T[::50],dWdT, function="quintic", smooth=1)(T[::50])
    #ddWddT_smooth=np.gradient(dWdT_smooth)
    #minima=argrelextrema(dWdT_smooth,np.greater)[0]
    #dWdT_smooth=median_filter(np.gradient(dWdT_smooth),100,mode="nearest")
    #dWdT_smooth[dWdT_smooth==0]=10
    #minimum=np.min(dWdT_smooth)
    #dWdT_smooth[dWdT_smooth==10]=minimum*0.1
    fig1,ax1=plt.subplots()
    fig2,ax2=plt.subplots()
    
    #ax1.plot(T[::50][minima],dWdT_smooth[minima],"kx")
    ax1.plot(Tred,ddWddT_smooth)
    ax2.plot(Tred,dWdT_smooth)
    idxmin=np.argmin(ddWddT_smooth[:int(len(ddWddT_smooth)/2)])
    
    idxtend=(idxmin+next(i for i,v in enumerate(ddWddT_smooth[idxmin::]>=0) if v))*spacing+1
    #idxtend=np.argmin(ddWddT_smooth)*spacing+1
    
    TSorption=T[:idxtend]
    TKrist=T[idxtend:]
    MSorption=W[:idxtend]
    MKrist=W[idxtend:]
    
    fig,ax=plt.subplots()
    ax.plot(TSorption,MSorption,'b-')
    ax.plot(TKrist,MKrist,'r--')
    
    
    Extracted=pd.DataFrame([TSorption,MSorption,TKrist,MKrist,T,W,np.asarray([M0])]).T
    import os
    Extracted.columns=["t[min]","w[-]","tK[min]","wK[-]","tges[min]","wges[-]","M0[mg]"]
    dateiname, file_extension = os.path.splitext(original)
    #dateiname=original.split(".")[-1]
    newfilename=dateiname+"Extracted"+".xlsx"
    Extracted.to_excel(newfilename,index=None)
    
if __name__=="__main__":
    original=askopenfilenames()
    [Extract(val) for i,val in enumerate(original)]