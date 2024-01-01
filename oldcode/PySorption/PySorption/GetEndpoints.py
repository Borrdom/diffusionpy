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
    Te=Time[indexofchange[0]+1:indexofchange[1]].astype(float)-Time[indexofchange[0]+1]
    Me=Mass[indexofchange[0]+1:indexofchange[1]].astype(float)
    ne=Te.shape[0]
    from scipy.signal import argrelextrema
    #from scipy.interpolate import Rbf
    
    We=(Me-M0)/Me
    spacing=10
    
    from interparc import interparc
    T,W=interparc(Te,We,5000)
    
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
    #fig1,ax1=plt.subplots()
    #fig2,ax2=plt.subplots()
    
    #ax1.plot(T[::50][minima],dWdT_smooth[minima],"kx")
    #ax1.plot(Tred,ddWddT_smooth)
    #ax2.plot(Tred,dWdT_smooth)
    idxmin=np.argmin(ddWddT_smooth[:int(len(ddWddT_smooth)/2)])
    
    idxtend=(idxmin+next(i for i,v in enumerate(ddWddT_smooth[idxmin::]>=0) if v))*spacing+1
    #idxtend=np.argmin(ddWddT_smooth)*spacing+1
    #idxmine=int(idxmin/5000*ne)
    idxtende=int(idxtend/5000*ne)
    
    TSorption=Te[:idxtende]
    TKrist=Te[idxtende:]
    MSorption=We[:idxtende]
    MKrist=We[idxtende:]
    
    fig,ax=plt.subplots()
    ax.plot(TSorption,MSorption,'b-')
    ax.plot(TKrist,MKrist,'r--')
    
    #Crank Fit
    XSorption=MSorption/(1-MSorption)
    Xe=We/(1-We)
    Mzero,Minfty=XSorption[0],XSorption[-1]

    Mzeroe,Minftye=Xe[0],Xe[-1]

    
    def crank(t,d):
        ns=30
        rho=1000
        A=np.pi/4*0.0145**2
        l=M0*(1E-6)/rho/A
        mt_minf = 1-sum([(8/np.pi**2)*(1/(2*n + 1)**2)*np.exp((-d*((n + 1/2)**2)*(np.pi**2)*t)/l**2) for n in range(ns)])
        return mt_minf

    mperc=(XSorption-Mzero)/(Minfty-Mzero)
    mperce=(Xe-Mzeroe)/(Minftye-Mzeroe)

    from scipy.optimize import curve_fit

    popt, pcov = curve_fit(crank, TSorption*60, mperc,p0=1E-13)
    popte, pcove = curve_fit(crank, Te*60, mperce,p0=1E-13)

    XSorptionsim=crank(TSorption*60,popt)*(Minfty-Mzero)+Mzero
    MSorptionsim=XSorptionsim/(XSorptionsim+1)
    
    Xesim=crank(Te*60,popte)*(Minftye-Mzeroe)+Mzeroe
    Wesim=Xesim/(Xesim+1)
    ###
    
    
    
    import os
    filename=os.path.basename(original)
    dirname=os.path.dirname(original)
    ax.set_title(filename)
    Extracted=pd.DataFrame([TSorption,MSorption,TKrist,MKrist,Te,We,np.asarray([M0]),TSorption,MSorptionsim,popt,Te,Wesim,popte]).T
    Extracted.columns=["t[min]","w[-]","tK[min]","wK[-]","tges[min]","wges[-]","M0[mg]","tsim[min]","wsim[-]","D[m^2/s]","tgessim[min]","wgessim[-]","Dges[m^2/s]"]
    dateiname, file_extension = os.path.splitext(filename)
    
    #dateiname=original.split(".")[-1]
    #if not os.path.exists(os.path.join(dirname,"Extracted_files")):
    #    os.makedirs(os.path.join(dirname,"Extracted_files"))
    newfilename=os.path.join(dirname,dateiname+"_Extracted"+".xlsx") 
    Extracted.to_excel(newfilename,index=None)
    return Extracted
    

    
    
if __name__=="__main__":
    original=askopenfilenames()
    Extractedlist=[Extract(val) for i,val in enumerate(original)]
    plt.show()