# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 07:09:38 2020

@author: domin
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from os import listdir
from os import getcwd
from os.path import join

# Values=listdir(path3)
# print(listdir(path3))
# molecule=input("Specify component\n")

# #Error if component not found
# if molecule not in Values:
#     print("Component not found abort\n")
#     import sys
#     sys.exit(0)

# files=listdir(path2)
# filesras=[s for s in files if "ras" in s]

from PySorption.dcc import window
window1=window(DataDir="PXRD")
molecule=window1.recent
fil=window1.result
filesras=[val for i,val in fil.items()] if isinstance(fil,dict) else [fil]
cwd=getcwd()
path3=join(cwd,"PXRD")
path2=join(path3,molecule)

def read_ras(path1,name):

    fobj1 = open(path1, errors='ignore').read()
    fobj2=fobj1.replace("\t","")
    fobj3=fobj2.replace("\n"," ")
    startidx1=fobj3.rfind("RAS_INT_START")+len("RAS_INT_START")
    endidx1=fobj3.rfind("*RAS_INT_END")
    Data1=fobj3[startidx1:endidx1]
    #Data1=Data1+("     ")
    Datared=[]
    j=0
    for i in range(len(Data1)):
        if i>j:
            if Data1[i]!=" ":
                k=0
                j=i
                if j<len(Data1):
                    while k==0:
                        j=j+1
                        if Data1[j]==" ":
                                k=1
                                Datared.append(float(Data1[i:j]))
    Datared=np.asarray(Datared)


    theta2=Datared[0::3]
    Signal=Datared[1::3]


    def PXRDPlot(theta2,Signal,name):
        dpi=96
        fs=(4604/(8*dpi),3602/(8*dpi))
        fig6, ax6 = plt.subplots(figsize=fs,dpi=dpi)
        ax6.plot(theta2,Signal,'k-',LineWidth=0.5)
        yax6=ax6.set_ylabel(r'Intensity'+'\n'+ r'      [-]')
        ax6.set_xlabel(r'2$\Theta$ [°]')
        ax6.set_xlim([np.min(theta2),np.max(theta2)])
        yax6.set_rotation(0)
        yax6.set_position([150/dpi,100/dpi,222/dpi,138.5/dpi])
        import PyPlorigin.PyOrigin as OriginVorlagePlot#added 03.03.2020 for origin
        xlab='Intensity%(CRLF)[-]'
        ylab='2\g(Q)[°]'
        x=[]
        y=[]
        xe=[]
        ye=[]
        xs=[]
        ys=[]
        x.append(np.asarray(theta2))
        y.append(np.asarray(Signal))
        xs.append(np.asarray(theta2))
        ys.append(np.asarray(Signal))
        xe.append(np.asarray([0]))
        ye.append(np.asarray([0]))

        xData,yData,xSim,ySim,xErr,yErr=OriginVorlagePlot.MusterDictOrigin()
        xu="°"
        yu="-"
        xc="Angle"
        yc="Signal"
        xData["values"],yData["values"],xSim["values"],ySim["values"],xErr["values"],yErr["values"]=x,y,xs,ys,xe,ye
        xData["unit"],yData["unit"],xSim["unit"],ySim["unit"],xErr["unit"],yErr["unit"]=xu,yu,xu,yu,xu,yu
        xData["comment"],yData["comment"],xSim["comment"],ySim["comment"],xErr["comment"],yErr["comment"]=xc,yc,xc,yc,xc,yc
        OriginVorlagePlot.Plot(xData,yData,xSim,ySim,xErr,yErr,ylab,xlab,os.path.basename(name))
        return
    PXRDPlot(theta2,Signal,name)

for i,val in enumerate(filesras):
    #path1=path2+"\\"+val
    path1=val#path2+"\\"+val
    read_ras(path1,val[:-4])
