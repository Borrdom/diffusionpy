import pandas as pd
import os
from os.path import join
import numpy as np
import lmfit
from lmfit import minimize,Parameters
import matplotlib.pyplot as plt
cwd=os.getcwd()
filename=join(cwd,"EGModul_PVPVA_Friederike.xlsx")
filename=join(cwd,"EGModul_PVPVAShifting2.xlsx")


def EModulPlot(logomega,logomegascan,storage,storagesim,loss,losssim,name):
    import OriginVorlagePlot#added 03.03.2020 for origin
    xlab='E*%(CRLF)/GPa'
    ylab='\g(t)/s'
    x=[]
    y=[]
    xe=[]
    ye=[]
    xs=[]
    ys=[]
    x.append(np.asarray(10**(-logomega)))
    x.append(np.asarray(10**(-logomega)))
    y.append(np.asarray(storage*1E-9))
    y.append(np.asarray(loss*1E-9))
    xs.append(np.asarray(10**(-logomegascan)))
    xs.append(np.asarray(10**(-logomegascan)))
    ys.append(np.asarray(storagesim*1E-9))
    ys.append(np.asarray(losssim*1E-9))
    xe.append(np.asarray(0))
    ye.append(np.asarray(0))
    xe.append(np.asarray(0))
    ye.append(np.asarray(0))
    
    xData,yData,xSim,ySim,xErr,yErr=OriginVorlagePlot.MusterDictOrigin()
    xu="s"
    yu="Pa"
    xc="relaxation time"
    yc="EModulus"
    xData["values"],yData["values"],xSim["values"],ySim["values"],xErr["values"],yErr["values"]=x,y,xs,ys,xe,ye
    xData["unit"],yData["unit"],xSim["unit"],ySim["unit"],xErr["unit"],yErr["unit"]=xu,yu,xu,yu,xu,yu
    xData["comment"],yData["comment"],xSim["comment"],ySim["comment"],xErr["comment"],yErr["comment"]=xc,yc,xc,yc,xc,yc
    OriginVorlagePlot.Plot(xData,yData,xSim,ySim,xErr,yErr,ylab,xlab,name)
    return
def ParameterPlot(tauJ,EJ,name):
    import OriginVorlagePlot#added 03.03.2020 for origin
    xlab='E%(CRLF)[GPa]'
    ylab='g(t)/s'
    x=[]
    y=[]
    xe=[]
    ye=[]
    xs=[]
    ys=[]
    x.append(np.asarray(10**(tauJ)))

    y.append(np.asarray(EJ*1E-9))
    etaJ=tauJ*EJ
    xs.append(np.asarray(etaJ*1E-9))

    ys.append(np.asarray(0))

    xe.append(np.asarray(0))


    ye.append(np.asarray(0))
    
    xData,yData,xSim,ySim,xErr,yErr=OriginVorlagePlot.MusterDictOrigin()
    xu="s"
    yu="Pa"
    xc="relaxation time"
    yc="modulus"
    xData["values"],yData["values"],xSim["values"],ySim["values"],xErr["values"],yErr["values"]=x,y,xs,ys,xe,ye
    xData["unit"],yData["unit"],xSim["unit"],ySim["unit"],xErr["unit"],yErr["unit"]=xu,yu,xu,yu,xu,yu
    xData["comment"],yData["comment"],xSim["comment"],ySim["comment"],xErr["comment"],yErr["comment"]=xc,yc,xc,yc,xc,yc
    OriginVorlagePlot.Plot(xData,yData,xSim,ySim,xErr,yErr,ylab,xlab,name)
    return

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
def GenMaxwell(par,logomega,storage,loss):
    omega=10**(logomega)
    tauJ=10**np.asarray([par["tauJ"+str(i)].value  for i in range(int(par["nArm"].value))])
    EJ=np.asarray([par["EJ"+str(i)].value  for i in range(int(par["nArm"].value))])
    etaJ=EJ*tauJ
    Eeq=par["Eeq"]
    etaEq=par["etaeq"]
    E=sum([Eeq+(omega*1j)/(1+val*omega*1j)*etaJ[k]+etaEq*omega*1j for k,val in enumerate(tauJ)])
    return E.real-storage,E.imag-loss
sheetvec=["70°C75RH","60°C75RH","50°C75RH","40°C75RH","80°C60RH","70°C60RH"]
sheetvec=["100°C","110°C","120°C"]
Data=[readEGModulxlsx(filename,val) for i,val in enumerate(sheetvec)]
#idxfreq=np.flatnonzero(book==tagfreq)
#idxstore=np.flatnonzero(book==tagstore)
#idxloss=np.flatnonzero(book==tagloss)
#book[idxloss]

par=Parameters()
nArm=4
tauJ0=np.linspace(-15,15,nArm+1)
[par.add('tauJ'+str(i), value=val, min=tauJ0[i-1], max=val,vary=1) for i,val in enumerate(tauJ0)]
[par.add('EJ'+str(i), value=0, min=0, max=1E20,vary=1) for i in range(nArm)]
par.add('Eeq',value=0,min=0, max=1E9,vary=0)
par.add('etaeq',value=0,min=0, max=1E11,vary=0)
par.add('nArm',value=nArm,vary=0 ,min=0, max=30)
par.add('etacross',value=1E8,vary=10 ,min=0, max=1E20)
par.add('taucross',value=1,vary=0 ,min=-10, max=20)
par.add('expcross',value=1,vary=0 ,min=0, max=2)


def FitData(par,logomega,storage,loss):
    def residual(par,logomega,storage,loss):
        res1,res2=GenMaxwell(par,logomega,storage,loss)
        res1=res1+storage
        res2=res2+loss
        vis=(res2**2+res1**2)**(1/2)
        visexp=(loss**2+storage**2)**(1/2)
        return vis-visexp
    # def residual(par,logomega,storage,loss):
    #     etacross=par["etacross"]
    #     taucross=par["taucross"]
    #     expcross=par["expcross"]
    #     eta=etacross/(1+(etacross*10**(logomega)/10**taucross)**(expcross-1))
    #     visexp=(loss**2+storage**2)**(1/2)
    #     res=(eta-visexp)/1E8
    #     print(res)
    #     return res
    
    o1 = minimize(residual, par, args=(logomega,storage,loss), method='powell')
    par=o1.params
    storagesim,losssim=GenMaxwell(par,logomega,storage,loss)
    return par,storagesim+storage,losssim+loss



paroptvec=[]
fig,ax=plt.subplots(2,len(Data),figsize=(9,6))
EJvec=[]
tauJvec=[]
etaeqvec=[]
viscosity=[]
logomegascan=np.linspace(-5,5,100)
for i,val in enumerate(Data):
    logomega=val[0]
    storage=val[1]
    loss=val[2]
    paropt,storagesim,losssim=FitData(par,logomega,storage,loss)
    #paropt["EJ0"].value=456704000
    #paropt["EJ1"].value=321794000
    #paropt["EJ2"].value=122470000
    #paropt["EJ3"].value=26876500
    #paropt["tauJ0"].value=np.log10(0.0774479630649733)
    #paropt["tauJ1"].value=np.log10(0.36881875713636)
    #paropt["tauJ2"].value=np.log10(1.37361808470257)
    #paropt["tauJ3"].value=np.log10(10.5053372740917)
    storagesim2,losssim2=GenMaxwell(paropt,logomegascan,np.zeros_like(logomegascan),np.zeros_like(logomegascan))
    #Try with values from Sorption
    # par["EJ0"].value=.7E7
    # par["EJ1"].value=300E7
    # par["EJ2"].value=13E7
    # par["EJ3"].value=.1E7
    # par["tauJ0"].value=np.log10(1.6E-3)# Works pretty good with sorption step
    # par["tauJ1"].value=np.log10(1.77E-1)
    # par["tauJ2"].value=np.log10(.376)
    # par["tauJ2"].value=np.log10(4)
    # storagesim,losssim=GenMaxwell(par,logomega,storage,loss)
    # storagesim,losssim=storagesim+storage,losssim+loss
    # paropt["EJ0"].value=240E7
    # paropt["EJ1"].value=.37E7
    # paropt["EJ2"].value=5E7
    # paropt["EJ3"].value=.01E7
    # paropt["tauJ0"].value=np.log10(1.77E-1)# Works pretty good with sorption step
    # paropt["tauJ1"].value=np.log10(1.59E-1)
    # paropt["tauJ2"].value=np.log10(2.5)
    # paropt["tauJ2"].value=np.log10(43)
    # storagesim,losssim=GenMaxwell(paropt,logomega,storage,loss)
    # storagesim,losssim=storagesim+storage,losssim+loss
    # par["EJ0"].value=299E7
    # par["EJ1"].value=2.86E7
    # par["tauJ0"].value=np.log10(1.19E-1)# Works pretty good with sorption step
    # par["tauJ1"].value=np.log10(5.59)
    # storagesim,losssim=GenMaxwell(par,logomega,storage,loss)
    # storagesim,losssim=storagesim+storage,losssim+loss
    # par["EJ0"].value=387E7
    # par["EJ1"].value=255E7
    # par["EJ1"].value=2.33E7
    # par["tauJ0"].value=np.log10(1.13E-1)# Works pretty good with sorption step
    # par["tauJ1"].value=np.log10(1.08E-1)
    # par["tauJ1"].value=np.log10(8.206)
    # storagesim,losssim=GenMaxwell(par,logomega,storage,loss)
    # storagesim,losssim=storagesim+storage,losssim+loss
    #Enter Your Model
    om=np.linspace(-15,15,300)
    #para=Parameters()
    #nArm=1
    #para.add('tauJ0', value=np.log10(6))
    #para.add('EJ0', value=2E10)
    #para.add('nArm',value=nArm)
    #para.add('Eeq',value=0,min=0, max=1E9,vary=0)
    #para.add('etaeq',value=0,min=0, max=1E9,vary=0)
    #a,b=GenMaxwell(para,om,np.zeros_like(om),np.zeros_like(om))
    aopt,bopt=GenMaxwell(paropt,om,np.zeros_like(om),np.zeros_like(om))
    #+2.2 emprical shifting
    #ax[0,i].plot(om+2.2,a,'r-.')
    #ax[0,i].plot(om+2.2,b,'b-.')
    
    ax[0,i].plot(logomega,np.log10(storage),'rx')
    ax[0,i].plot(logomega,np.log10(loss),'bx')
    ax[0,i].plot(om,np.log10(aopt),'r-')
    ax[0,i].plot(om,np.log10(bopt),'b-')
    #ax[0,i].set_xlim([logomega[0],logomega[-1]])
    ax[0,i].set_xlabel(r'$\log{\omega}$[1/s]')
    ax[0,i].set_ylabel(r'E[Pa]')
    #ax[0,i].set_ylim([0,3E9])
    #ax[0,i].set_xlim([-15,15])
    paroptvec.append(paropt)
    
    tauJ=np.asarray([paropt["tauJ"+str(j)].value  for j in range(int(paropt["nArm"].value))])
    EJ=np.asarray([paropt["EJ"+str(j)].value  for j in range(int(paropt["nArm"].value))])
    tauJvec.append(tauJ)
    EJvec.append(EJ)
    etaeqvec.append(paropt["etaeq"].value) 
    ax[1,i].plot(-tauJ,np.log10(EJ),'bx')
    #ax[1,i].set_xlim([logomega[0],logomega[-1]]) # -0.5645828264742274 shifting
    #ax[1,i].set_xlim([-5,5])
    #ax[1,i].set_ylim([0,1E9])
    ax[1,i].set_xlabel(r'$\log{1/\tau}$[1/s]')
    ax[1,i].set_ylabel(r'E[Pa]')
    
    #-0.5645828264742274
    fig4,ax4=plt.subplots()
    ax4.set_xlabel(r'$\log{\tau}$[s]')
    ax4.set_ylabel(r'E[Pa]')
    ax4.plot(tauJ,EJ,'ko')
    etaeq=paropt
#fig12,ax12=plt.subplots()    
#ax12.plot(np.log10(etaeqvec))
#EModulPlot(logomega,logomegascan,storage,storagesim2,loss,losssim2,"EGModul_PVPVA.opj")
#ParameterPlot(tauJ,EJ,"Parameter_PVPVA.opj")
#a,b=GenMaxwell(par,np.linspace(-10,10,300),np.zeros((300)),np.zeros((300)))
#viscosity=np.cumsum(a[:-1]+a[1:])/2*np.diff(-10**(-np.linspace(-10,10,300)))
#fig11,ax11=plt.subplots()
#ax11.plot(10**(-np.linspace(-10,10,300)[1:]),viscosity,'rx')
#ax.plot(omega,loss,'bx')

EJstr=str(EJvec[2].tolist())
taustr=str((10**tauJvec[2]).tolist())
np.savetxt("RheoPar.txt",np.asarray([EJstr,taustr]),fmt="%s")
    
    
