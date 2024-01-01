import os
from os.path import join
import numpy as np
#import PureDataBase
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
def ReadDiffcoeffFromCSV(molecule,molecule2,wPolyASD):
    cwd=os.getcwd()

    pathdiffcoeffp=join(cwd,"LeastFitNET",molecule+".csv") 

    pathdiffcoeffa=join(cwd,"LeastFit",molecule2+".csv") 

    pathdiffcoeffasd=join(cwd,"LeastFitNET",molecule+molecule2+str(wPolyASD)+".csv")
    
    import pandas as pd
    def DataBaseImport(*args):
        """Enter components as tuple, get tuple of dicts"""
        return tuple(getattr(PureDataBase,val) for i,val in enumerate(args))
    pol,water,api=DataBaseImport(molecule,"water",molecule2)
    
    def read_diffcoeff(filename):
        try:
            Lit=pd.read_csv(filename).dropna() 
        except:
            print("No Data for "+pol["name"]+" and "+api["name"]+" with DL of "+str(wPolyASD*100)+"%")
            Lit=pd.read_csv(join(cwd,"LeastFit","pvp"+".csv")).dropna() 
        Ð=Lit["DStefanH20bar[-]"].values
        D=Lit["DFickH20bar[-]"].values
        ww=Lit["wH20mittelbar[-]"].values
        return Ð,D,ww
    Ð,D,ww= read_diffcoeff(pathdiffcoeffasd)
    Ðp,Dp,wwp= read_diffcoeff(pathdiffcoeffp)
    Ða,Da,wwa= read_diffcoeff(pathdiffcoeffa)
    
    Tg0=np.asarray([val["Tg"] for i,val in enumerate([pol,water,api])])
    rho0=np.asarray([val["rho0"] for i,val in enumerate([pol,water,api])])
    from TgGordnonTaylor import TG2
    T=298.15
    Tga=T-TG2(wwa,1,Tg0,rho0)
    Tgp=T-TG2(wwp,0,Tg0,rho0)
    Tga=TG2(wwa,1,Tg0,rho0)
    Tgp=TG2(wwp,0,Tg0,rho0)
    
    
    Tgamax=T-TG2(0,1,Tg0,rho0)
    Tgpmax=T-TG2(0,0,Tg0,rho0)
    Tgamax=TG2(0,1,Tg0,rho0)
    Tgpmax=TG2(0,0,Tg0,rho0)
    #Tga=(Tga-np.min(Tga))/(abs(np.min(Tga))-np.min(Tga))
    #Tgp=(Tgp-np.min(Tgp))/(abs(np.min(Tgp))-np.min(Tgp))  
    #Tga=(Tga-np.min(Tga))/(2*abs(np.min(Tga)))
    #Tgp=(Tgp-np.min(Tgp))/(2*abs(np.min(Tgp))) 
    
    #Tga=(Tga-Tgamax)/(2*Tgamax)
    #Tgp=(Tgp-Tgpmax)/(2*Tgpmax)
    Tga=(Tgamax-Tga)/(abs(T-Tgamax))
    Tgp=(Tgpmax-Tgp)/(abs(T-Tgpmax))

    Ðafun=InterpolatedUnivariateSpline(Tga,Ða,ext=3)
    Ðpfun=InterpolatedUnivariateSpline(Tgp,Ðp,ext=3)
    #Ðafun=InterpolatedUnivariateSpline(wwa,Ða,ext=0)
    #Ðpfun=InterpolatedUnivariateSpline(wwp,Ðp,ext=0)
    
    TgASD=T-TG2(ww,1-wPolyASD,Tg0,rho0)
    TgASDmax=T-TG2(0,1-wPolyASD,Tg0,rho0)
    
        
    TgASD=TG2(ww,1-wPolyASD,Tg0,rho0)
    TgASDmax=TG2(0,1-wPolyASD,Tg0,rho0)
    #TgASD=(TgASD-np.min(TgASD))/(2*abs(np.min(TgASD))
    #TgASD=(TgASD-np.min(TgASD))/(abs(np.min(TgASD))-np.min(TgASD))
    TgASD=(TgASDmax-TgASD)/(abs(T-TgASDmax))
    #TgASDfun=InterpolatedUnivariateSpline(ww,TgASD,ext=3)



    def convert_fraction(*args):
        return tuple(j/sum(args) for i,j in enumerate(args))
    x1FEED,x2FEED,x3FEED=convert_fraction(wPolyASD/pol["Mi"],0,(1-wPolyASD)/api["Mi"]) # 19.08.20 selbe Ergebnisse wie bei Leusmann für PVPIND 80 TOP
    
    Ddict={"Ð":Ð,"D":D,"ww":ww,"Ðp":Ðpfun,"Dp":Dp,"wwp":wwp,"Ða":Ðafun,"Da":Da,"wwa":wwa,"x3FEED":x3FEED,"TgASD":TgASD}
    return Ddict

polylist=["pvpva64","pvp"]
apilist=["indometacin","ritonavir"]
DLlist=[0.8,0.5]

def PredictASDCoeffPlot(pol,api,wPolyASD):

    VR8=ReadDiffcoeffFromCSV(pol,api,wPolyASD)
    ÐVR8,DVR8,XVR8,wVR8=VR8["Ð"],VR8["D"],VR8["x3FEED"],VR8["ww"]
    TgASD=VR8["TgASD"]
    def DWASDFEED(Dp,Da,X3):
        return ((1-X3)/Dp+X3/Da)**-1
    ÐV,DV,wV=VR8["Ðp"],VR8["Dp"],VR8["wwp"]
    ÐR,DR,wR=VR8["Ða"],VR8["Da"],VR8["wwa"]
    ÐSIM=DWASDFEED(ÐV(TgASD),ÐR(TgASD),XVR8)
    #ÐSIM=DWASDFEED(ÐV(wVR8),ÐR(wVR8),XVR8)
    
    
    fig,ax1=plt.subplots()
    ax1.plot(wVR8,ÐVR8,"ko")
    ax1.plot(wVR8,ÐSIM,"kx")
    ax1.set_title(pol+"_"+api+"_"+str(wPolyASD))
    ax1.set_xlabel("w[-]")
    ax1.set_ylabel(r'Ð[$m^2/s$]')
    fig.savefig(pol+"_"+api+"_"+str(wPolyASD)+"NET"+".jpeg")


[[[PredictASDCoeffPlot(val1,val2,val3) for i,val1 in enumerate(polylist)] for j,val2 in enumerate(apilist)] for k,val3 in enumerate(DLlist)]

    
# VR8=ReadDiffcoeffFromCSV("pvpva64","ritonavir",0.8)
# ÐVR8,DVR8,XVR8,wVR8=VR8["Ð"],VR8["D"],VR8["x3FEED"],VR8["ww"]

# VR5=ReadDiffcoeffFromCSV("pvpva64","ritonavir",0.5)
# ÐVR5,DVR5,XVR5,wVR5=VR5["Ð"],VR5["D"],VR5["x3FEED"],VR5["ww"]

# PI8=ReadDiffcoeffFromCSV("pvp","indometacin",0.8)
# ÐPI8,DPI8,XPI8,wPI8=PI8["Ð"],PI8["D"],PI8["x3FEED"],PI8["ww"]

# PI5=ReadDiffcoeffFromCSV("pvp","indometacin",0.5)
# ÐPI5,DPI5,XPI5,wPI5=PI5["Ð"],PI5["D"],PI5["x3FEED"],PI5["ww"]

# VI8=ReadDiffcoeffFromCSV("pvpva64","indometacin",0.8)
# ÐVI8,DVI8,XVI8,wVI8=VI8["Ð"],VI8["D"],VI8["x3FEED"],VI8["ww"]

# VI5=ReadDiffcoeffFromCSV("pvpva64","indometacin",0.5)
# ÐVI5,DVI5,XVI5,wVI5=VI5["Ð"],VI5["D"],VI5["x3FEED"],VI5["ww"]

# ÐP,DP,wP=PI5["Ðp"],PI5["Dp"],PI5["wwp"]
# ÐV,DV,wV=VI5["Ðp"],VI5["Dp"],VI5["wwp"]
# ÐI,DI,wI=VI5["Ða"],VI5["Da"],VI5["wwa"]
# ÐR,DR,wR=VR5["Ða"],VR5["Da"],VR5["wwa"]

# def DWASDFEED(Dp,Da,X3):
#     return ((1-X3)/Dp+X3/Da)**-1
# ÐVR8S=DWASDFEED(ÐV,ÐR,XVR8)
# ÐVR5S=DWASDFEED(ÐV,ÐR,XVR5)
# ÐVI8S=DWASDFEED(ÐV,ÐI,XVI8)
# ÐVI5S=DWASDFEED(ÐV,ÐI,XVI5)
# ÐPI8S=DWASDFEED(ÐP,ÐI,XPI8)
# ÐPI5S=DWASDFEED(ÐP,ÐI,XPI5)


# fig,ax=plt.subplots(3,1)
# ax[0].plot(wVR8,DVR8,"ko")
# ax[0].plot(wVR5,DVR5,"go")
# ax[0].set_title("pvpva64"+"_"+"ritonavir")
# ax[0].set_xlabel("w[-]")
# ax[0].set_ylabel(r'D[$m^2/s$]')
# ax[1].plot(wVI8,DVI8,"ko")
# ax[1].plot(wVI5,DVI5,"go")
# ax[1].set_title("pvpva64"+"_"+"indometacin")
# ax[1].set_xlabel("w[-]")
# ax[1].set_ylabel(r'D[$m^2/s$]')
# ax[2].plot(wPI8,DPI8,"ko")
# ax[2].plot(wPI5,DPI5,"go")
# ax[2].set_title("pvp"+"_"+"indometacin")
# ax[2].set_xlabel("w[-]")
# ax[2].set_ylabel(r'D[$m^2/s$]')


# fig,ax1=plt.subplots(3,1)
# ax1[0].plot(wVR8,ÐVR8,"ko")
# #ax1[0].plot(wVR8,ÐVR8S,"k-")
# ax1[0].plot(wVR5,ÐVR5,"go")
# #ax1[0].plot(wVR5,ÐVR5S,"g-")
# ax1[0].set_title("pvpva64"+"_"+"ritonavir")
# ax1[0].set_xlabel("w[-]")
# ax1[0].set_ylabel(r'Ð[$m^2/s$]')
# ax1[1].plot(wVI8,ÐVI8,"ko")
# ax1[1].plot(wVI5,ÐVI5,"go")
# ax1[1].set_title("pvpva64"+"_"+"indometacin")
# ax1[1].set_xlabel("w[-]")
# ax1[1].set_ylabel(r'Ð[$m^2/s$]')
# ax1[2].plot(wPI8,ÐPI8,"ko")
# ax1[2].plot(wPI5,ÐPI5,"go")
# ax1[2].set_title("pvp"+"_"+"indometacin")
# ax1[2].set_xlabel("w[-]")
# ax1[2].set_ylabel(r'Ð[$m^2/s$]')
