import ProfilometerAuswertungFunktionen as PAF
import os
import numpy as np 
cwd = os.getcwd()
from os.path import join
import matplotlib.pyplot as plt


#M13NCD1.TXT PVP Film with strongly curved glass slide
#LEVLD1.TXT Indometacin Film with small to none curved glass slide
#window1(DataDir="Dektak")
#ASD Pure
#recent müsste Dateiname sein
from PySorption.dcc import window
window1=window(DataDir="Dektak")
path1=window1.result
Probe=window1.recent
molecule=window1.molecule+"_"+window1.molecule2 if hasattr(window1,"molecule2") else window1.molecule
wPolyASD=window1.wPolyASD
pathdata2=join(cwd,"Dektak","ASD",molecule,str(wPolyASD),Probe) if hasattr(window1,"molecule2") else join(cwd,"Dektak","Pure",molecule,Probe)

mass=float(input("Enter mass of film in mg\n"))
path1={0:path1} if isinstance(path1,str) else path1
#rhoreal=float(input("Enter real density of film in kg/m^3\n"))
#mass=2.24 #Indo Film
components=window1.components
rhoreal=(wPolyASD/components[0]["rho0"]+(1-wPolyASD)/components[2]["rho0"])**-1 if hasattr(window1,"molecule2") else components[1]["rho0"]
turns={0 : 0,
           1 : np.pi/2,
           2 : np.pi,
           3: np.pi/4,
           4: 3*np.pi/4}

Data=[]
Len=[]
path1
for kk in range(len(path1)):
    path2=join(pathdata2,path1[kk])
    Datared, Length= PAF.ReadData(path2)
    Data.append(Datared) # Measurement of height [nm]
    Len.append(Length) # Measurement Scan Length [micrometer]
  
name=Probe
rho=[]
D=[]
LD=[]
RD=[]
nR=[]
rvec=[]
shift1=300 #Do not consider the shift1 first and last Datapoints for jump detection
shift2=40 #Do not regress the shift2 first Datapoints for the correction polynomial
windowLength=500 #Only change if correction fails
plotrawdata=1 # 1 Shows curve correction of the glass slide 
for kk in range(len(path1)):
    rho1,D1,LD1,RD1,nR1,rvecR1=PAF.Dektak(Data[kk],Len[kk],mass,turns[kk],windowLength,plotrawdata,shift1,shift2,path1[kk][:-4],rhoreal)
    rho.append(rho1) #Returns the density in [kg/m^3]
    D.append(D1) #Returns the Diameter in [micrometer]
    LD.append(LD1) #Returns the Height in [nm]
    RD.append(RD1) #Returns the Height in [nm]
    nR.append(nR1) #Returns the number of Datapoints
    rvec.append(rvecR1) #Returns a equidistant vector for the radius in [micrometer]
Dav=np.average(D)
rhom=np.average(rho)
nAV=int(np.ceil(np.average(nR)))
rvecm=np.linspace(0,Dav/2,nAV)
Li,Ri=PAF.NormalizeOnSameRadius(Dav,D,rvecm,rvec,LD,RD)
string=["density[kg/m^3]","D[nm]"]
DatatoCSV=[rho,D]
PAF.WriteResultsCSV(string,DatatoCSV,name)
nPhi=100
import seaborn as sns
cmap=sns.dark_palette(color=[255/255,133/255,0],as_cmap="True") #orange
cmap=sns.dark_palette(color=[146/255,208/255,80/255],as_cmap="True") #green
cmap=sns.dark_palette(color=[236/255,175/255,0/255],as_cmap="True")#okka
#cmap='copper'
dens=100
PAF.ProfilometerPlot(Li,Ri,turns,nPhi,nAV,rvecm,cmap,dens,name)
