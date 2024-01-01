import ProfilometerAuswertungFunktionen as PAF
import os
import numpy as np 
cwd = os.getcwd()
from os.path import join

a =join(cwd,'Dektak')
print(os.listdir(a))


path1={0 : join(cwd,'Dektak','HB1P.TXT'), #HPMCAS MK2
           1 : join(cwd,'Dektak','HB1M.TXT')}
#path1={0 : join(cwd,'Dektak','HB2P.TXT'), #HPMCAS MK3 
#           1 : join(cwd,'Dektak','HB2M.TXT')}
#path1={0 : join(cwd,'Dektak','HB4M.TXT'), #Ritonavir Film unzureichend!
#           1 : join(cwd,'Dektak','HB4M.TXT')} #
#path1={0 : join(cwd,'Dektak','HB5M.TXT'), #Ritonavir Film ok aber zu eienr Seite !
#           1 : join(cwd,'Dektak','HB5P.TXT')} #
#path1={0 : join(cwd,'Dektak','HB6P.TXT'), #ASD FILM gut #M is corrupted
#           1 : join(cwd,'Dektak','HB6P.TXT')} #
path1={0 : join(cwd,'Dektak','HB7P.TXT'), #ASD FILM sehr gut
           1 : join(cwd,'Dektak','HB7M.TXT')} #
#path1={0 : join(cwd,'Dektak','HB8P.TXT'), #ASD FILM sehr gut
#           1 : join(cwd,'Dektak','HB8M.TXT')} #
#path1={0 : join(cwd,'Dektak','CALIB4.TXT'), #HPMCAS MK2
#           1 : join(cwd,'Dektak','CALIB4.TXT')}
pathA={0 : join(cwd,'Dektak','PVPA1D1.TXT'),
           1 : join(cwd,'Dektak','PVPA1D2.TXT'),
            2 : join(cwd,'Dektak','PVPA1D3.TXT'),
            3 : join(cwd,'Dektak','PVPA1D4.TXT')}
pathB={0 : join(cwd,'Dektak','PVPB1D1.TXT'),
           1 : join(cwd,'Dektak','PVPB1D2.TXT'),
            2 : join(cwd,'Dektak','PVPB1D3.TXT'),
            3 : join(cwd,'Dektak','PVPB1D4.TXT')}
pathC={0 : join(cwd,'Dektak','PVPC1D1.TXT'),
           1 : join(cwd,'Dektak','PVPC1D2.TXT'),
            2 : join(cwd,'Dektak','PVPC1D3.TXT'),
            3 : join(cwd,'Dektak','PVPC1D4.TXT')}
pathD={0 : join(cwd,'Dektak','PVPD1D1.TXT'),
           1 : join(cwd,'Dektak','PVPD1D2.TXT'),
            2 : join(cwd,'Dektak','PVPD1D3.TXT'),
            3 : join(cwd,'Dektak','PVPD1D4.TXT')}
pathE={0 : join(cwd,'Dektak','PVPE1D1.TXT'),
           1 : join(cwd,'Dektak','PVPE1D2.TXT'),
            2 : join(cwd,'Dektak','PVPE1D3.TXT'),
            3 : join(cwd,'Dektak','PVPE1D4.TXT')}
path14={0 : join(cwd,'Dektak','PVP13D1.TXT'),
           1 : join(cwd,'Dektak','PVP13D2.TXT'),
            2 : join(cwd,'Dektak','PVP13D3.TXT'),
            3 : join(cwd,'Dektak','PVP13D4.TXT')}
pathIND={0 : join(cwd,'Dektak','IND3D1.TXT'),
           1 : join(cwd,'Dektak','IND3D2.TXT'),
            2 : join(cwd,'Dektak','IND3D3.TXT'),
             3 : join(cwd,'Dektak','IND3D4.TXT')}
path13={0 : join(cwd,'Dektak','HB9M.TXT'),
           1 : join(cwd,'Dektak','HB9P.TXT')}
pathH9={0 : join(cwd,'Dektak','HB9M.TXT'),
           1 : join(cwd,'Dektak','HB9P.TXT')}
pathH10={0 : join(cwd,'Dektak','HB10M.TXT'),
           1 : join(cwd,'Dektak','HB10P.TXT')}
pathH11={0 : join(cwd,'Dektak','HB11M.TXT'),
           1 : join(cwd,'Dektak','HB11P.TXT')}
pathH12={0 : join(cwd,'Dektak','HB12M.TXT'),
           1 : join(cwd,'Dektak','HB12P.TXT')}

pathH13={0 : join(cwd,'Dektak','HB13M.TXT'),
           1 : join(cwd,'Dektak','HB13P.TXT')}
pathH14={0 : join(cwd,'Dektak','HB14M.TXT'),
           1 : join(cwd,'Dektak','HB14P.TXT')}
pathH15={0 : join(cwd,'Dektak','HB15M.TXT'),
           1 : join(cwd,'Dektak','HB15P.TXT')}
pathH16={0 : join(cwd,'Dektak','HB16M.TXT'),
           1 : join(cwd,'Dektak','HB16P.TXT')}
pathH17={0 : join(cwd,'Dektak','HB17M.TXT'),
           1 : join(cwd,'Dektak','HB17M.TXT')}





pathpvp={ 0: pathA,
      1: pathB,
      2: pathC,
      3: pathD,
      4 :pathE,
      5: path14,
      6: pathIND,
      7: path13}
#pathpvp={ 0: pathH9,
#      1: pathH10,
#      2: pathH11,
#      3: pathH12}
#
#pathpvp={ 0: pathH13,
#      1: pathH14,
#      2: pathH15,
#      3: pathH16,
#      4: pathH17}

masspvp={ 0: 1.69,
      1: 1.76,
      2: 1.86,
      3: 1.87,
      4 :1.88,
      5 : 1.88,
      6: 2.24,
      7: 2.2}

#masspvp={ 0: 1.69,
#      1: 1.76,
#      2: 1.86,
#      3: 1.87}
#masspvp={ 0: 1.69,
#      1: 1.76,
#      2: 1.86,
#      3: 1.87,
#      4: 1.87}

mass=1.69 #PVPA
mass=1.76 #PVPB
mass=1.86 #PVPC
mass=1.87 #PVPD
mass=1.88 #PVPE
#mass=1.88 # Mass Probe 13PVP die eigentlich 14 ist
#mass=2.2 # Masse der 13PVP Probe die bei 75RH lag
#mass=1.8
#mass=2.4
#mass=2.24 #Indo Film

#name="PVPFilm"
namepvp={ 0: "PVPFilmA",
      1: "PVPFilmB",
      2: "PVPFilmC",
      3: "PVPFilmD",
      4 :"PVPFilmE",
      5: "PVPFilm14",
      6: "IndFilm3",
      7: "PVPFilm13"}
#namepvp={ 0: "RITFilmA",
#      1: "RITFilmB",
#      2: "RITFilmC",
#      3: "RITFilmD"}
#
#namepvp={ 0: "RITFilmE",
#      1: "RITFilmF",
#      2: "RITFilmG",
#      3: "RITFilmH",
#      4: "RITFilmI"}
turns={0 : 0,
           1 : np.pi/2}

turns={0 : 0,
           1 : np.pi/2,
           2 : 1*np.pi/4,
           3 : 3*np.pi/4}
#M13NCD1.TXT PVP Film with strongly curved glass slide
#LEVLD1.TXT Indometacin Film with small to none curved glass slide
rhopvp=[]
for pp in range(len(pathpvp)):
    path1=pathpvp[pp]
    name=namepvp[pp]
    mass=masspvp[pp]
    Data=[]
    Len=[]
    for kk in range(len(path1)):
        Datared, Length= PAF.ReadData(path1[kk])
        Data.append(Datared) # Measurement of height [nm]
        Len.append(Length) # Measurement Scan Length [micrometer]
      
    
    rho=[]
    D=[]
    LD=[]
    RD=[]
    nR=[]
    rvec=[]
    shift1=300 #Do not consider the shift1 first and last Datapoints for jump detection
    shift2=40 #Do not regress the shift2 first Datapoints for the correction polynomial
    windowLength=300 #Only change if correction fails
    plotrawdata=1 # 1 Shows curve correction of the glass slide 
    for kk in range(len(path1)):
        rho1,D1,LD1,RD1,nR1,rvecR1=PAF.Dektak(Data[kk],Len[kk],mass,turns[kk],windowLength,plotrawdata,shift1,shift2,name)
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
    nPhi=100
    cmap='bone'
    dens=100
    PAF.ProfilometerPlot(Li,Ri,turns,nPhi,nAV,rvecm,cmap,dens,name)
    rhopvp.append(np.mean(rho))
rhoexp=np.mean(rhopvp)
rhostd=np.std(rhopvp)