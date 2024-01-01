import pandas as pd
from os.path import join
from os import getcwd
from PySorption.DiffusionFinal import Film
from PySorption.dcc import window




def ExtractIsotherm(RH=0.9,DL=0.68,T=25):
    cwd=getcwd()
    molecule1="pvpva64"
    molecule2="naproxen"
    print(molecule1)
    print(molecule2)
    RHstr=str(max(int(RH*100),RH*100))+"%"
    Tstr=str(max(int(T),T))+"°C"
    simpath=join(cwd,"LeastFit","FeedWasser",RHstr+"_"+Tstr,"_".join(["EXP",RHstr,molecule1,molecule2,Tstr]))
    #####
    SAFT=pd.read_excel(simpath+".xlsx")
    #EXP=pd.read_excel(exppath+".xlsx")
    wAPISIMraw=SAFT["wAPI[-]"].values[::-1]
    wwSIMraw=SAFT["ww[-]"].values[::-1]
    return wAPISIMraw,wwSIMraw
def SorptionAndKris(RH=0.9,DL=0.68,T=25):
    wAPISIMraw,wwSIMraw=ExtractIsotherm(RH,DL,T)
    alpharaw=(DL-wAPISIMraw)/(1-wAPISIMraw)
    window1=window()
    Data=window1.result
    Film1=Film(window1.result,window1.components,window1.wPolyASD)# if len(window1.components)==2 else ASDFilm(window1.result,*window1.components)
    molecule1=window1.components[0]["name"]
    molecule2=window1.components[2]["name"]

    rho1=window1.components[0]["rho0"]
    rho2=window1.components[2]["rho0"]
    Film1.alphavec=alpharaw
    #Film 90 oder Film 87
    wPolyASD=1-DL
    wAPIASD=DL
    rhoASD=(wPolyASD/rho1+wAPIASD/rho2)**-1

    Film1.rho2GGWvec=wwSIMraw/(1-wwSIMraw)*rhoASD
    Film1.wwSIMraw=wwSIMraw
    Film1.wAPISIMraw=wAPISIMraw
    Film1.RHraw=RH
    Film1.endpoint=False
    Film1.Mode="Kris"
    Film1.ideal=0
    Film1.SAFT=True
    Film1.NET=True
    Film1.exp_boundary()
    Film1.Origin=True
###

    sigma=1.97730286e-02
    DAPI=6.59010330e-17
    #DAPI=3.59010330e-17 #best?
    #DAPI=6.59010330e-20
    g=3.2
    kt=5.0778700e-12


    sigma=1.97730286e-02
    DAPI=2.59010330e-16
    #DAPI=6.59010330e-20
    g=2.8
    kt=7.0778700e-12

    sigma=1.97730286e-02
    DAPI=8.59010330e-17
    #DAPI=6.59010330e-20
    g=2.8
    kt=9.0778700e-12
##



    sigma=1.97730286e-02
    DAPI=6.59010330e-17
    #DAPI=3.59010330e-17 #best?
    #DAPI=6.59010330e-20
    g=3.2
    kt=5.0778700e-12
    #kt=2.0778700e-13
    Film1.solubility=0.03*(1-0.24968)
    Film1.kt=kt
    Film1.g=g
    Film1.DAPI=DAPI
    Film1.sigma=sigma

    #Maxwell
    #Film1.D12=8E-13 #Dünn Dünn
    Film1.D12=3E-13 #Dünn
    #Film1.D12=3E-16 #Dünn
    Film1.D12=2E-13 #Dick


    Film1.Mode="NEAT"
    #Film1.D12=Film1.CalcTerDiffcoeff() # now the  diffusion coefficient of naproxen is not available so self.D12 will be inserted
    Film1.Mode="Kris"
    

    # In the next step D12 is the new diffcoeff and now is the ternary diffusioncoefficent

    #Film1.kt=3E-5

    #Film1.D12=0.8E-13
    #Film1.kt=0.9E-5 #Film 90 fit #Prediction von Film 87 möglich # Film 76 ok #Film 71 ok #Film 66 ok #Film 55 exzellent
    #Film1.kt=0.7E-5
    Film1.wPolyASD=1-DL
    Film1.Diffusion1D()

    return Film1.D12
if __name__=="__main__":
    RH=0.9
    DL=0.68
    T=25
    D12=SorptionAndKris(RH=RH,DL=DL,T=T) # Aufräumen Dl=0
