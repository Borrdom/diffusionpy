import pandas as pd
from os.path import join
from os import getcwd
from tkinter.filedialog import askopenfilename
from ReadXandYdatafromExcel import sheet

def ExtractIsotherm(molecule1="pvpva64",molecule2="naproxen",RH=0.9,DL=0.7,T=25):
    cwd=getcwd()
    RHstr=str(max(int(RH*100),RH*100))+"%"
    Tstr=str(max(int(T),T))+"°C"
    simpath=join(cwd,"FeedWasser",RHstr+"_"+Tstr,"_".join(["EXP",RHstr,molecule1,molecule2,Tstr]))
    SAFT=pd.read_excel(simpath+".xlsx")
    wAPISIMraw=SAFT["wAPI[-]"].values[::-1]
    wwSIMraw=SAFT["ww[-]"].values[::-1]
    return wAPISIMraw,wwSIMraw
def ExtractSorption(filename):
    sheet1=sheet(filename)
    return sheet1.getrows(*tuple(["t[min]","w[-]","tK[min]","wK[-]","tges[min]","wges[-]","M0[mg]"]))

if __name__=="__main__":
    wAPISIMraw,wwSIMraw=ExtractIsotherm(molecule1="pvpva64",molecule2="naproxen",RH=0.9,DL=0.7,T=25)
    filename=askopenfilename()
    t1,w1,t2,w2,t3,w3,m0=ExtractSorption(filename)