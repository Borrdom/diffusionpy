import pandas as pd
from tkinter.filedialog import askopenfilename
import pandas.io.formats.excel
import numpy as np
import sys
import matplotlib.pyplot as plt
pd.options.mode.chained_assignment = None
pandas.io.formats.excel.header_style = None
original=askopenfilename()
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


try:
    indexofchange[1]
    MassDry1,TDry1=Mass[indexofchange[0]],Time[indexofchange[0]]
    MassDry2,TDry2=Mass[-1],Time[-1]
except:
    print("Are you sure there is a second drying step?")

x=np.hstack((TDry1,TDry2))
y=np.hstack((MassDry1,MassDry2))
fig,ax=plt.subplots()
ax.plot(Time,Mass,'k--',label="original data")
ax.plot(x,y,"r-",label="correction line")
CorrectionLine=lambda T:np.polyval(np.polyfit(x,y,1),T)
Difference=MassDry1-CorrectionLine(Time[indexofchange[0]:])
Mass[indexofchange[0]:]+=Difference

Datared2["m_Korr[mg]"]=Mass
Datared[Datanumberheader.columns[2]]=np.hstack(("m_Korr[mg]",Mass))


Time=Datared2["Time [minutes]"].values
Mass=Datared2["m_Korr[mg]"].values
ax.plot(Time,Mass,"k-",label="corrected data")
ax.legend()
DataNew=pd.concat([Dataupper,Datared])

filename=original.split(".")[-1]
newfilename=original.split(".")[0]+"_drift."+filename
DataNew.to_excel(newfilename,index=False)