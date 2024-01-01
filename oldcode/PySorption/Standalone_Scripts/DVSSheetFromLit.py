import pandas as pd
import os
import pandas.io.formats.excel
from tkinter.filedialog import askopenfilename
import interparc
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
Litname=askopenfilename()
Lit=pd.read_excel(Litname)
t=Lit["t/min"].dropna().values
w=Lit["w/-"].dropna().values
rho=Lit["rho/kg/m^3"].dropna().values
L=Lit["L0/mum"].dropna().values*1E-6
RH0=Lit["RH0/-"].dropna().values
RHINF=Lit["RHINF/-"].dropna().values
t12parc,wparc=interparc.interparc(t**(1/2),w,300)
#t12parc=np.linspace(t[0]**(1/2),t[-1]**(1/2),300)
tparc=t12parc**2
#wparc=InterpolatedUnivariateSpline(t,w)(tparc)

d0=0.0145
A0=np.pi/4*d0**2
m0=A0*rho*L*1E6
mparc=m0/(1-wparc)

pandas.io.formats.excel.header_style = None
cwd=os.getcwd()
TemplateName=os.path.join(cwd,"DVSSheetFromLitTemplate.xls")
Template=pd.read_excel(TemplateName)


Templateheader=Template[Template['DVS-Advantage-Plus-Data-File'].str.contains("minutes").fillna(False)]
idx=Templateheader.index[0]
Dataupper=Template.loc[:idx]
Datared2=Template.loc[idx+1:]
Datared=Template.loc[idx+1:]
Datared2.columns=Template.loc[idx:].iloc[0]

def Drying(RH1,RH0,m0,t,m):
    def mdummy(m0,m1,t):
            return m1*(1-np.exp(-t/t[-1]*10))+m0*(np.exp(-t/t[-1]*10))
    m1=m[0]
    m2=m[-1]
    tdry=np.linspace(0,300,31)
    adry=np.zeros(31)
    mdry=mdummy(m1,m0,tdry)
    mdry[-1]=m0
    tstep=tdry
    astep=np.ones(31)*RH0
    mstep=mdummy(m0,m1,tstep)
    mstep[-1]=m1
    if RH0!=0 :
        t=np.hstack((tdry,tstep+310,t+620))
        m=np.hstack((mdry,mstep,m))
        a=np.hstack((adry,astep,np.ones_like(t)*RH1))*100
    else:
        t=np.hstack((tdry,t+310))
        m=np.hstack((mdry,m))
        a=np.hstack((adry,np.ones_like(t)*RH1))*100
    return t,m,a
texp,mexp,aexp=Drying(RHINF,RH0,m0,tparc,mparc)


#append row to the dataframe
for i,val in enumerate(texp):
    new_row = {'DVS-Advantage-Plus-Data-File':texp[i], 'Unnamed: 1':mexp[i], 'Unnamed: 9':aexp[i]}
    Datared = Datared.append(new_row,ignore_index=True)





DataNew=pd.concat([Dataupper,Datared])
#DataNew=pd.concat([pd.DataFrame(DataNew.columns),DataNew],ignore_index=True,axis=1)
filename=Litname.split(".")[-1]
newfilename=Litname.split(".")[0]+"_DVS."+filename
DataNew.to_excel(newfilename,index=False)