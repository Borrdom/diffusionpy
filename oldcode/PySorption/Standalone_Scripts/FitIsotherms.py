from ReadXandYdatafromExcel import sheet
from lmfit import minimize,Parameters
import numpy as np
import Isotherms as Is

#filename='C:/Users/domin/OneDrive/Promotion/Wissen/Woche47_2020/ChalykhPVPGesamt.xlsx'
#filename='C:/Users/domin/OneDrive/Promotion/Wissen/Unsortierte_Literatur_Daten/ChalykhPVPGesamt.xlsx'
filename="example_FitIsotherms.xlsx"
import matplotlib.pyplot as plt
sheet1=sheet(filename)


strx3="RH[-]"
stry3="ww[-]"

x3=sheet1.getrow(strx3)
ww=sheet1.getrow(stry3)

wp=1-ww
rhop=1180
rhow=997
y3=ww/rhow*(ww/rhow+wp/rhow)**-1






xhi=0.894
phiwinfty=0.3
k=10
aw=np.linspace(1E-13,1,100)
fit_params=Parameters()
fit_params.add("xhi",value=xhi,min=-6, max=6,vary=1)
fit_params.add("k",value=0,min=0, max=10,vary=1)
fit_params.add("H",value=0,min=0, max=10,vary=1)
fit_params.add("phi",value=0,min=0, max=10,vary=1)
phiwf=Is.flory2(aw,xhi)
phiwl=Is.langmuir2(aw,phiwinfty,k)
#plt.plot(aw,phiwf,"k--")
#plt.plot(aw,phiwl,"r--")
floryisotherm=lambda aw,fit_params: Is.flory2(aw,fit_params["xhi"].value)
langmuirisotherm=lambda aw,fit_params: Is.langmuir2(aw,fit_params["phi"].value,fit_params["k"].value)+Is.henry2(aw,fit_params["H"].value)
chalykhisotherm=lambda aw,fit_params: Is.flory2(aw,fit_params["xhi"].value)+Is.langmuir2(aw,fit_params["phi"].value,fit_params["k"].value)
henryflisotherm=lambda aw,fit_params: Is.henry2(aw,fit_params["H"].value)

def FitIsotherm(fit_params,phiwisotherm,aw,phiwexp):
    residual=lambda fit_params: phiwisotherm(aw,fit_params)-phiwexp  
    opt=minimize(residual,fit_params,method="least_squares",ftol=1E-4)#["x"]
    return opt.params
floryfit=FitIsotherm(fit_params,floryisotherm,x3,y3)
langmuirfit=FitIsotherm(fit_params,langmuirisotherm,x3,y3)
chalykhfit=FitIsotherm(fit_params,chalykhisotherm,x3,y3)
henryflfit=FitIsotherm(fit_params,henryflisotherm,x3,y3)
phiwfopt=Is.flory2(aw,floryfit["xhi"].value)
phiwlopt=Is.langmuir2(aw,langmuirfit["phi"].value,langmuirfit["k"].value)+Is.henry2(aw,henryflfit["H"].value)
phiwcopt=Is.flory2(aw,chalykhfit["xhi"].value)+Is.langmuir2(aw,chalykhfit["phi"].value,chalykhfit["k"].value)
phiwfhopt=Is.henry2(aw,henryflfit["H"].value)
fig,ax=plt.subplots()
ax.plot(aw,phiwfopt,"k-")
ax.plot(aw,phiwlopt,"r-")
ax.plot(aw,phiwfhopt,"g-")
ax.plot(x3,y3,"kx")
