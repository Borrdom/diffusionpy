import xloil as xlo
from .Stefan_Maxwell_segmental import Diffusion_MS,D_Matrix,Diffusion_MS_iter
from .crank_and_other import crank,BHX
from .DasDennisSpacing import DasDennis
import numpy as np
import casadi as cs
import time
from tkinter import Tk,filedialog,simpledialog
from .read_componentdatabase import get_par
import xloil.pandas
import pandas as pd
#from epcsaftpy import pcsaft,component,mixture
from .PyCSAFT_nue import lngi,vpure,dlnai_dlnxi
import fnmatch
# mix=component()+component()
# eos=pcsaft(mix)


@xlo.func
def Diffusion_MS_xloil(t:xlo.Array(float,dims=1),L:float,Dvec:xlo.Array(float,dims=1),w0:xlo.Array(float,dims=1),w8:xlo.Array(float,dims=1),Mi:xlo.Array(float,dims=1),
mobile:xlo.Array(bool,dims=1),full_output:bool=False,dlnai_dlnwi:xlo.Array(float,dims=2)=None,swelling:bool=False,taui:xlo.Array(float,dims=1)=None,rho0i:xlo.Array(float,dims=1)=None):   
    return Diffusion_MS(t.copy(),L,Dvec,w0.copy(),w8.copy(),Mi,mobile,full_output,dlnai_dlnwi,swelling,taui,rho0i)
@xlo.func
def Diffusion_MS_iter_xloil(t:xlo.Array(float,dims=1),L:float,Dvec:xlo.Array(float,dims=1),w0:xlo.Array(float,dims=1),w8:xlo.Array(float,dims=1),Mi:xlo.Array(float,dims=1),
mobile:xlo.Array(bool,dims=1),full_output:bool=False,swelling:bool=False,taui:xlo.Array(float,dims=1)=None,rho0i:xlo.Array(float,dims=1)=None,
pure:xlo.Array(object,dims=2)=np.asarray([[]]),kij:xlo.Array(object,dims=2)=np.asarray([[]])):   
    T=298.15
    Mw=pure[1:,3].astype(float)
    mi=pure[1:,4].astype(float)
    sigi=pure[1:,6].astype(float)
    ui=pure[1:,7].astype(float)
    eAiBi=pure[1:,8].astype(float)
    kAiBi=pure[1:,9].astype(float)
    Na=pure[1:,10].astype(float)
    nc=len(Mw)
    kij1=np.char.replace(kij[1:,2].astype(str),",",".").astype(float)
    kij=D_Matrix(kij1,nc)
    vpures=vpure(1E5,T,mi,sigi,ui,eAiBi,kAiBi,Na)
    par={"mi":mi,
    "Mi":Mw,
    "si":sigi,
    "ui":ui,
    "eAi":eAiBi,
    "kAi":kAiBi,
    "NAi":Na,
    "kij":kij,
    "vpure":vpures}
    return Diffusion_MS_iter(t.copy(),L,Dvec,w0.copy(),w8.copy(),Mi,mobile,full_output,swelling,taui,rho0i,T,par)

@xlo.func
def gradient(x:xlo.Array(float,dims=1),y:xlo.Array(float,dims=1)):
    return np.gradient(x,y)

@xlo.func
def reduce_points(x,n:int):
    return x[::n]
@xlo.func
def crank_xl(t,L0,Ds,ws0,ws8):
    Xs0=ws0/(1-ws0)
    Xs8=ws8/(1-ws8)
    mt08=crank(t,Ds/L0**2,1)
    Xs=mt08*(Xs8-Xs0)+Xs0
    ws=Xs/(1+Xs)
    return ws
@xlo.func
def interp1d(x,xp,fp):
    return np.interp(x.flatten(),xp.flatten(),fp.flatten())

@xlo.func
def BHX_xloil(t:xlo.Array(float,dims=1),kf:float,kr:float,ws0:float,ws8:float,mfinfty:float,mrinfty:float):
    Xs0=ws0/(1-ws0)
    Xs8=ws8/(1-ws8)
    Xs=BHX(t,kf,np.asarray([kr]),1,Xs0,Xs8,mfinfty,np.asarray([mrinfty]))
    ws=Xs/(1+Xs)
    return ws

@xlo.func
def DasDennis_xloil(p:int,dim:int):
    return DasDennis(p,dim).T

@xlo.func
def add_custom(a,b):
    return a+b

#@xlo.func(command=True)
def get_path(ctrl):
    root = Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename()
    components= simpledialog.askstring("THexel", "Enter component names: Seperate each name via a comma").split(",")
    get_par_xloil(components,file_path)
    root.update()
    root.destroy()



#@xlo.func
def get_par_xloil(subst_input,path):
    wb=xlo.active_workbook()
    sheet1=wb.add(name="Parameters_PC_SAFT")
    xloil.Range("Parameters_PC_SAFT!A1").value="PC-SAFT parameters"
    xloil.Range("Parameters_PC_SAFT!A2").value="Number"

    pure,kij=get_par(subst_input,path)
    cell1=xloil.Range("Parameters_PC_SAFT!B2:Z2")
    cell2=xloil.Range("Parameters_PC_SAFT!A3:Z"+str(2+len(pure)))
    cell3=xloil.Range("Parameters_PC_SAFT!A14:F14")
    cell4=xloil.Range("Parameters_PC_SAFT!A15:F"+str(14+(len(pure)-1)*len(pure)//2))
    header1=["Substance","CAS-Nr","Molar Mass",	"Segment Number","Segment Number/Molar Mass","Segment Diameter","Energy Parameter",
             	"Associating Energy", "Associating Volume",	 "Associating Scheme","Associating Scheme",	"Dipole Moment",
                  "Polarizability", "Quadrupole Moment" ,"Polarizability", "Charge", "z_alpha",	"B_alpha_alpha",
                      	"T_melt",	"H_melt",	"Delta cp", "Delta cp slope","Tg","Free Param","Publication 1st"]
    header2=["Substance I", "Substance II","kij","kij T","kij T","kij T"]
    pure["Delta cp slope"]=pure["Molar Mass"].astype(float)*0
    pure["Free Param"]=pure["Molar Mass"].astype(float)*0

    pure['Associating Scheme'] = pure['Associating Scheme'].replace("2B","1/1").str.extract('(\d+)', expand=False)
    pure=pure[header1]
    kij["kij T"]=pure["Molar Mass"].astype(float)*0
    kij=kij[header2] # this datfarme containes the information about the kij 
    cell1.value=pure.columns[:,None].T
    pure=pure.reset_index(drop=True).reset_index().fillna(0).replace("none",0)
    
    mi=pure["Segment Number"].values.astype(float)
    sigi=pure["Segment Diameter"].values.astype(float)
    ui=pure["Energy Parameter"].values.astype(float)
    eAiBi=pure["Associating Energy"].values.astype(float)
    kAiBi=pure["Associating Volume"].values.astype(float)
    Na=pure["Associating Scheme"].T.drop_duplicates().T.values.astype(float).flatten()
    Mw=pure["Molar Mass"].values.astype(float).flatten()
    mi[mi==0.]=(pure["Segment Number/Molar Mass"].values.astype(float).flatten()*Mw)[mi==0.]
    pure["Segment Number"]=mi
    pure["Segment Number/Molar Mass"]=mi/Mw
    vpures=vpure(1E5,298.15,mi,sigi,ui,eAiBi,kAiBi,Na)
    rho0=(vpures/Mw*1000.)**-1
    pure["Free Param"]=rho0
    cell2.value=pure.values
    cell3.value=kij.columns[:,None].T
    cell4.value=kij.fillna(0).replace("none",0)[header2].values
    
    
    #cell1.color=(118,224,280)
    #cell4.color=(118,224,280)
    return 



#@xlo.func(command=True)
def create_header(ctrl):
    cell1=xloil.Range("Parameters_PC_SAFT!B3:B13")
    comps=cell1.trim().value.flatten()
    x=np.ones_like(comps).astype(float)
    x=x/x.shape[0]
    xloil.Range("A1:H1").value=np.asarray(["Fraction type"]+list(comps)+["Pressure [bar]","Temp. [K]","State (l/v) [-]","Density [kg/m³]"])[:,None].T
    xloil.Range("A2:G2").value=np.asarray(["x"]+list(x)+[1,298.15,"l"])[:,None].T
    xloil.Range("H2").formula="PC_SAFT_NpT(Parameters_PC_SAFT!$A$2:$Y$13;Parameters_PC_SAFT!$A$14:$F$59;A$1:Z$1; A2:G2)"

try:
    _excelgui = xlo.ExcelGUI(ribbon=r'''
    <customUI xmlns="http://schemas.microsoft.com/office/2009/07/customui">
        <ribbon>
            <tabs>
                <tab id="customTab" label="THexel">
                    <group id="customGroup" label="Starters">
                        <button id="Get_Parameter" imageMso="ImexRunImport" size="large" onAction="pressParameter" />
                        <button id="Get_Header" imageMso="PivotPlusMinusFieldHeadersShowHide" size="large" onAction="pressHeader" />
                    </group>
                </tab>
            </tabs>
        </ribbon>
    </customUI>
    ''', 
    funcmap={
        'pressParameter': get_path,
        'pressHeader': create_header
        })
except:
    pass

@xlo.func
def PC_SAFT_NpT2(pure,kij,header,inputs):
    # pure=pd.DataFrame(pure[1:,:],columns=pure[0,:])
    # kij=pd.DataFrame(kij[1:,:],columns=kij[0,:])
    # name=pure[1:,1].astype(str)
    Mw=pure[1:,3].astype(float)
    mi=pure[1:,4].astype(float)
    # mi_Mw=pure[1:,5].astype(float)
    sigi=pure[1:,6].astype(float)
    ui=pure[1:,7].astype(float)
    eAiBi=pure[1:,8].astype(float)
    kAiBi=pure[1:,9].astype(float)
    Na=pure[1:,10].astype(float)
    # Nd=pure[1:,11].astype(float)
    # dipol1=pure[1:,12].astype(float)
    # dipol2=pure[1:,13].astype(float)
    # quadro1=pure[1:,14].astype(float)
    # quadro2=pure[1:,15].astype(float)
    # charge1=pure[1:,16].astype(float)
    # charge2=pure[1:,17].astype(float)
    # charge3=pure[1:,18].astype(float)
    # TSL=pure[1:,19].astype(float)
    # HSL=pure[1:,20].astype(float)
    # dcp=pure[1:,21].astype(float)
    # dcpT=pure[1:,22].astype(float)
    # Tg=pure[1:,23].astype(float)
    free=pure[1:,24].astype(float)

    nc=len(Mw)
    #kij1=np.zeros((nc,nc))
    #kij1[np.triu_indices(nc,k=1)]=
    kij1=np.char.replace(kij[1:,2].astype(str),",",".").astype(float)
    fractiontype=inputs[0,0]
    xi=inputs[:,1:1+nc].astype(float).flatten()
    T=float(inputs[0,nc+2])
    p=float(inputs[0,nc+1])*1E5
    state=inputs[0,nc+3]
    kij=D_Matrix(kij1,nc)
    if np.any(free==0.): 
        vpures=vpure(p,T,mi,sigi,ui,eAiBi,kAiBi,Na)
    else:
        vpures=(free/Mw*1000.)**-1
    # if True:#nc!=eos.mixture.nc:
    #     a=[]
    #     for i in range(nc):
    #             a.append(component(name=name[i],ms=mi[i],Mw=Mw[i],sigma=sigi[i],eps=ui[i],eAB=eAiBi[i],kappaAB=kAiBi[i],sites=[0,Na[i],Nd[i]]))
    #             if i>0: 
    #                 pars+=a[i]
    #             else:
    #                 pars=a[i]
    #     eos = pcsaft(pars)



    #eos1.KIJ0saft=kij1
    results=np.asarray([])
    #rho0,Xass0=eos.density_aux(xi,eos.temperature_aux(T),p,state=state.upper())


    
    def add_entry(results,new):
        return np.hstack((results,new))
    def generate(var):
        return (val for val in var)

    if fractiontype=="w": Mi=Mw
    if fractiontype=="x": Mi=None
    Mw/np.sum(xi*Mw)
    for entry in header[0]:
        if "Density" in entry:
            if fractiontype=="x":
                rho0=(np.sum(vpures*xi))**-1
                results=add_entry(results,rho0*(xi*Mw/1000).sum())
            elif fractiontype=="w":
                v0i=vpures/Mw*1000
                rho0=(np.sum(v0i*xi))**-1
                results=add_entry(results,rho0)
        elif "Gamma" in entry:
            lnwx=np.log(Mi*np.sum(xi/Mi)) if fractiontype=="w" else 0.
            lngammai=generate(lngi(T,xi,mi,sigi,ui,eAiBi,kAiBi,Na,vpures,Mi=Mi,kij=kij).flatten()+lnwx) if 'lngammai' not in vars() else lngammai
            #lnphi=generate(eos.logfugef(xi,T,p,state=state.upper(),v0=1/rho0,Xass0=Xass0)[0]) if 'lnphi' not in vars() else lnphi
            results=add_entry(results,lngammai.send(None))
        elif "ln(a)" in entry:#fnmatch.fnmatchcase(entry.replace(" ", "").replace("[-]",""),"ln(a)?"):
            lnactivity=generate(lngi(T,xi,mi,sigi,ui,eAiBi,kAiBi,Na,vpures,Mi=Mi,kij=kij).flatten()+np.log(xi)) if 'lnactivity' not in vars() else lnactivity
            #lnphi=generate(eos.logfugef(xi,T,p,state=state.upper(),v0=1/rho0,Xass0=Xass0)[0]) if 'lnphi' not in vars() else lnphi
            results=add_entry(results,lnactivity.send(None))
        elif "Mass fraction" in entry:
            if fractiontype=="x":
                fracw=generate(xi/Mw/(xi/Mw).sum()) if 'fracw' not in vars() else fracw
            elif fractiontype=="w":
                fracw=generate(xi) if 'fracw' not in vars() else fracw
            results=add_entry(results,fracw.send(None))
        elif "Mole fraction" in entry:
            if fractiontype=="x":
                fracx=generate(xi) if 'fracx' not in vars() else fracx
            elif fractiontype=="w":
                fracx=generate(xi*Mw/(xi*Mw).sum()) if 'fracx' not in vars() else fracx
            results=add_entry(results,fracx.send(None))
        elif "di" in entry:
            THFaktor=generate(dlnai_dlnxi(T,xi,mi,sigi,ui,eAiBi,kAiBi,Na,vpures,Mi=Mi,kij=kij).flatten()) if 'THFaktor' not in vars() else THFaktor
            #lnphi=generate(eos.logfugef(xi,T,p,state=state.upper(),v0=1/rho0,Xass0=Xass0)[0]) if 'lnphi' not in vars() else lnphi
            results=add_entry(results,THFaktor.send(None))
        elif "M [" in entry:
            if fractiontype=="x":
                results=add_entry(results,(xi*Mw).sum())
            elif fractiontype=="w":
                results=add_entry(results,1/(xi/Mw).sum())
        elif "Z [" in entry:
            pass
            #results=add_entry(results,eos.Z(xi,T,p,state=state.upper()))
        elif "gres [" in entry:
            pass
        elif "hres [" in entry:
            pass
        elif "cpres [" in entry:
            pass
        elif "Stability" in entry:
            pass
        elif "hE" in entry:
            pass
        elif "gE" in entry:
            pass
        elif "ln(f)" in entry:
            pass
        elif "ln(phi)" in entry:
            results=add_entry(results,np.asarray([1.]))
        elif "ln(a)" in entry:
            pass
        elif "ln(Gamma)" in entry:
            pass
        elif "pmV" in entry:
            pass
        elif "pLV" in entry:
            pass
        elif "rhoLV" in entry:
            pass
        elif "Monomer frac." in entry:
            pass
        else:
            pass


    return results[:,None].T

