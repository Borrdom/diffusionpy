import xloil as xlo
from .Stefan_Maxwell_segmental import Diffusion_MS,Diffusion1D
from .crank_and_other import crank,BHX
from .DasDennisSpacing import DasDennis
import numpy as np
import casadi as cs
import time
from tkinter import Tk,filedialog,simpledialog
from .read_componentdatabase import get_par
import xloil.pandas

@xlo.func
def Diffusion_MS_xloil(t:xlo.Array(float,dims=1),L:float,Dvec:xlo.Array(float,dims=1),w0:xlo.Array(float,dims=1),w8:xlo.Array(float,dims=1),Mi:xlo.Array(float,dims=1),
volatile:xlo.Array(bool,dims=1),full_output:bool,Gammai:xlo.Array(float,dims=2)=None,swelling:bool=False,taui:xlo.Array(float,dims=1)=None):   
    return Diffusion_MS(t,L,Dvec,w0,w8,Mi,volatile,full_output,Gammai,swelling,taui)

@xlo.func
def Diffusion1D_xloil(t:xlo.Array(float,dims=1),L0:float,Ds:float,ws0:float,ws8:float):
    return Diffusion1D(t,L0,Ds,ws0,ws8)
    
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
    pure=pure.reset_index(drop=True).reset_index()
    cell2.value=pure.fillna(0).replace("none",0).values
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
    xloil.Range("H2").formula="PC_SAFT_NpT(Parameters_PC_SAFT!$A$2:$Y$13;Parameters_PC_SAFT!$A$14:$F$59;A$1:H$1; A2:G2)"


_excelgui = xlo.ExcelGUI(ribbon=r'''
   <customUI xmlns="http://schemas.microsoft.com/office/2009/07/customui">
       <ribbon>
           <tabs>
               <tab id="customTab" label="THexel">
                   <group id="customGroup" label="MyButtons">
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