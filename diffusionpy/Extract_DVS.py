import os
import time
import numpy as np
from tkinter.filedialog import askopenfilenames
import xlrd
from scipy.ndimage import median_filter
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
import matplotlib
import tkinter as tk
from functools import partial
import matplotlib
# matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
cwd=os.getcwd()
import pandas as pd
from collections import OrderedDict
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import interp1d
import matplotlib as mpl
from scipy.optimize import curve_fit
from numba import njit

import xlwings as xw
# from xlutils.copy import copy as xl_copy


def averagelist(lis):
    ave=sum(lis)/len(lis)
    std=(sum([(ave-val)**2 for val in lis])/(len(lis)-1))**0.5
    return ave,std

def averagelisttimeseries(t,x):
    tarr=np.array(t)
    xarr=np.array(x)
    nt=tarr.shape[2]
    for i in range(nt):
        tarr[:,:,i]=tarr[:,:,i]-tarr[:,:,0]
    tendn=tarr[:,:,-1]
    amin=np.argmin(tendn,axis=0)
    for i in range(len(t)):
        for j in range(len(t[i])):
            x_fun=interp1d(tarr[i,j,:],xarr[i,j,:])
            xarr[i,j,:]=x_fun(tarr[amin[j],j,:])
    xave=np.average(xarr,axis=0)
    xstd=np.std(xarr,axis=0)
    return xave,xstd

def averagedict(x):
    df= pd.DataFrame(x,dtype=object)
    listdf=[df[i].dropna().tolist() for i in df.keys()]
    mean=np.asarray([np.mean(i, axis=0) for i in listdf])
    std=np.asarray([np.std(i, axis=0) for i in listdf])
    return mean, std

def unnesting(df, explode, axis):
    if axis==1:
        idx = df.index.repeat(df[explode[0]].str.len())
        df1 = pd.concat([
            pd.DataFrame({x: np.concatenate(df[x].values)}) for x in explode], axis=1)
        df1.index = idx

        return df1.join(df.drop(explode, axis=1), how='left')
    else :
        df1 = pd.concat([pd.DataFrame(df[x].tolist(), index=df.index).add_prefix(x) for x in explode], axis=1)
        return df1.join(df.drop(explode, axis=1), how='left')
#from sklearn.gaussian_process import GaussianProcessRegressor
#Least Square Regression for data provided by lis




class Measurement:
    """ The class Measurement inherits methods to read and process the DVS excel files from the Apparatus of Surface Measurement Systems.
    The endpoints of each RH step and their kinetics are extracted and the crank equation is fitted to the kinetics.
    The extracted data is added via two new sheets to all selected .xls files. The reading functions are specialized for the format of the .xls files. 
    There is little robustness so large changes cause the methids not to work . 
    The class accepts more than one excel file and then performs averaging of the kinetics. For that purpose, it is assumed that the all excel sheets are 
    replicates of each other, meaning that their RH steps must match 
    """
    def __init__(self,DATA_FILE,desorption=False,actual=False,tend=-1,arc=True,nt=15):
        """Initialization"""
        self.DATA_FILE=DATA_FILE  #if isinstance(DATA_FILE, str) else DATA_FILE[0]#join("Pure",self.molecule,DATA_FILE) if isinstance(DATA_FILE, str) else join("Pure",self.molecule,DATA_FILE[0])
        self.nt=nt
        self.tend=tend
        self.arc=arc
        self.des=desorption
        self.actual=actual
        self.indJ=0
        for val in self.DATA_FILE:
            self.read_excel_file(val)
        self.FitDiffusionData()
    @staticmethod
    def FasterXlrdRead(filename):
        book = xlrd.open_workbook(filename, encoding_override = "utf-8")
        sheet = book.sheet_by_name("DVS Data")
        nstart=[i for i,val in enumerate(sheet.col_values(0)) if "Time" in str(val)][0]+1
        nrow=sheet.nrows
        ncol=sheet.ncols
        nrow2=nrow-nstart
        statusnames=[sheet.cell(i, 0).value for i in range(35)]
        status=[sheet.cell(i, 1).value for i in range(35)]
        statusinfo=dict(zip(statusnames,status))
        header=np.asarray([sheet.cell(nstart-1, i).value for i in range(ncol)])
        idxtime=[i for i,val in enumerate(header) if "min" in val]
        idxmg=[i for i,val in enumerate(header) if "mg" in val]
        idxmass=idxmg[1] if len(idxmg)>1 else idxmg[0]
        idxPartialA=[i for i,val in enumerate(header) if "Partial Pressure (Solvent A)" in val]
        idxPartialB=[i for i,val in enumerate(header) if "Partial Pressure (Solvent B)" in val]
        idxTargetA=idxPartialA[0]
        idxActualA=idxPartialA[1]
        idxTargetB=idxPartialB[0]
        idxActualB=idxPartialB[1]
        def get_col(idx):
            col1=np.asarray(sheet.col_values(idx)[nstart:])
            # col1[col1=='']='0'
            return np.asarray(col1,dtype=float)
        time=get_col(idxtime[0])
        mass=get_col(idxmass)
        TargetA=get_col(idxTargetA)
        ActualA=get_col(idxActualA)
        TargetB=get_col(idxTargetB)
        ActualB=get_col(idxActualB)
        return time,mass,TargetA,ActualA,TargetB,ActualB,statusinfo

    def read_excel_file(self,filename):
        nt=self.nt

        
        
        # a2_exp,a1_exp,a_exp,t_exp,m_exp=self.FasterXlrdRead(filename) 
        t_exp,m_exp,TargetA,ActualA,TargetB,ActualB,statusinfo=self.FasterXlrdRead(filename)
        if statusinfo["Solvent A:"]:
            a_exp=TargetA
            aa_exp=ActualA
        elif statusinfo["Solvent B:"]:
            a_exp=TargetB
            aa_exp=ActualB

        # aa_exp=a2_exp if (a1_exp==np.zeros_like(a1_exp)).all() else a1_exp

        m_exp=m_exp*1E-6
        t_exp=t_exp*60

        Feuchte, indices=np.sort(np.unique(a_exp,return_index=True))

        dRHdt=np.diff(a_exp)
        indch=np.hstack((1,np.where(dRHdt!=0)[0]+1))
        des=dRHdt[indch-1]<=0
        abso=dRHdt[indch-1]>=0
        nJump=indch.shape[0]
        idxdes=np.where(des)[0]
        indices=np.append(indch,-2) #Anfügen des Endes der Messung
        mnull=m_exp[indices[1]] #Masse der Probe vor dem Experiment [mg]
        w_exp=(m_exp-mnull)/m_exp
        self.t_exp=t_exp
        self.w_exp=w_exp
        wreal=[]
        wanf=[]
        wrealmittel=[]
        Feuchtereal=[]
        mlist=[]
        tlist=[]
        wlist=[]
        key=[]
        wrealdes=[]
        wanfdes=[]
        wrealmitteldes=[]
        Feuchterealdes=[]
        mlistdes=[]
        tlistdes=[]
        wlistdes=[]
        keydes=[]
        Feuchteactual=[]
        Feuchteactualdes=[]
        Feuchteanf=[]
        Feuchteanfdes=[]
        sratel=[]
        srateldes=[]
        t0list=[]
        t0listdes=[]

        def Appending(wreal,wanf,wrealmittel,Feuchtereal,Feuchteactual,Feuchteanf,mlist,tlist,wlist,sratel,t0list,key,indices,k,des=False):
                jump=str(a_exp[indices[k]-1])+"_"+str(a_exp[indices[k]]) #if des==False else str(a_exp[indices[k]])+"_"+str(a_exp[indices[k]-1])
                key.append(jump)
                toffset=0 
                mapp=m_exp[indices[k]+toffset:indices[k+1]]
                mapp=mapp 
                tapp=t_exp[indices[k]+toffset:indices[k+1]]
                wapp=(mapp-mnull)/mapp 
                mmax=np.max(mapp)
                mpercent=(mapp-mapp[0])/(mmax-mapp[0]) 
                srate=np.gradient(wapp,tapp)
                srate_smooth=median_filter(srate,50,mode="nearest")
                t0list.append(tapp)
                tapp=(tapp-tapp[0])
                
                if self.tend==-1:
                    indmmax=self.tend
                else:
                    indmmax=np.where(tapp>(self.tend*60))[0][0]
                mapp=mapp[0:indmmax]
                tapp=tapp[0:indmmax]
                wapp=wapp[0:indmmax]
                srate_smooth=srate_smooth[0:indmmax]

                if self.arc:

                    tapp2, mapp2=self.interparc(tapp, mapp, nt+1)
                else:

                    tapp2=np.linspace(0,tapp[-1],nt+1)
                    mapp2=np.interp(tapp2,tapp,mapp) 

                wapp2=(mapp2-mnull)/mapp2 #if des==False else (mapp2+mnull)/mapp2
                wapp2=np.fmax(wapp2,1E-4)
                mapp2=mapp2 #if des==False else -mapp2
                tlist.append(tapp2)
                mlist.append(mapp2)
                wlist.append(wapp2) #if des==False else wlist.append(-wapp2+wapp2[0]+wapp2[-1])
                wanf.append(wapp2[0]) #if des==False else wanf.append(wapp2[-1])
                wreal.append(wapp2[-1]) #if des==False else wreal.append(wapp2[0])
                srate2=np.gradient(wapp2,tapp2)
                srate2_smooth=median_filter(srate2,50,mode="nearest")

                wrealmittel.append(0.7*wapp2[-1]+0.3*wapp2[0])
                FFeuchte=np.mean(aa_exp[indices[k]+toffset:indices[k+1]-1])*0.01 #if des==False else np.mean(aa_exp[indices[k-1]+toffset:indices[k]-1])*0.01
                Feuchtean=a_exp[indices[k-1]]*0.01
                Feuchte=a_exp[indices[k]]*0.01 #if des==False else a_exp[indices[k-1]]*0.01
                Feuchtereal.append(Feuchte)
                Feuchteanf.append(Feuchtean)
                Feuchteactual.append(FFeuchte)
                sratel.append(srate2_smooth)
                #plt.plot(tapp,mapp)
                return
        for k in range(1,nJump):
            if all(k!=idxdes):
                Appending(wreal,wanf,wrealmittel,Feuchtereal,Feuchteactual,Feuchteanf,mlist,tlist,wlist,sratel,t0list,key,indices,k)
            elif any(k==idxdes) and len(m_exp[indices[k]:indices[k+1]])>0:
                Appending(wrealdes,wanfdes,wrealmitteldes,Feuchterealdes,Feuchteactualdes,Feuchteanfdes,mlistdes,tlistdes,wlistdes,srateldes,t0listdes,keydes,indices,k,des=True)
        self.nHum=len(wreal) if not self.des==True else len(wrealdes)
        #self.nHumdes=len(wrealdes)
        self.dickeys=key if not self.des==True else keydes
        #endhums=[val.split("_")[1] for i,val in enumerate(self.dictkeys)]
        #self.wrealdes=wrealdes#if not self.des==True else keydes
        #self.wanfdes=wanfdes
        self.wrealmittelref=wrealmittel if not self.des==True else wrealmitteldes
        #self.wrealmitteldes=wrealmitteldes
        #self.Feuchterealdes=Feuchterealdes#[1:]
        self.tlist=tlist if not self.des==True else tlistdes
        self.mlist=mlist if not self.des==True else mlistdes
        self.wlist=wlist if not self.des==True else wlistdes
        self.t0list=t0list if not self.des==True else t0listdes
        self.wreal=np.asarray(wreal) if not self.des==True else np.asarray(wrealdes)
        self.wanf=np.asarray(wanf) if not self.des==True else np.asarray(wanfdes)
        self.Feuchtereal=np.asarray(Feuchtereal) if not self.des==True else np.asarray(Feuchterealdes)
        self.Feuchteactual=np.asarray(Feuchteactual) if not self.des==True else np.asarray(Feuchteactualdes)
        self.mnull=mnull
        self.Feuchteanf=np.asarray(Feuchteanf) if not self.des==True else np.asarray(Feuchteanfdes)
        self.sratel=sratel if not self.des==True else srateldes
        Feuchteges=np.hstack((Feuchtereal,Feuchterealdes))
        wges=np.hstack((wreal,wrealdes))
        wges=wges[Feuchteges!=0]
        Feuchteges=Feuchteges[Feuchteges!=0]


        self.Feuchteiso,idxunique=np.unique(Feuchteges,True)
        self.wiso=np.asarray([np.average(wges[Feuchteges==val]) for i,val in enumerate(self.Feuchteiso)])

        # If Absorbtion and Desorption have differences do not take the average between them even if its the same RH
        Feuchteun,idxun=np.unique(self.Feuchtereal,True)
        wun=np.asarray([np.average(np.asarray(self.wreal)[self.Feuchtereal==val]) for i,val in enumerate(Feuchteun)])
        self.wiso=np.asarray([np.average(wun[Feuchteun==val]) if len(wun[Feuchteun==val])>0 else np.average(wges[Feuchteges==val]) for i,val in enumerate(self.Feuchteiso)])

        # If the actual RH is required we take it instead
        Feuchteactualges=np.hstack((Feuchteactual,Feuchteactualdes))
        self.Feuchtereal=self.Feuchteactual if self.actual==True else self.Feuchtereal
        self.Feuchteiso=np.asarray([np.average(self.Feuchtereal[idxun][Feuchteun==val]) if len(self.Feuchtereal[idxun][Feuchteun==val])>0 else np.average(Feuchteactualges[Feuchteges==val]) for i,val in enumerate(self.Feuchteiso)]) if self.actual==True else self.Feuchteiso
    
    def FitDiffusionData(self):
        nJumplist=[]
        keylist=[]
        DlistCranc,wreallist,Errorlist,Feuchtelist,tveclist,wtlist,wtreflist,wrealmittellist,mveclist,sratelist,sraterlist,sraterrlist,ARDlist,ARDrlist,m0list,filenamelist=([] for i in range(16))
        backend_=mpl.get_backend()
        mpl.use("Agg")
        Data=self.DATA_FILE
        mol="_"
        Data={0:Data,1:Data} if isinstance(Data,str) else Data
        for ppp in range(len(Data)):
            DCranc,tvecsol,wtrefsol,wtsol,Dsol=([] for i in range(5))
            self.read_excel_file(Data[ppp])
            wreal=self.wreal
            Feuchtereal=self.Feuchtereal
            nJumplist.append(self.nHum)
            mvecsol=[]
            sratesol=[]
            sratersol=[]
            sraterrsol=[]
            ARDsol=[]
            ARDrsol=[]
            for iii in range(nJumplist[ppp]):
                plt.close("all")
                self.indJ=iii
                self.Cranc() 
                DCranc=np.append(DCranc,self.D12) 
                tvecsol.append(self.tlist[iii]/60)
                mvecsol.append(self.mlist[iii]) #if not self.des else mvecsol.append(self.mlist[iii][0]+self.mlist[iii][-1]-self.mlist[iii])
                wtrefsol.append(self.wtref) #if not self.des else wtrefsol.append(self.wtref[0]+self.wtref[-1]-self.wtref)
                wtsol.append(self.wt) #if not self.des else wtsol.append(self.wt[0]+self.wt[-1]-self.wt)
                ttemp=self.tlist[iii]/60
                ARDtemp=np.abs((self.wtref-self.wt))/(self.wtref)
                ARDsol.append(ARDtemp)
                ARDrsol.append(np.average(ARDtemp))
                sratesol.append(self.sratel[iii])
                sratersol.append(np.average(self.sratel[iii][ttemp>60]))
                sraterrsol.append(np.average(self.sratel[iii][0]))
                
            wreallist.append(np.asarray(wreal))
            wrealmittellist.append(np.asarray(self.wrealmittelref))
            Feuchtelist.append(np.asarray(Feuchtereal))
            DlistCranc.append(DCranc)
            tveclist.append(tvecsol)
            mveclist.append(mvecsol)
            sratelist.append(sratesol)
            sraterlist.append(sratersol)
            sraterrlist.append(sraterrsol)
            ARDlist.append(ARDsol)
            ARDrlist.append(ARDrsol)
            wtreflist.append(wtrefsol)
            wtlist.append(wtsol)
            keylist.append(self.dickeys)
            m0list.append(self.mnull*np.ones_like(DCranc))
            filenamelist.append(np.array([Data[ppp]]*len(DCranc)))

        tvecdict=[]
        wrealdict=[]
        wrealmitteldict=[]
        Feuchtedict=[]
        wtrefdict=[]
        wtdict=[]
        DlistCrancdic=[]
        sratedict=[]
        sraterdict=[]
        sraterrdict=[]
        ARDdict=[]
        ARDrdict=[]
        mpl.use(backend_)
        for ppp in range(len(Data)):
            wrealdict.append(OrderedDict(zip( keylist[ppp],wreallist[ppp])))
            tvecdict.append(OrderedDict(zip( keylist[ppp],tveclist[ppp])))
            wrealmitteldict.append(OrderedDict(zip( keylist[ppp],wrealmittellist[ppp])))
            DlistCrancdic.append(OrderedDict(zip( keylist[ppp],DlistCranc[ppp])))
            Feuchtedict.append(OrderedDict(zip( keylist[ppp],Feuchtelist[ppp])))
            wtrefdict.append(OrderedDict(zip( keylist[ppp],wtreflist[ppp])))
            wtdict.append(OrderedDict(zip( keylist[ppp],wtlist[ppp])))
            sratedict.append(OrderedDict(zip( keylist[ppp],sratelist[ppp])))
            sraterdict.append(OrderedDict(zip( keylist[ppp],sraterlist[ppp])))
            sraterrdict.append(OrderedDict(zip( keylist[ppp],sraterrlist[ppp])))
            ARDdict.append(OrderedDict(zip( keylist[ppp],ARDlist[ppp])))
            ARDrdict.append(OrderedDict(zip( keylist[ppp],ARDrlist[ppp])))
        
        wreallist=[]
        Feuchtelist=[]
        tveclist=[]
        wtlist=[]
        wtreflist=[]
        wrealmittellist=[]
        ARDlist=[]
        ARDrlist=[]
        sratelist=[]
        sraterlist=[]
        sraterrlist=[]
        for ppp in range(len(Data)):
            orderdkeys=keylist[ppp]
            wreallist.append(np.asarray([wrealdict[ppp][val] for i,val in enumerate(orderdkeys)]))
            tveclist.append(np.asarray([tvecdict[ppp][val] for i,val in enumerate(orderdkeys)]))
            wrealmittellist.append(np.asarray([wrealmitteldict[ppp][val] for i,val in enumerate(orderdkeys)]))
            Feuchtelist.append(np.asarray([Feuchtedict[ppp][val] for i,val in enumerate(orderdkeys)]))
            wtreflist.append(np.asarray([wtrefdict[ppp][val] for i,val in enumerate(orderdkeys)]))
            wtlist.append(np.asarray([wtdict[ppp][val] for i,val in enumerate(orderdkeys)]))
            sratelist.append(np.asarray([sratedict[ppp][val] for i,val in enumerate(orderdkeys)]))
            sraterlist.append(np.asarray([sraterdict[ppp][val] for i,val in enumerate(orderdkeys)]))
            sraterrlist.append(np.asarray([sraterrdict[ppp][val] for i,val in enumerate(orderdkeys)]))
            ARDlist.append(np.asarray([ARDdict[ppp][val] for i,val in enumerate(orderdkeys)]))
            ARDrlist.append(np.asarray([ARDrdict[ppp][val] for i,val in enumerate(orderdkeys)]))

        Crancave,Crancstd=averagedict(DlistCrancdic)
        wave, wstd=averagedict(wrealdict)
        ARDrave, ARDrstd=averagedict(ARDrdict)
        sraterave, sraterstd=averagedict(sraterdict)
        sraterrave, sraterrstd=averagedict(sraterrdict)
        wmittelave, wmittelstd=averagedict(wrealmitteldict)
        Feuchteave, Feuchtestd=averagedict(Feuchtedict)
        wtave, wtstd=averagelisttimeseries(tveclist,wtlist)
        wtrefave, wtrefstd=averagelisttimeseries(tveclist,wtreflist)
        srateave, sratestd=averagelisttimeseries(tveclist,sratelist)
        ARDave, ARDstd=averagelisttimeseries(tveclist,ARDlist)
        # mvecave, mcvecstd=averagelisttimeseries(tveclist,mveclist)

        tarr=np.array(tveclist)
        nt=tarr.shape[2]
        for i in range(nt):
            tarr[:,:,i]=tarr[:,:,i]-tarr[:,:,0]
        tendn=tarr[:,:,-1]
        amin=np.argmin(tendn,axis=0)
        tveca=[tarr[val,i,:] for i,val in enumerate(amin)]
        nHum=int(DlistCranc[0].shape[0])
        listofdicts={"w2/- ":wreallist,
        "w2ave/- ":[wave],
        "w2std/- ":[wstd],
        "RH/- ":Feuchtelist,
        "RHave/- ":[Feuchteave],
        "RHstd/- ":[Feuchtestd],
        "w2bar/- ":wrealmittellist,
        "w2barave/- ":[wmittelave],
        "w2barstd/- ":[wmittelstd],
        "(DCranc/L^2)/s^-1":DlistCranc,
        "(DCrancave/L^2)/s^-1 ":[Crancave],
        "(DCrancstd/L^2)/s^-1 ":[Crancstd],
        "m0/kg":m0list,
        "ARD/-":ARDrlist,
        "ARDave/-":[ARDrave],
        "ARDstd/-":[ARDrstd],
        "srate/-":sraterlist,
        "srateave/-":[sraterave],
        "sratestd/-":[sraterstd],
        "sraterr/-":sraterrlist,
        "sraterrave/-":[sraterrave],
        "sraterrstd/-":[sraterrstd],
        "filenames":filenamelist}

        listofdicts2={"texpave/min":[tveca],
        "w2expave/-":[wtrefave],
        "w2expstd/-":[wtrefstd],
        "tave/min":[tveca],
        "w2ave/-":[wtave],
        "texp/min":tveclist,
        "w2exp/-":wtreflist,
        "t/min":tveclist,
        "w2/-":wtlist}
        # Plotlist(tveca,wtrefave,tveca,wtave,wtrefstd,wtrefstd,Origin=Film1.Origin,filename="Kinetics_ave_"+mol,xla="t/min",yla="w2/-")
        #WriteToExcel(listofdicts,["Tabelle1"],mol=mol)
        df1=pd.json_normalize(listofdicts, sep='_')
        df2=unnesting(df1,df1.columns.tolist(),axis=0)
        df3=df2.explode(df2.columns.tolist())

        df12=pd.json_normalize(listofdicts2, sep='_')
        df22=unnesting(df12,df12.columns.tolist(),axis=0)
        df32=df22.explode(df22.columns.tolist())
        df42=df32.explode(df32.columns.tolist())

        
        for val in Data:
           with xw.App(visible=False) as app:
                wb = xw.Book(val)
                names=[val.name for val in wb.sheets]
                if "Extracted_Isotherm" not in names: wb.sheets.add("Extracted_Isotherm",after=len(wb.sheets))
                wb.sheets["Extracted_Isotherm"].clear_contents()
                wb.sheets["Extracted_Isotherm"]["A1"].value = df3
                if "Extracted_Kinetics" not in names: wb.sheets.add("Extracted_Kinetics",after=len(wb.sheets))
                wb.sheets["Extracted_Kinetics"].clear_contents()
                wb.sheets["Extracted_Kinetics"]["A1"].value = df42
                wb.save(val)
                wb.close()

    def Cranc(self):
        X=self.wlist[self.indJ]/(1-self.wlist[self.indJ])
        mperc=(X-X[0])/(X[-1]-X[0])
        @njit("f8[::1](f8[::1],f8)",cache=True)
        def Cranc1(x, D):
            n=20
            Forier=np.zeros_like(x)
            LCranc=2.
            for cc in range(n):
                Forier=Forier+1/(2*cc+1)**2*np.exp(-((2*cc+1)**2*np.pi**2*D*x)/(LCranc**2))
                mperc=1-8/np.pi**2*Forier
            return mperc
        popt, pcov = curve_fit(Cranc1, self.tlist[self.indJ], mperc,p0=1E-13)
        mopt=Cranc1(self.tlist[self.indJ],popt[0])
        Xopt=mopt*(X[-1]-X[0])+X[0]
        wopt=Xopt/(1+Xopt)
        self.wt=wopt
        self.wtref=self.wlist[self.indJ]
        self.D12=popt


    @staticmethod
    def interparc(x,y,N):
        """ function that interpolates between a given data series to provide datapoints that are equidistant along the arc of the data.
        Hence the name interp(olate)arc. This is quite handy for kinetic data as the most change in concentration is at earlier times and the least
        change is observed at later times. As a result, usually more data points are in the later stages where nothing happens, since measurements are
        usually performed at equidistant time point """
        xmax=np.max(x)/3
        ymax=np.max(y)
        xmin=np.min(x)/3
        ymin=np.min(y)
        x=(x-xmin)/(xmax-xmin)
        y=(y-ymin)/(ymax-ymin)
        data=np.vstack((x,y)).T
        xd = np.diff(x)
        yd = np.diff(y)
        dist = np.sqrt(xd**2+yd**2)
        u = np.cumsum(dist)
        u = np.hstack([[0],u])

        t = np.linspace(0,u[-1],N)
        xn = InterpolatedUnivariateSpline(u, x,k=1)(t)*(xmax-xmin)+xmin
        yn = InterpolatedUnivariateSpline(u, y,k=1)(t)*(ymax-ymin)+ymin
        return xn,yn

class Checkbar(tk.Frame):
   def __init__(self, parent=None, picks=[],commands= lambda : None, side="top", anchor="w"):
      tk.Frame.__init__(self, parent)
      self.vars = []
      for pick,command in zip(picks,commands):
         var = tk.IntVar()
         chk = tk.Checkbutton(self, text=pick, variable=var,command=command)
         chk.pack(side=side, anchor=anchor, expand="yes")
         self.vars.append(var)
   def state(self):
      return map((lambda var: var.get()), self.vars)


class window3:
    def __init__(self,Film1):
        self.Film=Film1
        self.Film.indF=0
        self.main()

    def enter(self):
        self.Film.indJ=int(self.Box0.get())
    def indexjump(self,indJ):
        indJ=int(self.spinbox.get())-1 if self.spinbox.get()!='' else 0
        self.Film.indJ=indJ
        #self.label3.config(text="Current Jump is "+str(indJ))
        self.PlotRHSteps()
        plt.close('all')
    def indexFile(self,indF):
        indF=int(self.spinbox2.get())-1 if self.spinbox.get()!='' else 0
        self.Film.indF=indF
        self.Film.read_excel_file(self.Film.DATA_FILE[self.Film.indF])
        #self.label3.config(text="Current Jump is "+str(indJ))
        t=self.Film.t_exp
        m=self.Film.w_exp 
        # [val.remove() for i,val in enumerate(self.line4)]
        self.line4.set_data(t,m)
        tbefore=self.Film.t0list[self.Film.indJ][0]
        self.line5.set_data((self.Film.tlist[self.Film.indJ]+tbefore)/60,self.Film.wlist[self.Film.indJ])
        self.PlotRHSteps()
        plt.close('all')
    def PlotRHSteps(self):
        matplotlib.use("TkAgg")
        x=np.arange(0,self.Film.nHum+1,1)
        y=np.hstack((0,self.Film.Feuchtereal)) 
        self.line1.set_data(x,y)
        self.line2.set_data(x[self.Film.indJ:self.Film.indJ+2],y[self.Film.indJ:self.Film.indJ+2])
        [val.remove() for i,val in enumerate(self.line3)]
        self.line3=[]
        for index in range(len(x)-1):
            self.line3.append(self.ax.text(x[index]+0.1, y[index]-0.02, x[index+1], size=7))
        tbefore=self.Film.t0list[self.Film.indJ][0]
        self.line5.set_data((self.Film.tlist[self.Film.indJ]+tbefore)/60,self.Film.wlist[self.Film.indJ])
        self.canvas.draw()
        self.canvas2.draw()
        plt.close('Figure 1')

    def Desorption(self):
        self.Film.des=True if self.Film.des==False else False
        self.Film.read_excel_file(self.Film.DATA_FILE[self.Film.indF])
        self.spinbox.destroy()
        self.spinbox = tk.Spinbox(self.frame, from_ = 1, to = self.Film.nHum,increment = 1,command=partial(self.indexjump, 1))
        self.spinbox.pack(side="top", fill="x")
        #self.label4.config(text="Current Mode: \tDesorption") if self.Film.des==True else self.label4.config(text="Mode is Absorption")

        self.PlotRHSteps()

    def main(self):
     """The main function initializes the GUI-window."""
     self.root3 = tk.Tk()
     self.root3.attributes("-topmost", True)
     #root.geometry("350x150+450+280")  # "width x height + position right + position down"
     self.root3.grid_columnconfigure((0,1), weight=1)


     #Btn3 = tk.Button(self.root3, text="Diffusionsmodell mit mechanischer Zustandsgleichung", padx=0, pady=2, font='Helvetica 7 bold', command=self.Diffusion1DMEOS)
     self.frame=tk.LabelFrame(self.root3,text="Sorption Steps",font='Helvetica 7 bold')

     self.frame.pack(side="top", fill="x")

     self.BtnCh=tk.Button(self.frame, text="Desorption", padx=0, pady=2, font='Helvetica 7 bold', command=self.Desorption)
     self.BtnCh.pack(side="top", fill="x")

     self.spinbox = tk.Spinbox(self.frame, from_ = 1, to = self.Film.nHum,increment = 1,command=partial(self.indexjump, 1))
     self.spinbox.pack(side="top", fill="x")

     self.spinbox2 = tk.Spinbox(self.frame, from_ = 1, to = len(self.Film.DATA_FILE),increment = 1,command=partial(self.indexFile, 1))
     self.spinbox2.pack(side="top", fill="x")

     self.fig,self.ax=plt.subplots(figsize=(4, 2))
     self.fig2,self.ax2=plt.subplots(figsize=(4, 2))
     self.ax.set_ylabel("RH/-")
     self.ax.set_xlabel("Number/-")
     self.ax2.set_ylabel("w/-")
     self.ax2.set_xlabel("t/min")
     self.ax.set_ylim([0,1])
     self.ax.set_yticks([0,0.2,0.4,0.6,0.8,1])
     self.ax.grid()
     self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame)
     self.canvas2 = FigureCanvasTkAgg(self.fig2, master=self.frame)
     
     
     x=np.arange(0,self.Film.nHum+1,1)
     y=np.hstack((0,self.Film.Feuchtereal))
     t=self.Film.t_exp
     m=self.Film.w_exp   
     self.line1, =self.ax.step(x,y,'b',where='pre')
     self.line2, =self.ax.step(x[self.Film.indJ:self.Film.indJ+2],y[self.Film.indJ:self.Film.indJ+2],'r',where='pre')
     self.line3=[]
     for index in range(len(x)-1):
        self.line3.append(self.ax.text(x[index]+0.1, y[index]-0.02, x[index+1], size=7))
     #for index in range(len(t)-1):
     self.line4,=self.ax2.plot(t/60,m,'b')
     tbefore=self.Film.t0list[self.Film.indJ][0]
     self.line5, =self.ax2.plot((self.Film.tlist[self.Film.indJ]+tbefore)/60,self.Film.wlist[self.Film.indJ],'r')
     #self.line2X, =self.ax.step(self.tlist[self.Film.indJ:self.Film.indJ+2],y[self.Film.indJ:self.Film.indJ+2],'r',where='pre')
     self.canvas.get_tk_widget().pack(side="left")
     self.canvas2.get_tk_widget().pack(side="right")
     self.root3.mainloop()



if __name__=="__main__":
    filenames=askopenfilenames()
    Measurement(filenames,tend=-1,nt=15,arc=False)
