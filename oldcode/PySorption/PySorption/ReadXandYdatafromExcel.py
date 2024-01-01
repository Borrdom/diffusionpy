import pandas as pd
from tkinter.filedialog import askopenfilename
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
import PySorption.interparc
from scipy.ndimage import median_filter
import xlrd
"""nifty little class that lets you read columns of excel sheets. Sheet is search by a string. If no string is
provided instead a askopenfilename dialog pops up. The columns is searched by the getrow method. It requires
a string of the header of the sought-after column of the excel-sheet"""
class sheet:
    def __init__(self,filename=askopenfilename,sheet_name=0):
        self.filename=filename if not callable(filename) else filename()
        self.Data=pd.read_excel(self.filename,sheet_name=sheet_name,header=None) if ".xlsx" in self.filename else pd.read_csv(self.filename,header=None)

    def getrows(self,*args):
        return tuple([self.getrow(val) for i,val in enumerate(args)])

    def getrow(self,strx):
        boolalpha=any([val.isalpha() for val in strx])
        if boolalpha:
            Datanumberheader=[any(self.Data[val].astype(str).isin([strx]).fillna(False)) for i,val in enumerate(self.Data)]
            ncol=self.Data[self.Data==strx].dropna(how='all').index
            #ncol=self.Data[self.Data.columns[np.asarray(Datanumberheader)]].astype(str).isin([strx]).fillna(False)
            Datanumberheader2=[self.Data[val][ncol[0]+1:].dropna() for i,val in enumerate(self.Data) if Datanumberheader[i]]
            defindex=0
        else:
            Datanumberheader2=[]
            defindex=int(strx)-1
        Datanumberheader2=[self.Data[defindex][1:].dropna()] if Datanumberheader2==[] else Datanumberheader2
        return Datanumberheader2[0].values.astype(float)



def FasterXlrdRead(filename,start,*args):
    book = xlrd.open_workbook(filename, encoding_override = "utf-8")
    sheet = book.sheet_by_index(0)

    col1=sheet.col(0)
    col1arr=np.asarray([val.value for i,val in enumerate(col1)])
    nstart=np.where(col1arr==start)[0][0]+1 if np.where(col1arr==start)[0].size else 0
    nrow=sheet.nrows
    ncol=sheet.ncols
    nrow2=nrow-nstart
    #if any
    class MyException(Exception):
        pass


    header=np.asarray([sheet.cell(nstart-1, i).value for i in range(ncol)])
    soughtrows=np.asarray([np.where(val==header)[0][0] if val in header else -1 for i,val in enumerate(args) ])
    if soughtrows[-1]==-1:
        raise MyException(args[-1]+" not found")
    soughdata=[np.asarray([sheet.cell(i+nstart, valj).value if valj>=0 else 0 for i in range(nrow2)]) for j,valj in enumerate(soughtrows)] # [s]
    return tuple(soughdata)
def CombineKinetikData(strt,strx,filenames):



    def averagelist(lis):
        ave=sum(lis)/len(lis)
        std=(sum([(ave-val)**2 for val in lis])/(len(lis)-1))**0.5
        return ave,std


    t=[]
    x=[]
    tend=[]
    for i,val in enumerate(filenames):
        sheet1=sheet(val)
        t.append(sheet1.getrow(strt))
        x.append(sheet1.getrow(strx))
        tend.append(t[i][-1])
    amin,amax=np.argmin(tend),np.argmax(tend)
    #tvec,xvec=interparc.interparc(t[amin],median_filter(x[amin],100,mode="nearest"),200)
    tvec,xvec=interparc.interparc(t[amin],x[amin],200)
    #tvec,xvec=t[amin],x[amin]
    #tvec=tvec[tvec>0]
    tvec=tvec[tvec>=0]
    #tvec=np.linspace(0,tend[amin]**(1/2))**2
    xvec=[InterpolatedUnivariateSpline(valt,valx)(tvec) for (valt,valx) in zip(t,x)]
    fig,ax=plt.subplots()

    [ax.plot(valt,valx,'x') for (valt,valx) in zip(t,x)]
    [ax.plot(tvec,valx,'o') for valx in xvec]

    xave,xstd=averagelist(xvec)
    ax.errorbar(x=tvec,y=xave,yerr=xstd,fmt='ko')
    return [tvec]*len(xvec),xvec,xave,xstd

def CombineData(strt,strx,filenames):

    t=[]
    x=[]
    tend=[]
    for i,val in enumerate(filenames):
        sheet1=sheet(val)
        t.append(sheet1.getrow(strt))
        x.append(sheet1.getrow(strx))
    Dataheader=tuple([sheet1.Data[val][0] for i,val in enumerate(sheet1.Data)])
    #allrows=sheet1.getrows(*Dataheader)
    return t,x,x[0],np.zeros_like(x[0])
