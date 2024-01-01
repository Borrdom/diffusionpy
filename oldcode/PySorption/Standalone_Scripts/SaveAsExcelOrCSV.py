import pandas as pd
import os


def SaveAsExcelOrCSV(pathname,headerlist,datalist):
    """To import it write << from SaveAsExcelAndCSV import SaveAsExcelAndCSV >>
    To call it <<SaveAsExcelAndCSV(pathname,headerlist,datalist)>> where pathname is the full path to the file
    headerlist is a list of the rows that will be above each column of data
    datalist is a list of 1D Arrays which will be presented as columns
    headerlist and datalast must be the same length"""
    Fileextension=(pathname).split(".")[1]
    dicoflist=dict(zip(headerlist,datalist))
    dataframeofdict=pd.DataFrame(dicoflist)
    dataframeofdict.to_excel(pathname,sheet_name="Sheet1",index=None) if (Fileextension=="xlsx" or Fileextension=="xls") else None
    dataframeofdict.to_csv(pathname,index=None) if Fileextension=="csv" else None
    

if __name__=="__main__":
    import numpy as np
    x=np.linspace(0,1,10)
    y=np.linspace(1,0,10)
    z=x*y
    headerlist=["x/m","y/1/s","z/m/s"]
    datalist=[x,y,z]
    Foldername="Testnames"
    Filename="Test"
    cwd=os.getcwd()
    pathname=os.path.join(cwd,Foldername,Filename)
    dirname=os.path.join(cwd,Foldername)
    os.mkdir(dirname) if not os.path.exists(dirname) else None
    SaveAsExcelOrCSV(pathname+".xlsx",headerlist,datalist)
    SaveAsExcelOrCSV(pathname+".csv",headerlist,datalist)
    