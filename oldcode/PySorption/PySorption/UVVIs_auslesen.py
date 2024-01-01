import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np

from datetime import datetime
filename="Stammlösung PVA_0,00024%.csv"
filename="Auflösung_PVPVA_Krummnow.csv"
def Sotax(filename):
    foldername="DissolutionData"
    substancename="pvp"
    solventname="water"
    path=os.path.join(os.getcwd(),foldername,substancename,solventname,filename)

    data=pd.read_csv(path,delimiter=";")


    Wellenlänge=data.columns[5:-1]
    Signal=np.asarray(data[Wellenlänge])
    lambdas=np.asarray(Wellenlänge,dtype=float)

    form1='%d.%m.%Y %H:%M:%S'
    form2='%d.%m.%Y %H:%M'
    def gettime(format):
        Zeit0=datetime.strptime(data["Datum/Zeit"][0],format) 
        Zeit=np.asarray([(datetime.strptime(val,format)-Zeit0).total_seconds() for i,val in enumerate(np.asarray(data["Datum/Zeit"]))])
        return Zeit
    try:
        Zeit=gettime(format=form1)    
    except:
        Zeit=gettime(format=form2)  
    print(Zeit)

    fig1,ax1=plt.subplots()
    fig2,ax2=plt.subplots()
    [ax1.plot(lambdas,Signal[i,:]) for i,val in enumerate(Signal[:,0])]
    ax1.set_xlabel("Wellenlänge/nm")
    ax1.set_ylabel("Extinktion/-")
    [ax2.plot(Zeit,Signal[:,i]) for i,val in enumerate(Signal[0,:])]
    ax2.set_xlabel("Zeit/s")
    ax2.set_ylabel("Extinktion/-")
    mw=703
    m0=data["Notiz"][0]
    print(m0/(mw+m0))

    plt.show()      
    return Zeit,Signal,m0

if __name__=="__main__":
    Zeit,Signal,m0=Sotax(filename)


#Wellenlänge=data[0,:]
#header=data[:,0]
