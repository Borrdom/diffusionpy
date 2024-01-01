import numpy as np
import matplotlib.pyplot as plt
#import OriginVorlagePlot#added 03.03.2020 for origin
from PyPlorigin.PyOrigin import Contour,MusterDictOrigin,Plot
def WriteHistorycsv(tspan,wt,wtref,wtstd=None,wtrefstd=None,tsim=None,msim=None):
    import pandas as pd
    wtstd=np.zeros_like(wt) if wtstd is None else wtstd
    wtrefstd=np.zeros_like(wtref) if wtrefstd is None else wtrefstd
    tsim=np.zeros_like(tsim) if tsim is None else tsim
    msim=np.zeros_like(msim) if msim is None else msim
    pd.DataFrame(data=np.asarray([tspan/60,wt*100,wtref*100,wtrefstd*100,wtstd*100,tsim/60,msim*100]).T,columns=["t[min]","wSim[%]","wExp[%]","wSimStd[%]","wExpStd[%]","tsim[min]","msim[%]"]).to_csv("history.csv",sep=";")
def Plotlist(x,y,xs,ys,xe=0,ye=0,Origin=True,filename="default",xla="?",yla="?"):
    fig,ax=plt.subplots()
    [ax.plot(xs[i],ys[i],'k-',linewidth=2) for i,val in enumerate(x)]
    [ax.plot(x[i],y[i],'o',linewidth=0.5,markersize=10,markeredgecolor=(0,0,0)) for i,val in enumerate(xs)]
    xData,yData,xSim,ySim,xErr,yErr=MusterDictOrigin()
    xe=[ np.zeros_like(val) for i,val in enumerate(x)] if isinstance(xe,int) else xe
    ye=[ np.zeros_like(val) for i,val in enumerate(x)] if isinstance(ye,int) else ye
    ax.set_xlabel(xla)
    ax.set_ylabel(yla)
    xData["values"],yData["values"],    xSim["values"],ySim["values"],xErr["values"],yErr["values"]=x,y,xs,ys,xe,ye
    Plot(xData,yData,xSim,ySim,xErr,yErr,"","",filename) if Origin==True else None
    return fig,ax
def Plotcontour(xx,yy,zz,Origin=True,filename="default",xla="?",yla="?",zla="?"):
    X=xx.flatten()
    Y=yy.flatten()
    Z=zz.flatten()
    fig,ax=plt.subplots()
    cont=ax.contourf(xx,yy,zz)
    ax.set_xlabel(xla)
    ax.set_ylabel(yla)
    cbar=plt.colorbar(cont)
    cbar.set_label(zla, rotation=270)
    Contour([np.asarray([X,Y,Z])],"Contour") if Origin else None

if __name__=="__main__":
    a = np.array(range(20))
    b = a * 2
    c = a * 3
    d = a *4
    xl1=[a,a,a,a]
    yl1=[a,b,c,d]
    xlr=[a,a,a,a]
    ylr=[a,b,c,d]
    xe=[a-a,a-a,a-a,a-a]
    ye=[a-a,a-a,a-a,a-a]
    x=y=np.linspace(-3,3,100)
    xx,yy=np.meshgrid(x,y)
    zz=xx+yy**2
    Plotcontour(xx,yy,zz,xla="x/-",yla="y/-",zla="z/-")

    #Plotlist(xl1,yl1,xlr,ylr,xe,ye)
    plt.show()
