from diffusionpy import vpure,lngi
import numpy as np
import matplotlib.pyplot as plt
def testit(Mw):
    T=298.15
    p=1E5
    #Water
    #XAi0 =0.16 und 0.6
    par={"mi":np.asarray([1.20469,0.0489*Mw]),
    "si":np.asarray([2.797059952,2.889]),
    "ui":np.asarray([353.95,298.047]),
    "eAi":np.asarray([2425.67,1602.3]),
    "kAi":np.asarray([0.04509,0.02]),
    "NAi":np.asarray([1.,931/7335*0.0489*Mw]),
    "kij":np.asarray([[0,0],[0,0]])}
    xi=np.asarray([1.-1E-8,1E-8])
    vpures=vpure(p,T,**par)

    par["vpure"]=vpures
    return lngi(T,xi,**par)

def initialize():
    T=298.15
    p=1E5
    npoint=2
    #Water
    #XAi0 =0.16 und 0.6

    par={"mi":np.asarray([1.20469,2.38269789]),
    "si":np.asarray([2.797059952,3.1771]),
    "ui":np.asarray([353.95,198.24]),
    "eAi":np.asarray([2425.67,2653.4]),
    "kAi":np.asarray([0.04509,0.032384]),
    "NAi":np.asarray([1.,1.])}
    x1=np.linspace(0,1,npoint)
    x2=1-x1
    xi=np.vstack((x1,x2))
    vpures=vpure(p,T,**par)
    par["vpure"]=vpures
    lngiammai=np.asarray([lngi(T,xi[:,i],**par).flatten() for i,val in enumerate(xi[0,:])])
    return lngiammai

lngiammai=initialize()
Mw=np.linspace(149990,150000,11)

lngis=np.asarray([testit(val) for val in Mw])
lng1=lngis[:,0]
lng2=lngis[:,1]
plt.plot(Mw,lng1,"ko")
plt.show()
import pandas as pd
pd.DataFrame((Mw,lng1)).T.to_clipboard(excel=True, sep=None, index=False, header=None)