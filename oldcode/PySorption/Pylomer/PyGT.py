import numpy as np

def Tggt(wi,rhoi,Tgi,q):
    nc=rhoi.shape[0]
    #ki=rhoi[idx]*Tgi[idx]/(rhoi*Tgi)
    qmat=np.zeros((nc,nc))
    qmat[np.triu_indices(nc, k=1)]=q
    #qmat=qmat[0]
    Excess=np.asarray([np.sum(np.prod(np.meshgrid(wi[i,:],wi[i,:]),axis=0)*qmat) for i,val in enumerate(wi[:,0])])
    Ideal=np.sum(wi*1/rhoi,1)/np.sum(wi*1/rhoi/Tgi,1)
    #Ideal=np.sum(wi*Tgi*ki,1)/np.sum(wi*ki,1)
    return Ideal+Excess
    #Triangular numbers for binary interaction parameters
if __name__=="__main__":
    import matplotlib.pyplot as plt
    w2=np.linspace(0,1,100)
    w1f=0
    w1=(1-w2)*w1f
    w3=(1-w2)*(1-w1f)
    #w3=np.zeros_like(w1)
    #rhoi=np.asarray([1320,997,1216]) #ind w pvp
    rhoi=np.asarray([1320,997,1190]) #ind w pvac
    rhoi=np.asarray([1150,997,1216]) #rit w pvp
    #rhoi=np.asarray([1150,997,1250])
    Tgi=np.asarray([317.6,136,441.51])#ind w pv

    #Tgi=np.asarray([317.6,136,383.9])
    #Tgi=np.asarray([317.6,136,383.9])
    #Tgi=np.asarray([317.6,136,383.9])
    writ=np.asarray([0,2.069,3.1149,3.8391,4.8046])/100
    writ=writ/(1+writ)
    Trit=np.asarray([50.718,32.456,27.911,23.923,19.856])+273.15



    #q=np.asarray([-600,0,0])
    #w1=np.linspace(0,1,100)
    #w2f=np.linspace(0,1,100)
    #w2=w2f*(1-w1)
    #w3=(1-w2f)*(1-w1)
    wi=np.asarray([w1,w2,w3])

    rhoi[0]=(Tgi[1]*rhoi[1])/(0.11*Tgi[0]) #ind
    #rhoi[2]=(Tgi[1]*rhoi[1])/(0.33*Tgi[2]) #pvp
    #rhoi[0]=(Tgi[1]*rhoi[1])/(0.244*Tgi[0]) #rit
    #Tgw1=Tggt(wi.T,rhoi,Tgi,q=q)
    Tgw2=Tggt(wi.T,rhoi,Tgi,q=0)
    import pandas as pd
    pd.DataFrame(w2,Tgw2).to_excel("fil.xlsx")
    plt.plot(w2,Tgw2)

    #y=(Tgi[0]-Tgw2)/(Tgw2-Tgi[1])
    #x=w2/(1-w2)
    #plt.plot(x,y,'kx')
    plt.plot([0,1],[298.15,298.15])
    plt.show()
    #plt.plot(writ,1/Trit,'kx')
