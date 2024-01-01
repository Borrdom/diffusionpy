from scipy.interpolate import interp1d
import numpy as np
import pandas as pd



# def averagelist(lis):
#     ave=0
#     for i,val in enumerate(lis):
#         ave+=val/len(lis)
#     var=0
#     for i,val in enumerate(lis):
#         var+=(ave-val)**2/(len(lis)-1)
#     std=var**0.5
#     x, y=ave,std if len(lis)>0 else (lis[0],lis[0]*1E-14)
#     return x,y

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
