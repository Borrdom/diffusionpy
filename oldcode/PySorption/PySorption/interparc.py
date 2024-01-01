import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.ndimage import median_filter
from numpy import inf

def interparc(x,y,N):
    #x=median_filter(x,500,mode='nearest')
    #y=median_filter(y,500,mode='nearest')

    # xlog=np.log10(np.abs(x))
    # ylog=np.log10(np.abs(y))
    # xlog=xlog[xlog != -inf]
    # ylog=ylog[ylog != -inf]
    # xmag=np.floor(np.average(xlog))
    # ymag=np.floor(np.average(ylog))
    #xmax=50*60
    #xmax=50*60


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

if __name__=="__main__":
    N=10
    Np=100
    x=np.linspace(0.1,1,Np)
    y=1-np.exp(-10*x)
    xn,yn=interparc(x,y,N)
    f = plt.figure()
    ax = f.add_subplot(111)
    ax.plot(x,y,'o')
    ax.plot(xn,yn,'ro')
    plt.show()



#For multidimensional data with rows x dimension
# diffs = data[1:, :] - data[:-1, :]
# dist = np.linalg.norm(diffs, axis=1)
# u = np.cumsum(dist)
# u = np.hstack([[0], u])
# t = np.linspace(0, u[-1], 10)
# resampled = scipy.interpolate.interpn((u,), pts, t)
