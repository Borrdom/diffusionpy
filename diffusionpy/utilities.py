import numpy as np
from math import comb
import itertools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
import mpltern
from matplotlib.pyplot import Axes
import matplotlib.projections as proj
from mpltern.ternary import TernaryAxes

def DasDennis(p,dim):
    """create a equidistance n dimensional spacing which satisfies the mass balance constraint
    
    Examples:
        >>> p=30
        >>> dim=3
        >>> spacevec=DasDennis(p, dim)
        >>> pd.DataFrame(spacevec.T).to_excel("test.xlsx")
    """
    co=(p+dim-2)/p
    nco=int(round(co*p,0))
    covec=np.linspace(0,co,nco+1)
    complexity=comb(dim+p-1,p)
    combinations=np.asarray(list(itertools.combinations(covec, dim-1))).T
    for i in range(dim-1):
        for j in range(complexity):
            combinations[i,j]=combinations[i,j]-i/p
    Nlinspaces=complexity//p
    spacevec=np.zeros((dim,complexity))
    spacevec[0,:]=combinations[0,:]
    for i in range(1,dim-1):
        spacevec[i,:]=combinations[i,:]-combinations[i-1,:]
    spacevec[dim-1,:]=1-combinations[dim-2,:]    
    spacevec[spacevec<=1E-8]=0
    return spacevec


class ternary(TernaryAxes):
    name = 'tern'
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.set_position([0.18,1-0.1416-0.6603,0.6803,0.6803])
        self.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False, bottom=True, top=False, left=True, right=True, direction="in",length=4, width=0.5,labelrotation='horizontal')
        self.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False, bottom=True, top=False, left=True, right=True, direction="in",length=2, width=0.5,labelrotation='horizontal',which="minor")
        self.tick_params(axis="t",pad=10)
        self.tick_params(axis="l",pad=10)
        self.tick_params(axis="r",pad=10)
        self.taxis.set_label_rotation_mode( 'horizontal')
        self.laxis.set_label_rotation_mode( 'horizontal')
        self.raxis.set_label_rotation_mode( 'horizontal')

    def set_labels(self,label="mass fractions / -",title="T = 298.15 K \np = 1 bar",xlabel='solvent',ylabel='polymer',zlabel="API"):
        self.text(s=label,x=470, y=80)
        self.text(s=title ,x=50, y=700)
        self.set_tlabel(xlabel)
        self.set_llabel(ylabel)
        self.set_rlabel(zlabel)
        self.taxis.set_minor_locator(AutoMinorLocator(2))
        self.laxis.set_minor_locator(AutoMinorLocator(2))
        self.raxis.set_minor_locator(AutoMinorLocator(2))

    def filled_line(self,x,y,z,Formatstring,legend):
        p=self.plot(x, y, z,Formatstring,linewidth=1,label=legend+"_filled")
        color = p[0].get_color()
        self.fill(x, y, z, alpha=0.2,color=color,label=legend+"_filled")

    def conodes(self,RBx,RBy,RBz,LBx,LBy,LBz,Formatstring,legend):
        self.plot(RBx,RBy,RBz,Formatstring,linewidth=1,label=legend)
        self.plot(LBx,LBy,LBz,Formatstring,linewidth=1,label=legend)
        
        for i,(rt,rl,rr,lt,ll,lr) in enumerate(zip(RBx,RBy,RBz,LBx,LBy,LBz)):
                self.plot([rt,lt],[rl,ll],[rr,lr],"k-",linewidth=0.5,label=f"Konode {i}")

def get_line_from_tern(cs):
    p = cs.collections[0].get_paths()[0]
    v = p.vertices
    x = v[:,0]
    y = v[:,1]
    a=y
    b=0.5*(1-a)-x*np.sqrt(3)/2
    c=1-a-b
    return a,b,c

proj.register_projection(ternary)
