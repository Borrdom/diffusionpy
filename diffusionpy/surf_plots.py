
import numpy as np
import matplotlib.pyplot as plt
def circular(t,zvec,wtz,Lt=None,instances=6,comp=0,cmap="Blues"):
    L=zvec[-1]
    expansion=Lt[:,None]/L if Lt is not None else np.ones_like(t)
    phi=np.linspace(0,2*np.pi,41)
    Rad,Phi=np.meshgrid(zvec*1E6,phi)
    fig, axes = plt.subplots(2,instances//2, constrained_layout=True,subplot_kw={'projection': 'polar'})
    # axes=[]
    axes=axes.flatten()
    pls=[]
    delt=len(t)//instances
    for i in range(instances):
        axes[i].grid(False)
        # axes.append(fig.add_subplot(2,instances//2,i+1, polar=True))
        pls.append(axes[i].contourf(Phi,Rad*expansion[delt*i],np.meshgrid(wtz[delt*i,comp,:],phi)[0],cmap=cmap,vmin=0, vmax=1))
        axes[i].grid(False)
        axes[i].set_xticklabels([])
        axes[i].set_yticklabels([])
        axes[i].set_ylim(0, np.max(zvec*1E6*expansion))
        axes[i].spines['polar'].set_visible(False)
        axes[i].set_title(f'{t[delt*i]/60:.2f}'+" min", va='bottom')

    axes=np.asarray(axes)
      
    # ax30 = fig10.add_subplot(3,instances//2,instances+1)
    # sm = plt.cm.ScalarMappable(cmap="Blues", norm=plt.Normalize(vmin=0, vmax=1))
    # fig10.colorbar(sm,ax30,orientation='horizontal')
    fig.colorbar(pls[0], ax=axes.ravel().tolist(),orientation="horizontal")
    fig.subplots_adjust(hspace=0,wspace=0)  