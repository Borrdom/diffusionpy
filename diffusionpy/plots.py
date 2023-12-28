
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
import mpltern
from scipy.interpolate import interp1d


def circular(t,zvec,wtz,Lt=None,instances=6,comp=0,cmap="Blues",vmin=None,vmax=None,label=None,tinterp=None):
    matplotlib.rcParams['mathtext.fontset'] = 'custom'
    matplotlib.rcParams['mathtext.rm'] = 'Calibri'
    matplotlib.rcParams['mathtext.it'] = 'Calibri'
    matplotlib.rcParams['mathtext.bf'] = 'Calibri'
    matplotlib.rcParams['xtick.major.pad']='5'
    matplotlib.rcParams['ytick.major.pad']='5'
    matplotlib.rcParams['axes.linewidth'] = 0.5
    # matplotlib.rcParams["toolbar"] = "toolmanager"
    # plt.rcParams['axes.autolimit_mode'] = 'round_numbers'
    
    font = {'weight' : 'normal',
            'size'   : 16,
            'family' : "calibri"}
    plt.rc('font', **font)
    plt.rc('axes', titlesize=font["size"]) 
    L=zvec[-1]
    expansion=Lt[:,None]/L if Lt is not None else np.ones_like(t)[:,None]
    phi=np.linspace(0,2*np.pi,41)
    Rad,Phi=np.meshgrid(zvec*1E6,phi)

    fig, axes = plt.subplots(2,instances//2, constrained_layout=True,subplot_kw={'projection': 'polar'},figsize=(7,5),dpi=250)
    # axes=[]
    axes=axes.flatten()
    pls=[]
    delt=len(t)//instances
    if vmin is None: vmin=np.min(wtz) 
    if vmax is None: vmax=np.max(wtz)
    wtz_fun=interp1d(t,wtz,axis=0)
    if tinterp is None: tinterp=np.linspace(t[0],t[-1],instances) 
    for i in range(instances):
        axes[i].grid(False)
        # axes.append(fig.add_subplot(2,instances//2,i+1, polar=True))
        pls.append(axes[i].contourf(Phi,Rad*expansion[delt*i],np.meshgrid(wtz_fun(tinterp[i])[comp,:],phi)[0],cmap=cmap,vmin=vmin,vmax=vmax))
        axes[i].grid(False)
        axes[i].set_xticklabels([])
        axes[i].set_yticklabels([])
        axes[i].set_ylim(0, np.max(zvec*1E6*expansion))
        axes[i].spines['polar'].set_visible(False)
        axes[i].set_title(f'{tinterp[i]/60:.2f}'+" min", va='bottom')

    axes=np.asarray(axes)
      
    cbar=fig.colorbar(pls[-1], ax=axes.ravel().tolist(),orientation="horizontal")
    if label is None: label='$w_i / -$'
    cbar.set_label(label)
    fig.subplots_adjust(hspace=0,wspace=0)  
    return fig, axes

def basic_colors(Formatstring):
    if "g" in Formatstring: return "#99CC00" #green
    if "c" in Formatstring: return "#99CDE9" #cyan
    if "b" in Formatstring: return "#246FE2" #blue
    if "r" in Formatstring: return "#FF8500" #orange
    if "m" in Formatstring: return "#FFCCCC" #magenta
    if "y" in Formatstring: return "#FFD67E" #yellow
    if "a" in Formatstring: return "#666666" #gray
    if "k" in Formatstring: return "#000000" #black
    return "#000000"
class origin_like:

    def subplots(sharex=False):
        matplotlib.rcParams['mathtext.fontset'] = 'custom'
        matplotlib.rcParams['mathtext.rm'] = 'Calibri'
        matplotlib.rcParams['mathtext.it'] = 'Calibri'
        matplotlib.rcParams['mathtext.bf'] = 'Calibri'
        matplotlib.rcParams['xtick.major.pad']='5'
        matplotlib.rcParams['ytick.major.pad']='5'
        matplotlib.rcParams['axes.linewidth'] = 0.5
        matplotlib.rcParams['axes.axisbelow'] = True
        # matplotlib.rcParams["toolbar"] = "toolmanager"
        # plt.rcParams['axes.autolimit_mode'] = 'round_numbers'
        
        font = {'weight' : 'normal',
                'size'   : 16,
                'family' : "calibri"}
        plt.rc('font', **font)
        fig=plt.figure(figsize=(4 , 4), dpi = 250)
        
        if sharex :
            ax=fig.add_axes([0.2667,1-0.2042-0.5017,0.5908,0.5017/2]) 
            ax2=fig.add_axes([0.2667,1-0.2042-0.5017,0.5908,0.5017/2])
            ax.get_shared_x_axes().join(ax, ax2) 

        else:
            ax=fig.add_axes([0.2667,1-0.2042-0.5017,0.5908,0.5017])
        ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False, bottom=True, top=True, left=True, right=True, direction="in",length=4, width=0.5)
        return fig,ax
    def set_xlabel(ax,xlabel1,xunit1=None):
        ax.set_xlabel(f'$\mathrm{{{xlabel1}}}$ / $\mathrm{{{xunit1}}}$') if xunit1 is not None else ax.set_xlabel(f'$\mathrm{{{xlabel1}}}$') 
    def set_ylabel(ax,ylabel1,yunit1=None):
        ax.set_ylabel(f'$\mathrm{{{ylabel1}}}$ / $\mathrm{{{yunit1}}}$',linespacing=1.5) if yunit1 is not None else ax.set_ylabel(f'$\mathrm{{{ylabel1}}}$',linespacing=1.5)
    def plot(ax,x,y,Formatstring,label=None,order=1,yerr=None,z=None):
        
        if z is not None:
            ax.plot(x,y,z,Formatstring , zorder=order,linewidth = 1.5,label=label,markersize=5, markeredgecolor='k',markeredgewidth=0.5,color=basic_colors(Formatstring)) 
        else:
            if yerr is not None:
                if sum(yerr)==0 : yerr=None
                ax.errorbar(x,y,yerr,None,"ko" , zorder=order,label=label,markersize=0, markeredgecolor='k',markeredgewidth=0.5,capsize=5, elinewidth=0.5,ecolor="k")
            ax.plot(x,y,Formatstring , zorder=order,linewidth = 1.5,label=label,markersize=5, markeredgecolor='k',markeredgewidth=0.5,color=basic_colors(Formatstring))

        
        # ax.plot(x,y,Formatstring , zorder=order,linewidth = 1.5,label=label,markersize=5, markeredgecolor='k',markeredgewidth=0.5)
    def set_ticks(ax,x0=None,x1=None,y0=None,y1=None):
        if x0 is not None: ax.axis([x0, x1, y0, y1]) 
        start, end = ax.get_xlim()
        ax.xaxis.set_ticks(np.linspace(start, end, 5))
        start, end = ax.get_ylim()
        ax.yaxis.set_ticks(np.linspace(start, end, 5))

    def ternary():
        matplotlib.rcParams['mathtext.fontset'] = 'custom'
        matplotlib.rcParams['mathtext.rm'] = 'Calibri'
        matplotlib.rcParams['mathtext.it'] = 'Calibri'
        matplotlib.rcParams['mathtext.bf'] = 'Calibri'
        # matplotlib.rcParams['xtick.major.pad']='5'
        # matplotlib.rcParams['ytick.major.pad']='5'
        matplotlib.rcParams['axes.linewidth'] = 0.5
        # matplotlib.rcParams["toolbar"] = "toolmanager"
        # plt.rcParams['axes.autolimit_mode'] = 'round_numbers'

        font = {'weight' : 'normal',
                'size'   : 16,
                'family' : "calibri"}
        plt.rc('font', **font)

        fig = plt.figure(figsize=(6, 5),dpi=200)
        fig.subplots_adjust( wspace=0,hspace=0)
        ax = fig.add_axes(projection="ternary",rect=[0.18,1-0.1416-0.6603,0.6803,0.6803])
        ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False, bottom=True, top=False, left=True, right=True, direction="in",length=4, width=0.5,labelrotation='horizontal')
        ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False, bottom=True, top=False, left=True, right=True, direction="in",length=2, width=0.5,labelrotation='horizontal',which="minor")
        ax.tick_params(axis="t",pad=10)
        ax.tick_params(axis="l",pad=10)
        ax.tick_params(axis="r",pad=10)
        ax.taxis.set_label_rotation_mode( 'horizontal')
        ax.laxis.set_label_rotation_mode( 'horizontal')
        ax.raxis.set_label_rotation_mode( 'horizontal')
        return fig,ax
    def set_labels(ax,label="mass fractions / -",title="T = 298.15 K \np = 1 bar",xlabel='solvent',ylabel='polymer',zlabel="API"):
        ax.text(s=label,x=470, y=80)
        ax.text(s=title ,x=50, y=700)
        ax.set_tlabel(xlabel)
        ax.set_llabel(ylabel)
        ax.set_rlabel(zlabel)
        ax.taxis.set_minor_locator(AutoMinorLocator(2))
        ax.laxis.set_minor_locator(AutoMinorLocator(2))
        ax.raxis.set_minor_locator(AutoMinorLocator(2))

    def filled_line(ax,x,y,z,Formatstring,legend):
        p=ax.plot(x, y, z,Formatstring,linewidth=1,label=legend+"_filled",color=basic_colors(Formatstring))
        color = p[0].get_color()
        ax.fill(x, y, z, alpha=0.2,color=color,label=legend+"_filled")

    def conodes(ax,RBx,RBy,RBz,LBx,LBy,LBz,Formatstring,legend):
        ax.plot(RBx,RBy,RBz,Formatstring,linewidth=1,label=legend,color=basic_colors(Formatstring))
        ax.plot(LBx,LBy,LBz,Formatstring,linewidth=1,label=legend,color=basic_colors(Formatstring))
        
        for i,(rt,rl,rr,lt,ll,lr) in enumerate(zip(RBx,RBy,RBz,LBx,LBy,LBz)):
                ax.plot([rt,lt],[rl,ll],[rr,lr],"-",linewidth=0.5,label=f"Konode {i}",color=basic_colors(Formatstring))
