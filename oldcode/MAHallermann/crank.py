# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 19:22:50 2021

@author: Moritz
"""
# A script to plot the Crank equation with the desired diffusion coefficient.

import numpy as np
import matplotlib.pyplot as plt

pi = np.pi

def crank(t,d,l):
    ns=30
    mt_minf = 1-sum([(8/pi**2)*(1/(2*n + 1)**2)*np.exp((-d*((n + 1/2)**2)*(pi**2)*t)/l**2) for n in range(ns)])
    return mt_minf
# All settings to generate a nice plot which is then safed.    
def PlotCrank(time_vec,mt_minf):
    plot_font = {'fontname':'Calibri'}
    plt.rc('xtick',labelsize=18)
    plt.rc('ytick',labelsize=18)
    plt.ioff()
    fig1,ax1 = plt.subplots()
    ax1.plot(time_vec,mt_minf,'g', linewidth=2.5)
    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(1.5)
    plt.ylim(top=1.05,bottom=0) 
    ax1.set_xlabel(r'$\mathregular{t \ / \ h}$',**plot_font,size=18)
    ax1.set_ylabel(r'$\mathregular{\dfrac{m_t}{m_{\infty}} \ / -}$',**plot_font,size=18,loc='top',rotation=0)  
    ax1.tick_params(direction='in', right=True, top=True, labelright=False, labeltop=False, width=1.5, length=7)
    ax1.grid()
    plt.tight_layout()
    plt.show()
    plt.savefig('Crank_Plot_DifCoe_'+str(float(dif_coe_input)))
    #return mt_minf

# The user can either type in the desired value for the diffusion coefficient or it is set to a fix value.
if __name__=="__main__":
    dif_coe_input = input("Please enter the value of the desired diffusion coefficient [m^2/s]: ")
    tend_raw = 216
    tend = tend_raw*3600
    time_vec = np.linspace(0,tend_raw,300)
    l = 5*10**-6
    if dif_coe_input == []:
        mt_minf=crank(tend_raw,2.255*10**-12,l)
    else:
        mt_minf=crank(tend_raw,float(dif_coe_input),l)
    PlotCrank(time_vec,mt_minf)
    
        