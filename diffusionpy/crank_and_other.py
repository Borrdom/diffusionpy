import numpy as np
def crank(t,kf,mfinfty):
    """cranks equation for diffusion"""
    ninfty=100
    mtdiff=mfinfty*(1-sum([8/np.pi**2*1/(2*n+1)**2*np.exp(-(2*n+1)**2*kf*t) for n in range(ninfty)]))
    return mtdiff
def Relaxation(t,kr,nr,mrinfty):
    """"used by the Hophenberg modification """
    return sum([(1-np.exp(-kr[n]*t))*mrinfty[n] for n in range(nr)])
def BHModel(t,kf,kr,nr,mfinfty,mrinfty):
    """"Berens and Hopfenberg modification of crank equation for relaxation"""
    return crank(t,kf,mfinfty)+Relaxation(t,kr,nr,mrinfty)
def BHX(t,kf,kr,nr,X0,Xinfty,mfinfty,mrinfty):
    """"used by the Hophenberg modification """
    M_M=BHModel(t,kf,kr,nr,mfinfty,mrinfty)
    return (1-M_M)*X0+M_M*Xinfty