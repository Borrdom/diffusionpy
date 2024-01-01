# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 12:46:12 2021

@author: Moritz
"""

import pandas as pd
import numpy as np
import casadi as cs
import scipy as sp
import os
from scipy import constants
import matplotlib.pyplot as plt
import PyCSAFTNil
from PyCSAFTNil import Mixture
import get_pcpar
import crank
from ExtractWeidemannData import ExtractIsotherm,ExtractSorption
from scipy.optimize import curve_fit


# PC-SAFT is used to calculate the equilibrium chemical potemtial as well as a variety of chemical potentials belonging 
# to a specific mass fraction. First, it is validated whether the data of the desired component has already been calculated 
# by checking the existence of the file.
def PCSAFT(temp,polstr,apistr,solstr,Dw,dl):

    excelstr = [polstr,apistr] 
    if Dw == 0:
        excelstr = [polstr,apistr]
    elif dl == 1:
        excelstr = [apistr,solstr]
    else:
        excelstr = [polstr,apistr,solstr] 
    
    currentpath = os.getcwd()
    foldername = "PyCSAFTSheets"
    tablename_1 = "sle_napvp"
    tablename_2 = "GammaScan"# DB: changed to the new function call for the PCSAFT Scan
    
    filename_1 ="_".join([tablename_1]+excelstr)+".xlsx"    
    filename_2 ="_".join([tablename_2]+excelstr)+".xlsx"   
    filepath = os.path.join(currentpath,foldername)
    pathname_1 = os.path.join(currentpath,foldername,filename_1)
    pathname_2 = os.path.join(currentpath,foldername,filename_2)
    pure, kij = get_pcpar.get_par(excelstr,T=temp)
    napvp = Mixture(*pure,dikij=kij)
    if Dw == 0:
        napvp.idx = 1
        napvp.idxp = 0
        napvp.idxa = 1
        nScan=100
    elif Dw > 0 and len(excelstr)<3:
        napvp.idx = 1
        napvp.idxp = 1
        napvp.idxa = 0
        nScan=100
    else:
        napvp.idx = 2
        napvp.idxp = 0
        napvp.idxa = 1
        nScan=30
        
    
    if not os.path.isfile(pathname_1):
        napvp.pH2OLV = np.asarray([1.013E5])
        
        result_napvp_sle = [napvp.SLE(psys=1.013E5,T=temp,wwASD=0)] 

        napvp.WriteToExcel([result_napvp_sle,result_napvp_sle],["sle_napvp","sle_napvp"],Datatype="sle_napvp")
    
    else: 
        None
    
    if not os.path.isfile(pathname_2):
        napvp.pH2OLV = np.asarray([1.013E5])
        result_napvp_gammabi = napvp.CalcSpaceDrugload(psys=napvp.pH2OLV,T=temp,n=nScan)
    
    else:
        None

    # The calculated PC-SAFT-data is safed in an Excel-File and here transfered to a dataframe for further processing.
    df_eq = pd.read_excel(pathname_1,tablename_1)
    df_scan = pd.read_excel(pathname_2,tablename_2)
    df_scan = df_scan.sort_values("wi"+str(napvp.idxa)) if len(excelstr)<3 else df_scan # makes sure everything is ascending for interpolant
    
    mu_s_SAFT_eq = df_eq.loc[0, "mui"+str(napvp.idxa)]
    mu_la_SAFT = np.fmax(df_scan.loc[:, "mui"+str(napvp.idxa)].values,mu_s_SAFT_eq)#DB:might get to -inf
    w_la_SAFT = df_scan.loc[:, "wi"+str(napvp.idxa)].values 
    
    if len(excelstr)>2:
        dl_SAFT=np.linspace(0,1,nScan+1)
        ww_SAFT=np.linspace(0,1,nScan+1)
        mu_la_SAFT=np.reshape(mu_la_SAFT,(nScan+1,nScan+1)).T.flatten() 
        # w1,w2=np.meshgrid(dl_SAFT,ww_SAFT)
        # from mpl_toolkits.mplot3d import Axes3D
        # fig = plt.figure()
        # ax = plt.axes(projection="3d")
        # ax.scatter(w1,w2,mu_la_SAFT)
        dmu_int = cs.interpolant("mu_int","bspline",[dl_SAFT,ww_SAFT],mu_la_SAFT-mu_s_SAFT_eq) #2D interpolation
    else:
        dmu_int = cs.interpolant("mu_int","linear",[w_la_SAFT],mu_la_SAFT-mu_s_SAFT_eq)

    return dmu_int,excelstr

def cryst(rho, M, temp, sigma, kt, g, D, Dw=0):
    
    # tend = 9.6E4
    # tsteps = 1200
    tend = 78950.4
    tsteps = 1000
    
    # Defining all constants:
    exp = cs.exp
    ln = cs.log
    pi = np.pi
    kB = sp.constants.Boltzmann
    NA = sp.constants.Avogadro   
    R = kB*NA
    
    # Crystallization model parameters:
    v_a_0 = M/rho #Unit: m^3/mol
    c_a_0 = 1/v_a_0 #Unit: mol/m^3
    # m_ges_0 = 0.1 #Unit: kg
    m_ges_0= 3.12712E-06
    #m_ges_0=4.78189353E-6
    ww0 = 0
    dl =0.68
    #dl = 1
    m_la_0 = dl*m_ges_0
    mp = m_ges_0 - m_la_0
    w_la_0 = m_la_0/m_ges_0
    Vs = m_ges_0/rho #Unit: m^3 
    
    a_sphere = (36*pi*v_a_0**2)**(1/3) #Unit: m^(2/3)
    a = a_sphere
    
    # Defining all symbolic casadi variables:
    dmu_sym = cs.MX.sym("dmu_sym")
    n_sym = cs.MX.sym("n_sym")
    nstar_sym = cs.MX.sym("nstar_sym")
    N = cs.MX.sym("N")
    N_strich = cs.MX.sym("N_strich")
    r = cs.MX.sym("r")
    r_strich = cs.MX.sym("r_strich")
    T = cs.MX.sym('T',1)   
    scale_r=(m_ges_0/0.1)**0 
    scale_N=(m_ges_0/0.1)**0  
    r_exp = scale_r*r
    N_exp = scale_N*N
    polstr,apistr,solstr = "pvpva64","naproxen","water"
    
    
    dmu_int,excelstr=PCSAFT(temp,polstr,apistr,solstr,Dw,dl)
    
    dmu_sla_0 = dmu_int(cs.vertcat(w_la_0,ww0)) if len(excelstr)>2 else dmu_int(w_la_0)

    r0 = ((2*sigma*v_a_0)/(R*temp*dmu_sla_0))/scale_r
    
    # Defining the main variable 'm_sa'which corresponds to the crystal mass (s) of the API (a). 'm_sa' is then used to calculate alpha.
    # Depening on the shape of the forming crystals, the Volume of 'm_sa_strich' needs to be changed.
    
    # Needle crystals
    m_sa_strich =rho*pi*(r0**2)*(r_strich)*N_strich*NA  #Unit: kg*molecule/mol
    
    # Spherical crystals
    # m_sa_strich =rho*(4/3)*pi*(r_strich**3)*N_strich*NA
    
    m_sa_fun = cs.Function("ms_a_fun",[N_strich,r_strich],[m_sa_strich])
    m_sa = m_sa_fun(N_exp,r_exp)
    alpha = m_sa/m_ges_0

    w_la_w=np.asarray([0])
    if Dw>0:
        #####
        #Wenn die Masse steigt
        durchmesser = 0.0145
        A=pi*(durchmesser/2)**2
        l=m_ges_0/A/rho 

        #mw_infty = 2E-7
        cw_infty=0.017/(1-0.017)#mw_infty/m_ges_0

        beta = crank.crank(T,Dw,l)
        beta = 1
        kt = beta*kt
        D = beta*D
        if len(excelstr)>2:
            
            dlvec_exp,w_la_w=ExtractIsotherm()
            cw_infty_exp=w_la_w/(1-w_la_w)
            cw_infty_fun = cs.interpolant("cw_infty_fun",'bspline',[dlvec_exp[::-1]],cw_infty_exp[::-1])
            
            dl_la = (m_la_0-m_sa)/(m_la_0-m_sa+mp)
            cw_la = cw_infty_fun(dl_la)
            
        ######
            #mw_infty = m_ges_0*cw_la #Wenn die Masse steigt   
            cw=beta*cw_la*(1-alpha)
        else:
            #mw_infty = m_ges_0*cw_infty #Wenn die MAsse nicht steigt
            cw=beta*cw_infty*(1-alpha)
        
        ww=cw/(1+cw)
        ww_fun = cs.Function("ww_fun",[N,r,T],[ww])
        mw =cw*m_ges_0 #
    else:
        mw=0
        ww_fun=cs.Function("ww_fun",[N,r,T],[0])
    
    m_la=m_la_0 - m_sa
    w_la = m_la/(m_la + mp + mw) #DB Muss sich auf die flüssige Phase beziehen. Sonst ist die Übersättigung im 1 Komponenten Fall nicht konstant
    dmu_sla =dmu_int(cs.vertcat(dl_la,ww)) if len(excelstr)>2 else dmu_int(w_la) #Variable Übersättigung "springt" rauf und runter
        
    n_la =  m_la/M
    c_la = n_la 
    
    nstar = (8/27)*(a*sigma/(dmu_sla*R*temp))**3
    nstar_fun = cs.Function("nstar_fun",[N,r,T],[nstar])
    S = cs.exp(dmu_sla)
    S_fun = cs.Function("S_fun",[N,r,T],[S])
    #mu_fun = cs.Function("mu_fun",[N,r,T],[mu_la])
    dmu_fun = cs.Function("dmu_fun",[N,r,T],[dmu_sla])
    
    alpha_fun = cs.Function("alpha_fun",[N,r,T],[alpha])
    
    W_sym = -n_sym*dmu_sym*R*temp + a*sigma*n_sym**(2/3)
    W_fun = cs.Function("Wfun",[n_sym,dmu_sym],[W_sym])
    dWdn_sym = cs.gradient(W_sym,n_sym)
    dWdn_sym_2 = cs.gradient(dWdn_sym,n_sym)
    dWdn_fun = cs.Function("dWdn_fun",[n_sym],[dWdn_sym])
    dWdn_2_fun = cs.Function("dWdn_fun",[n_sym],[dWdn_sym_2])
    Wstar = W_fun(nstar,dmu_sla)
    dWstardnstar = dWdn_fun(nstar)
    dWstardnstar_2 = dWdn_2_fun(nstar)
    
    c_fstar_dif = (48*(pi**2)*v_a_0)**(1/3)
    fstar_dif_sym = c_fstar_dif*D*c_la*(nstar_sym**(1/3))
    fstar_dif_fun = cs.Function("fstar_dif_fun",[nstar_sym],[fstar_dif_sym])
    fstar_dif = fstar_dif_fun(nstar)
    fstar_dif_fun2 = cs.Function("fstar_dif_fun2",[N,r,T],[fstar_dif])
    
    c_fstar_inter = (6*(pi**2)*v_a_0)**(1/3)
    fstar_inter_sym = c_fstar_inter*D*c_la*(nstar_sym**(2/3))
    fstar_inter_fun = cs.Function("fstar_inter_fun",[nstar_sym],[fstar_inter_sym])
    fstar_inter = fstar_inter_fun(nstar)
    fstar_inter_fun2 = cs.Function("fstar_inter_fun2",[N,r,T],[fstar_inter])
    
    z = (-dWstardnstar_2/(2*pi*R*temp))**(1/2)
    z_fun = cs.Function("z_fun",[N,r,T],[z])
    
    # Defining the ODE System
    c_r = kt*M/rho
    # f*_inter
    # dNdt = cs.fmax(1/scale_N*z*fstar_inter*c_a_0*exp(dmu_sla)*exp(-Wstar/(R*temp)),0)
    # f*_dif
    dNdt = cs.fmax(1/scale_N*z*fstar_dif*c_a_0*exp(dmu_sla)*exp(-Wstar/(R*temp)),0)
    # Power law
    drdt = cs.fmax(1/scale_r*c_r*(dmu_sla)**g,0)
    # BaS
    # drdt = cs.fmax(1/scale_r*kt*((dmu_sla)**(5/6))*exp(-g/dmu_sla),0)

    
    ode=cs.vertcat(dNdt,drdt)
    
    # In case alpha exceeds the drugload (dl), more API would have crystallized than there was initially. Thus, a casadi conditional sets both
    # dN/dt and dr/dt to zero as soon as alpha exceeds dl.
    ode = cs.conditional(alpha<(dl),[cs.vertcat(0,0)],cs.vertcat(dNdt,drdt))
    X = cs.vertcat(N,r)
    dummy = cs.MX.sym("dummy")
   
    ode_def = {'x':X,'t':T,'ode':ode}   
    
    dNdt_fun = cs.Function("dNdt_fun",[N,r,T],[dNdt])
    drdt_fun = cs.Function("drdt_fun",[N,r,T],[drdt])
    
    # Intitial Conditions Definition
    N0=1/NA/scale_N#ms0/M
    # r0 = ((2*sigma*v_a_0)/(R*temp*dmu_sla_0))/scale_r
    N0 = cs.DM(N0)
    r0 = cs.DM(r0)
    
    # DGL Solving Properties
    ts = np.linspace(0,tend,tsteps+1)
    x0 = [N0.full().T[0], r0.full().T[0]]
    opts = {}

    opts ['abstol'] = 1e-6
    opts['reltol'] = 1e-6
    opts['grid'] = ts
    # opts['max_num_steps'] = 10
    #opts['Timeout']=200
    opts['output_t0'] = True

    F = cs.integrator('sim', 'idas', ode_def, opts)
    sol = F(x0=x0)
    sol = sol['xf'].full().T
    N_vec = sol[:,0]
    r_vec = sol[:,1]
    r_vec_m = r_vec*scale_r
    r_vec_mum = r_vec*(10**6)*scale_r
    r_vec_nm = r_vec*(10**9)*scale_r
    
    # Defining solution vectors of all desired outputs
    alpha_vec = alpha_fun(N_vec,r_vec,ts).full().T[0]
    S_vec = S_fun(N_vec,r_vec,ts).full().T[0]
    #mu_vec = mu_fun(N_vec,r_vec,ts).full().T[0]
    dmu_vec = dmu_fun(N_vec,r_vec,ts).full().T[0]
    dNdt_vec = dNdt_fun(N_vec,r_vec,ts).full().T[0]*NA
    drdt_vec = drdt_fun(N_vec,r_vec,ts).full().T[0]*10**9
    nstar_vec = nstar_fun(N_vec,r_vec,ts).full().T[0]*NA
    z_vec = z_fun(N_vec,r_vec,ts).full().T[0]*NA**(-1/2)
    fstar_dif_vec = fstar_dif_fun2(N_vec,r_vec,ts).full().T[0]*NA
    fstar_inter_vec = fstar_inter_fun2(N_vec,r_vec,ts).full().T[0]*NA
    ww_vec=ww_fun(N_vec,r_vec,ts).full().T[0]
    N_vec=N_vec*scale_N*NA
    r_vec=r_vec*scale_r
    t_vec = ts/60 # in min
    
    
    cryst_output_dic = {}
    pa_cryst_output_keys = ['t_vec','N_vec','r_vec_mum','alpha_vec','S_vec','dmu_vec','dNdt_vec','drdt_vec','nstar_vec','z_vec','fstar_dif_vec','fstar_inter_vec']
    aw_cryst_output_keys = ['t_vec','N_vec','r_vec_mum','alpha_vec','S_vec','dmu_vec','dNdt_vec','drdt_vec','nstar_vec','z_vec','fstar_dif_vec','fstar_inter_vec','ww_vec']
    
    if Dw == 0:
        for i in pa_cryst_output_keys:
            cryst_output_dic.update({i:eval(i)})
    elif Dw > 0:
        for i in aw_cryst_output_keys:
            cryst_output_dic.update({i:eval(i)})
    
    return cryst_output_dic

if __name__ == "__main__":
    
    from scipy.interpolate import InterpolatedUnivariateSpline
    df_exp = pd.read_excel("exp_cryst_data_nap.xlsx","Sheet1") 
    df_t_exp = df_exp['time[min]']
    t_exp_arr = np.asarray(df_t_exp)
    t_vec_exp_raw = t_exp_arr.T
    t_vec_exp = df_exp['time[min]'].dropna().values
    df_y_exp = df_exp['alphaExp[-]']
    y_exp_arr = np.asarray(df_y_exp)
    y_vec_exp_raw = y_exp_arr.T
    y_vec_exp = df_exp['alphaExp[-]'].dropna().values
    
    def crystfit(t_vec_exp,sigma,kt,g,D):
        rho = 1200
        M = 0.23
        temp = 298.15
        out = cryst(rho, M, temp, sigma, kt, g, D, Dw=1E-13)
        tsim,alphasim = out['t_vec'],out['alpha_vec']
        alpha_fun = InterpolatedUnivariateSpline(tsim, alphasim)
        return alpha_fun(t_vec_exp)
    
    # fit sigma,kt,g,D 
    
    # BaS + f*_inter
    # popt,pcov=curve_fit(crystfit,t_vec_exp,y_vec_exp,p0=[0.004,1*10**-6,10,10**-9],bounds=([0.003,1*10**-10,0.1,5*10**-13],[0.005,1*10**-4,10**3,1*10**-8]))
    # print(popt)
    # popt = [4.07286572e-03 9.43590091e-07 1.90042559e+01 1.15505801e-09]
    
    # BaS + f*_dif
    # popt,pcov=curve_fit(crystfit,t_vec_exp,y_vec_exp,p0=[0.004, 5e-08, 5, 5e-17],bounds=([0.001,1*10**-9,1,1*10**-18],[0.02,1*10**-7,10**3,9*10**-16]))
    # print(popt)
    # popt = [3.19035369e-03, 5.61940583e-08, 2.19045109e+01, 2.21937506e-16]
    
    # Power law + f*_inter
    # popt,pcov=curve_fit(crystfit,t_vec_exp,y_vec_exp,p0=[0.005,1*10**-6,2.5,3*10**-9],bounds=([0.004,8*10**-7,2,9*10**-11],[0.006,1*10**-5,3,1*10**-8]))
    # print(popt)
    # popt = [4.19E-03,1.03E-06,2.11E+00,7.53E-10]

    # Power law + f*_diff
    popt,pcov=curve_fit(crystfit,t_vec_exp,y_vec_exp,p0=[0.004,5*10**-9,2,5*10**-17],bounds=([0.001,1*10**-10,1,1*10**-18],[0.01,1*10**-7,6,1*10**-15]))
    print(popt)
    # Fit One
    # popt = [3.60265454e-03,2.54601218e-08,2.62022881e+00,1.06388186e-16]
    # Fit Two 27 müm
    # popt = [4.71222829e-03 8.50686280e-08 3.74572157e+00 7.43239621e-18]
    # Fit Three 8.22 müm
    # popt = [4.71222829e-03, 2.50686280e-08, 3.74572157e+00, 2.43239621e-17]
    # Fit Three 4.57 müm
    # popt = [4.71222829e-03, 1.40686280e-08, 3.74572157e+00, 4.43239621e-17]
    #Fit Four 
    #popt=[0.00451222829, 1.8068628e-08, 3.74572157, 4.43239621e-17]
    #popt = [4.51222829e-03, 1.80686280e-08, 3.74572157e+00, 4.43239621e-17]
    # popt = [0.005512228,1.650686E-08,3.72572157,4.624E-17,1.45082147e-13]
    # popt = [5.11e-03, 1.6e-08, 3.74e+00, 5.3e-17]
    # popt = [0.005512228,1.650686E-08,3.72572157,4.624E-17]
    #popt=np.asarray([4.99351746e-03, 8.21705026e-09, 3.99364705e+00, 8.21947580e-17])
    cryst_output_dic_1=cryst(1200, 0.23, 298.15, *popt,Dw=1E-13) 
    # cryst_output_dic_1=cryst(1200, 0.23, 298.15, 0.004, 3E-8, 2.7, 1E-16,100)
    # cryst_output_dic_3=cryst(1200, 0.23, 298.15, 0.0035, 2.2E-9, 2.8, 1.2E-15, 100)
    # cryst_output_dic_4=cryst(1200, 0.23, 298.15, 0.0035, 2.2E-9, 2.9, 1.3E-15, 100)
    # cryst_output_dic_5=cryst(1200, 0.23, 298.15, 0.0035, 2.2E-9, 3.0, 1.4E-15, 100)             
    print(cryst_output_dic_1["r_vec_mum"])
    print(cryst_output_dic_1["S_vec"])
    fig,ax=plt.subplots()
    # ax.plot(cryst_output_dic["t_vec"],cryst_output_dic["alpha_vec"],'r-')
    ax.plot(cryst_output_dic_1["t_vec"],cryst_output_dic_1["alpha_vec"],'r-')
    # ax.plot(cryst_output_dic_2["t_vec"],cryst_output_dic_2["alpha_vec"],'b-')
    # ax.plot(cryst_output_dic_3["t_vec"],cryst_output_dic_3["alpha_vec"],'r-')
    # ax.plot(cryst_output_dic_4["t_vec"],cryst_output_dic_4["alpha_vec"],'r-')
    # ax.plot(cryst_output_dic_5["t_vec"],cryst_output_dic_5["alpha_vec"],'r-')
    ax.plot(t_vec_exp,y_vec_exp,'kx', linewidth=2.5)
    plt.show()
    
    # alphafit = cryst_output_dic_1["alpha_vec"]
    # t_fit = cryst_output_dic_1["t_vec"]
    # alphafit_fun = InterpolatedUnivariateSpline(t_fit, alphafit)
    # crank_df = pd.DataFrame({'t_exp':t_vec_exp})
    # crank_df['alpha_exp'] =  y_vec_exp.tolist()
    # crank_df['alpha_fit'] =  alphafit_fun(t_vec_exp).tolist()
    # crank_df['sigma'] = ''
    # crank_df['sigma'][0] = popt[0]
    # crank_df['kt'] = ''
    # crank_df['kt'][0] = popt[1]
    # crank_df['g'] = ''
    # crank_df['g'][0] = popt[2]
    # crank_df['D'] = ''
    # crank_df['D'][0] = popt[3]
    # crank_df['l_max'] = ''
    # crank_df['l_max'][0] = cryst_output_dic_1["r_vec_mum"][-1]
    # writer = pd.ExcelWriter('crystfit_dg68_f71.xlsx')
    # crank_df.to_excel(writer,'fitting_data')#,float_format='%.5f')
    # writer.save()

    # from tkinter.filedialog import askopenfilename
    # t1,w1,t2,w2,t3,w3,m0=ExtractSorption(askopenfilename())
    # tend=(t3[-1]-t3[0])
    # ax.plot((t3-t3[0]),w3,'ko')
    
    # Power law
    # Good initial fit
    # cryst_output_dic=cryst(1200, 0.23, 298.15, 0.004, 3E-8, 2.7, 1E-16,100) 
    # Good initial and end fit
    # cryst_output_dic=cryst(1200, 0.23, 298.15, 0.0035, 2E-9, 2.6, 1.2E-15, 100)

                