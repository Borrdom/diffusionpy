# -*- coding: utf-8 -*-
"""
Created on Sat Apr  3 14:16:06 2021

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
#import Crank

def par_study():
    
    NA = sp.constants.Avogadro
    mges0 = 0.1
    
    # case = input('Please type in which information should be plotted (n, r, alpha, s, delta_mu, nt, rt, z, fstar_dif or fstar_inter): ')
    
    # All necessary parameters are defined and the corresponding excel sheet is transfered to a data frame (DF).
    # pars_org = ['M', 'rho', 'temp', 'kt', 'sigma', 'g', 'D']
    # pars = ['M', 'rho', 'temp', 'kt', 'sigma', 'g', 'D']
    # pars_org = ['M', 'rho', 'temp', 'kt', 'sigma', 'g', 'D','N0']
    # pars = ['M', 'rho', 'temp', 'kt', 'sigma', 'g', 'D','N0']
    pars_org = ['rho', 'M', 'temp','sigma', 'kt', 'g', 'D','N0']
    pars = ['rho', 'M', 'temp','sigma', 'kt', 'g', 'D','N0']
    
    # df = pd.read_excel("Case study prep.xlsx","Tabelle3")
    df = pd.read_excel("par_bounds.xlsx","values")
    
    # A for loop interates through the DF and safes the values belonging to the desired parameters in a dictionary. 
    dict_par_raw = {}
    dict_par = {}
    for i in pars:
        df_pars =  df.loc[df['Property'] == i]
        dict_par_raw = df_pars.to_dict('list')
        dict_par[i] = dict_par_raw
    
    # PC-SAFT is used to calculate the equilibrium chemical potemtial as well as a variety of chemical potential belonging to
    # a specific mass fraction

    
    # The created dictionary is modified in a way that the parameter names are the keys on the top level, the 
    # specific case (lower bound, upper bound, middle) the key on the lower level and belonging to this key the corresponding value. This value is
    # extracted from a list and transfered into a float.
    dict_par_org = {}
    for k0,v0 in dict_par.items():
        for index, (k1,v1) in enumerate(v0.items()):
            if index >= 1 :
                dict_par_org.setdefault(k0,{}).update({k1:float(v1[0])})
            else:
                True
    
        # To enable the user to keep certain parameters constant, the user has the possiblity to determine those parameters which are kept constant.
    while True:
        try:
            #par_input_raw = input("Please enter the parameter you would like to keep constant during the case study (out of: M,rho,temp,sigma,kt, D and g): ")
            par_input_raw="M" #DB the amount of inputs in your code is getting out of control
            par_input = par_input_raw.split()
        except ValueError:
            print('\n')
            print('The input was not appropriate. Please try again.')
            continue
        else:
            break   
# In a for loop those parameters which should be kept constant are appended with a "_c" in order to distinguish between constant parameters and parameters 
# which will be varied betweeen its physically reasonable lower and upper boundary. Then the user is asked to either enter a specific value for those parameters
# which should be kept constant or to decide if those parameters should be set to the lower (l) or upper (u) physically possible boundary or in between them (middel:m)
    pars_c = {}
    for i in par_input:
        if i in pars:
            pars[pars.index(i)] = i + "_c"
            var = i + "_c"
            while True:
                try:
                    #set_par = input("Please enter a value for the parameters " + var.upper() + " where C stands for constant or choose " + var.upper() + " to be in the middle (m), lower boundary (l) or upper boundary (u) of the physically possible range: ")
                    set_par="m"
                except ValueError:
                    print('\n')
                    print('The input was not appropriate. Please try again.')
                    continue
                else:
                    break 

            if set_par == 'm':
                locals()[var] = dict_par_org[i]['middle']
                pars_c[var] = locals()[var]
            elif set_par == 'l':
                locals()[var] = dict_par_org[i]['lower_bound']
                pars_c[var] = locals()[var]
            elif set_par == 'u':
                locals()[var] = dict_par_org[i]['upper_bound']
                pars_c[var] = locals()[var]
            else:
                locals()[var] = set_par
                pars_c[var] = locals()[var]

# The user can specify over which period of time the crystal mass fraction should be calculated (in hours). 
    # while True:
    #     try:
    #         tend_raw = float(input("Please enter the number of hours over which the crystal mass fraction should be calulated: "))
    #     except ValueError:
    #         print('\n')
    #         print('The input was not appropriate. Please try again.')
    #         continue
    #     else:
    #         break
    
    tend = 280
    
    # tend = tend_raw
    # tend = tend_raw*3600
    # if tend_raw%1 == 0:
    #     time_frame = str(int(tend_raw/1))
    # else:
    #     time_frame = str(tend_raw)

# The user has also the possibilty to specify how may parametersets should be genereated (later used as input to calculate the crystal mass fraction). 
    # while True:
    #     try:
    #         number_datasets = int(input("Please enter the number of datasets that should be generated: "))
    #     except ValueError:
    #         print('\n')
    #         print('The input was not appropriate. Please try again.')
    #         continue
    #     else:
    #         break      
    
    number_datasets = 4
    
    # A for loop iterates through the new list "pars" and if an element has not been replaced by a constant parameter it is checked whether the difference of 
    # magnitude of the lower and upper bound is higher than 3. To do so "Type Specifying" of 'format()' is used. 0 stands for the first argrument (value) and 1 
    # for the numbers after the dot. 'e' converts the value to the exponential notation with 10 numbers afer the dot. 
    # If the magnitude difference is higher than 3, the function "linspace" is used to create NOT evenly spaced values between the lower and
    # upper bound. They are NOT evenly spaced due to the high difference of the magnitudes. That is why first the log10() of the lower bound is given to
    # the function "linspace()". If the difference of magnitude of the lower and upper bound is equal or lower than 3 the function "linspace()" is used naturally.
    # If an element of pars has been replaced by a constant parameter an array is generated which only contains the constant value earlier defined by the user.
    # In general, all generated arrays are added to a matrix called "arr_pars".
    arr_pars = []
    for i in pars:
        if i in pars_org:
    
            if abs(int(("{0:.{1}e}".format(dict_par_org[i]['lower_bound'],10))[-3:])-int(("{0:.{1}e}".format(dict_par_org[i]['upper_bound'],10))[-3:])) > 3:
                lin_pars = 10**np.linspace(np.log10(dict_par_org[i]['lower_bound']),np.log10(dict_par_org[i]['upper_bound']),number_datasets+1)
                arr_pars.append(lin_pars)
            
            else:  
                lin_pars = np.linspace(dict_par_org[i]['lower_bound'],dict_par_org[i]['upper_bound'],number_datasets+1)
                arr_pars.append(lin_pars)
    
        else:
    
            lin_pars = np.linspace(float(locals()[i]),float(locals()[i]),number_datasets+1)
            arr_pars.append(lin_pars)    
    
    # Due to the fact that an increasing density decreases the molecular volume (v0) but an increasing molar mass INcreases the molecular volume 
    # a prerequisite needs to be made before calculating v0. The array of the density is reversed so that the array starts at the highest value
    # and ends at the lowest value. Afterwards, a for loop iterates through the reversed density array and saves the first value of the density
    # array (highest value) in the variable "rho_calc"; the corresponding molar mass "M_calc" equals the (not reversed) first value of the zero row of 
    # the array "arr_pars". Both are used to calculate v0. v0 is then appended to a list "v0_lst". This is repeated until all v0s are calculated.
    # This list is then converted into an array and appended to the array "arr_pars".
    rho_lst_reverse = reversed(arr_pars[0]) 
    
    v0_lst = []
    for index, i in enumerate(rho_lst_reverse):
        rho_calc = i
        M_calc = arr_pars[1][index]
        v0 = 1/(rho_calc*(1/M_calc)*NA)
        print('v0: '+str(v0))
        v0_lst.append(v0)
    
    arr_v0 = np.asarray(v0_lst)
    arr_pars.append(arr_v0)  
    
    # The matrix "arr_pars" is transposed so that now each row can serve as an input for the function "cryst()". The function "cryst()" needs seven arguments 
    # to calculate the crystal massfraction and to plot the result.
    # global A
    data_lst = [None]
    A = np.asarray(arr_pars).T
    m, n = A.shape
    data_lst = data_lst*len(range(0,m,1))
      
    # A for loop iterates through the just created matrix and uses one row with 7 columns (each of the 7 parameters) as input
    # for the function. Therefore, the row is converted into a list and passed to the function "cryst()" to calculate the
    # crystal fraction. Furthermore, an index is passed to the fuction due to saving the plots as PNGs and the index is for clean 
    # namening.
    for index,i in enumerate(range(0,m,1)):
        B = A[[i],:].tolist()
        print(*B[0],mges0,tend)
        data_lst[index] = cryst(*B[0],mges0,tend)
        # try:
        #     data_lst[index] = cryst(*B[0],mges0,tend)
        # except:
        #     True
    
    name_prep = []        
    for k,v in pars_c.items():
        name_prep.append(k) 

    name_prep.sort()
    
    name = ''
    for i in name_prep:
        name += str(i) + '=' + str(pars_c[i])
    
    # data_lst_max = []
    # for i in data_lst:
    #     data_lst_max.append(i[:,0].max())
  
    # Plotting Properties
    fig_lst = []
    ax_lst = []
    case_lst = ['N','r','S',r'\alpha',r'\Delta \mu',r'\dfrac{dN}{dt}',r'\dfrac{dr}{dt}','z','f^*_{dif}','f^*_{inter}']
    safename_lst = ['n','r','s','alpha','dmu','nt','rt','z','fdif','finter']
    unit_lst = ['-',r'\mu m','-','-','J','s^{-1}',r'\mu m \cdot s^{-1}','-','s^{-1}','s^{-1}']
    
    ts = np.linspace(1,tend,300)
    plot_font = {'fontname':'Calibri'}
    plt.rc('xtick',labelsize=18)
    plt.rc('ytick',labelsize=18)
    plt.ioff()
        
    for index1,val1 in enumerate(case_lst):
        fig,ax = plt.subplots()
        ax.set_ylabel(r'$\mathregular{ '+case_lst[index1]+r' } \ / \;'+unit_lst[index1]+r'}$',**plot_font,size=20,loc='top',rotation=0)
        ax.set_xlabel(r'$\mathregular{time \ / \; h}$',**plot_font,size=20)
    
        for index2,val2 in enumerate(data_lst):
            ax.plot(ts,val2[index1],'k-', linewidth=3, alpha=float(1/number_datasets*(index2)))

        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(1.5)
        
        ax.tick_params(direction='in', right=True, top=True, labelright=False, labeltop=False, width=1.5, length=7)
        ax.grid()
        plt.tight_layout()
        
        # savepath = r"C:\Users\Moritz\Documents\TU Dortmund\Master\Masterarbeit\Dominik_Kristallisation_Modellierung\Woche 22 (12.04.21-16.04.21)\Real Binary System\Plots - Parameter Study"
        savepath = r"C:\Users\Moritz\Documents\TU Dortmund\Master\Masterarbeit"
        savefile = 'mges0='+str(mges0)+name+'t='+str(tend)+'_'+str(safename_lst[index1])+'.png'
        fig.savefig(os.path.join(savepath, savefile),bbox_inches='tight')
        plt.show()
        
# def cryst(M, rho, temp, kt, sigma, g, D, v0, mges0, tend):
# def cryst(M, rho, temp, kt, sigma, g, D, N0, v0, mges0, tend):
def cryst(rho, M, temp, sigma, kt, g, D, N0, v0, mges0, tend):

    # Defining all constants:
    exp = cs.exp
    ln = cs.log
    pi = np.pi
    kB = sp.constants.Boltzmann
    NA = sp.constants.Avogadro   
    R=kB*NA
    
    # Crystallization Model Parameters
    v0 = M/rho
    c0 = 1/v0
    mges0 = 0.1
    dl = 0.67
    ml_a0 = dl*mges0
    mp = mges0 - ml_a0
    w0 = ml_a0/mges0
    Vs = mges0/rho
    a_sphere = (36*pi*v0**2)**(1/3)
    a = a_sphere
    
    # Defining all symbolic casadi variables
    dmu_sym = cs.MX.sym("dmu_sym")
    n_sym = cs.MX.sym("n_sym")
    nstar_sym = cs.MX.sym("nstar_sym")
    N = cs.MX.sym("N")
    N_strich=cs.MX.sym("N_strich")
    r = cs.MX.sym("r")
    r_strich=cs.MX.sym("r_strich")
    T = cs.MX.sym('T',1)   
    scale_r=1E-0
    scale_N=1E-0
    r_exp=scale_r*r
    N_exp=scale_N*N
    
    # Defining all casadi expressions and their respective functions
    ms_strich = rho*(4/3)*pi*(r_strich**3)*N_strich*NA
    ms_fun=cs.Function("ms_fun",[N_strich,r_strich],[ms_strich])
    ms=ms_fun(N_exp,r_exp)
    mw=0
    Dw=1E-13
    #mw=Crank.Crank(T,Dw)
    ws = (ml_a0 - ms)/(ml_a0 + mp+ mw)
    
    # PC-SAFT is used to calculate the equilibrium chemical potemtial as well as a variety of chemical potential belonging to
    # a specific mass fraction
    currentpath = os.getcwd()
    foldername = "PyCSAFTSheets"
    filename_1 = "sle_napvp_pvpva64_naproxen.xlsx"
    filename_2 = "GammaScan_pvpva64_naproxen.xlsx" 
    tabename_1 = "sle_napvp"
    tabename_2 = "GammaScan"# DB: changed to the new function call for the PCSAFT Scan
    filepath = os.path.join(currentpath,foldername)
    pathname_1 = os.path.join(currentpath,foldername,filename_1)
    pathname_2 = os.path.join(currentpath,foldername,filename_2)
    
    pure, kij=get_pcpar.get_par(["pvpva64","naproxen"],T=temp)
    napvp = Mixture(*pure,dikij=kij)
    napvp.idx=0
    napvp.idxp=1
    napvp.idxa=1
    
    if not os.path.isfile(pathname_1):
        napvp.pH2OLV= np.asarray([1.013E5])
        result_napvp_sle = [napvp.SLE(psys=1.013E5,T=temp)] 
        napvp.WriteToExcel([result_napvp_sle,result_napvp_sle],["sle_napvp","sle_napvp"],Datatype="sle_napvp")
    else: 
        None
    
    if not os.path.isfile(pathname_2):
        #wbin=np.linspace(0,1,100) #new function call for the PCSAFT Scan
        result_napvp_gammabi=napvp.CalcSpace(psys=napvp.pH2OLV,T=temp,n=100)
        #napvp.WriteToExcel([result_napvp_gammabi,result_napvp_gammabi],["gammabi_napvp","gammabi_napvp"],Datatype="gammabi_napvp")
    else:
        None

    df_eq = pd.read_excel(pathname_1,tabename_1)
    df_scan = pd.read_excel(pathname_2,tabename_2).sort_values("wi1") #makes sure everything is ascending for interpolant
    muSAFT_eq = df_eq.loc[0, "mui1"]

    nSAFT =  (ml_a0 - ms)/M
    cSAFT = nSAFT/Vs
    
    muSAFT = df_scan.loc[:, "mui1"].tolist()
    wSAFT = df_scan.loc[:, "wi1"].tolist()
    
    mu_int = cs.interpolant("mu_int","linear",[wSAFT],muSAFT)
   
    mu0 = mu_int(w0)
    dmu0 = mu0 - muSAFT_eq
    # SSAFT0 = aSAFT0/aSAFT_eq
    mu = mu_int(ws)
    dmu = mu - muSAFT_eq

    nstar = (8/27)*(a*sigma/(dmu*R*temp))**3
    nstar_fun = cs.Function("nstar_fun",[N,r],[nstar])
    S=cs.exp(dmu)
    S_fun = cs.Function("S_fun",[N,r],[S])
    dmu_fun = cs.Function("mu_fun",[N,r],[dmu])
    alpha = ms/mges0
    alpha_fun = cs.Function("alpha_fun",[N,r],[alpha])
    
    W_sym = -n_sym*dmu_sym*R*temp + a*sigma*n_sym**(2/3)
    W_fun = cs.Function("Wfun",[n_sym,dmu_sym],[W_sym])
    dWdn_sym = cs.gradient(W_sym,n_sym)
    dWdn_sym_2 = cs.gradient(dWdn_sym,n_sym)
    dWdn_fun = cs.Function("dWdn_fun",[n_sym],[dWdn_sym])
    dWdn_2_fun = cs.Function("dWdn_fun",[n_sym],[dWdn_sym_2])
    Wstar = W_fun(nstar,dmu)
    dWstardnstar = dWdn_fun(nstar)
    dWstardnstar_2 = dWdn_2_fun(nstar)
    
    c_fstar_dif = (48*(pi**2)*v0)**(1/3)
    fstar_dif_sym = c_fstar_dif*D*cSAFT*(nstar_sym**(1/3))
    fstar_dif_fun = cs.Function("fstar_dif_fun",[nstar_sym],[fstar_dif_sym])
    fstar_dif = fstar_dif_fun(nstar)
    fstar_dif_fun2=cs.Function("fstar_dif_fun2",[N,r],[fstar_dif])
    z = (-dWstardnstar_2/(2*pi*R*temp))**(1/2)
    z_fun = cs.Function("z_fun",[N,r],[z])
    
    # ce_sym = c0*exp(-W_sym/(kB*temp))
    # ce_fun = cs.Function("ce_fun",[n_sym,dmu_sym],[ce_sym])
    
    # Defining the ODE System
    c_r = kt*M/rho
    #fstar_dif=D
    dNdt = 1/scale_N*cs.fmax(z*fstar_dif*c0*exp(dmu)*exp(-Wstar/(R*temp))*Vs,0)
    drdt = 1/scale_r*cs.fmax(c_r*(dmu)**g,0)
    ode = cs.vertcat(dNdt,drdt)
    X = cs.vertcat(N,r)
    dummy = cs.MX.sym("dummy")

    
    #alg_def = {'x':dummy,'z':X,'alg':cs.vertcat(ode[0],0),'ode':1}   
    ode_def = {'x':X,'t':T,'ode':ode}   
    
    dNdt_fun = cs.Function("dNdt_fun",[N,r],[dNdt])
    drdt_fun = cs.Function("drdt_fun",[N,r],[drdt])
    
    # Intitial Conditions Definition
    #ms0 = (10**-10)*ml_a0
    # S0 = mges0/ms0
    #if (2*sigma*v0)/(kB*temp*ln(SSAFT0)) < (3*v0/(4*pi))**(1/3):

    #r0 = ((3*V0/(4*pi))**(1/3))
    N0=1/NA/scale_N#ms0/M
    #ms0=rho*(4/3)*pi*(r**3)*N*NA
    #r0=(ms0/(N0*NA*rho*(4/3)*pi))**(1/3)
    #else:
    r0 = ((2*sigma*v0)/(R*temp*dmu0))/scale_r
    
    
    #N0 = ms0/(rho*(4/3)*pi*(r0*r_scale**-1)**3)

    N0 = cs.DM(N0)
    r0 = cs.DM(r0)
    
    # print('N0: '+str(N0))
    # print('r0: '+str("{0:.{1}e}".format(r0,10)))
    # print('r0: '+str(r0*10**-13))
    # print('r0_geo: '+str(((3*v0/(4*pi))**(1/3))*10**9))
    # print('r0_th: '+str(((2*sigma*v0)/(kB*temp*ln(S0)))*10**9))
    
    # DGL Solving Properties
    ts = np.linspace(1,tend,300)
    # x0 = [N0,r0] 
    x0 = [N0.full().T[0], r0.full().T[0]]
    opts = {}
    opts['fsens_err_con'] = True
    opts['quad_err_con'] = True
    opts ['abstol'] = 1e-6
    opts['reltol'] = 1e-6
    opts['grid'] = ts
    # opts['max_num_steps'] = 10
    opts['output_t0'] = True
    # F = cs.integrator('sim', 'collocation', alg_def, {'grid':ts})
    # sol_r_scale = F(x0=0,z0=x0)
    # sol_r_scale = sol_r_scale['zf'].full().T

    F = cs.integrator('sim', 'idas', ode_def, opts)
    sol = F(x0=x0)
    sol = sol['xf'].full().T
    N_vec = sol[:,0]
    r_vec = sol[:,1]
    r_vec_m = r_vec
    r_vec_mum = r_vec*10**6
    
    alpha_vec = alpha_fun(N_vec,r_vec).full().T[0]
    S_vec = S_fun(N_vec,r_vec).full().T[0]
    dmu_vec = dmu_fun(N_vec,r_vec).full().T[0]#dmu_fun(N_vec,r_vec).full().T
    dNdt_vec = dNdt_fun(N_vec,r_vec).full().T[0]
    drdt_vec = drdt_fun(N_vec,r_vec).full().T[0]*10**6
    nstar_vec = nstar_fun(N_vec,r_vec).full().T[0]
    z_vec = z_fun(N_vec,r_vec).full().T[0]
    fstar_dif_vec = fstar_dif_fun2(N_vec,r_vec).full().T[0]#fstar_dif_fun(S_vec,nstar_vec).full().T
    fstar_inter_vec = fstar_dif_vec#fstar_inter_fun(S_vec,nstar_vec).full().T
    # ce_vec = ce_fun()
    N_vec=N_vec*scale_N
    r_vec=r_vec*scale_r
    
    return N_vec, r_vec_mum, S_vec, alpha_vec, dmu_vec, dNdt_vec, drdt_vec, z_vec, fstar_dif_vec, fstar_inter_vec


if __name__=="__main__":
    result=cryst(1200, 0.23, 293.15, 0.02, 3.16E-8, 1, 1E-11,1, 0.23/1200, 0.1, 280)
    par_study()            
# ListOfNamesOfInputs=['rho', 'M', 'temp',"sigma","kt","g","D","N0","v0","mges0","tend"]
# rho=1200
# M=0.23
# temp=293.15
# sigma=0.02
# kt=10**np.linspace(-4.5,-6.5,3)
# g=1
# D=10**np.linspace(-12,-10,3)
# v0=M/(rho)
# N0=1
# mges0=0.1 
# tend=280
                
                
                