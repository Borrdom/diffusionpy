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
from PyCSAFT.get_pcpar import get_par
import PyCSAFT.get_pcpar
from PyCSAFT.PyCSAFTNil import Mixture

import crank
from ExtractWeidemannData import ExtractIsotherm,ExtractSorption


def PCSAFT(temp,polstr,apistr,solstr,Dw,dl):
    # polstr,apistr,solstr="pvpva64","lactose","water"
    # Dw=0

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
    pure, kij = get_par(excelstr,T=temp)
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
        nScan=100#30


    if not os.path.isfile(pathname_1):
        napvp.pH2OLV = np.asarray([1.013E5])

        result_napvp_sle = [napvp.SLE(psys=1.013E5,T=temp,wwASD=0)] #DB Ich dachte, dass sich im ternären das chem. Pot ändern muss da sich die Löslichkeit mit der Wasserkonzentration ändern!
        #DB Interessanterweise bleibt dieses aber immer konstant

        napvp.WriteToExcel([result_napvp_sle,result_napvp_sle],["sle_napvp","sle_napvp"],Datatype="sle_napvp")
    else:
        None

    if not os.path.isfile(pathname_2):
        # wbin=np.linspace(0,1,100) #new function call for the PCSAFT Scan

        napvp.pH2OLV = np.asarray([1.013E5])
        result_napvp_gammabi = napvp.CalcSpaceDrugload(psys=napvp.pH2OLV,T=temp,n=nScan)
        #napvp.WriteToExcel([result_napvp_gammabi,result_napvp_gammabi],["gammabi_napvp","gammabi_napvp"],Datatype="gammabi_napvp")
    else:
        None

    # The calculated PC-SAFT-data is safed in an Excel-File and here transfered to a dataframe for further processing.
    df_eq = pd.read_excel(pathname_1,tablename_1)
    df_scan = pd.read_excel(pathname_2,tablename_2)
    df_scan = df_scan.sort_values("wi"+str(napvp.idxa)) if len(excelstr)<3 else df_scan # makes sure everything is ascending for interpolant

    mu_s_SAFT_eq = df_eq.loc[0, "mui"+str(napvp.idxa)]
    #w_la_SAFT_eq = df_eq.loc[0, "wi"+str(napvp.idxa)]
    mu_la_SAFT = np.fmax(df_scan.loc[:, "mui"+str(napvp.idxa)].values,mu_s_SAFT_eq)#DB:might get to -inf
    #mu_la_SAFT =df_scan.loc[:, "mui"+str(napvp.idxa)].values
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
        dmu_int = cs.interpolant("mu_int","linear",[dl_SAFT,ww_SAFT],mu_la_SAFT-mu_s_SAFT_eq) #2D interpolation
    else:
        dmu_int = cs.interpolant("mu_int","linear",[w_la_SAFT],mu_la_SAFT-mu_s_SAFT_eq)

    return dmu_int,excelstr

def cryst(rho, M, temp, sigma, kt, g, D, N,r,T,Dw=0,m_ges_0=3.127118983E-6):

    tend = 89040.64959
    tsteps = 450

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
    #m_ges_0 = 0.1 #Unit: kg
    #m_ges_0=1.8E-6
    #m_ges_0=4.46199942E-6 #55
    #m_ges_0=3.127118983E-6 #71
    #m_ges_0=6.8658E-6#76
    #m_ges_0=6.2984E-6 #103

    ww0 = 0
    dl =0.68
    #dl = 1
    m_la_0 = dl*m_ges_0
    mp = m_ges_0 - m_la_0
    w_la_0 = m_la_0/m_ges_0
    Vs = m_ges_0/rho #Unit: kg/m^3 #Eventuell mit mges ersetzen, da Volumen steigt

    #Vs = m_ges_0/rho#Das darf sich glaube ich nicht auf den gesamte phase beziehen
    a_sphere = (36*pi*v_a_0**2)**(1/3) #Unit: m^(2/3)
    a = a_sphere



    # Defining all symbolic casadi variables:
    dmu_sym = cs.MX.sym("dmu_sym")
    n_sym = cs.MX.sym("n_sym")
    nstar_sym = cs.MX.sym("nstar_sym")

    N_strich = cs.MX.sym("N_strich")

    r_strich = cs.MX.sym("r_strich")

    # scale_r =1E5
    # scale_N =1E-4
    scale_r=(m_ges_0/0.1)**0 #DB: Das was ich bei Whatsapp meinte. Du musst das variable mit der Masse Skalieren
    scale_N=(m_ges_0/0.1)**0 #DB: Denk dir das so. Der SweetSpot war bei mges=0.1 hier waren rscale und Nscale 0.
    r_exp = scale_r*r# DB: Dementsprechend guckst du dir an wie die beiden mit mges skalieren
    N_exp = scale_N*N
    polstr,apistr,solstr = "pvpva64","naproxen","water"



    dmu_int,excelstr=PCSAFT(temp,polstr,apistr,solstr,Dw,dl)

    dmu_sla_0 = dmu_int(cs.vertcat(w_la_0,ww0)) if len(excelstr)>2 else dmu_int(w_la_0)



    r0 = ((2*sigma*v_a_0)/(R*temp*dmu_sla_0))/scale_r

    # Defining the two main casadi expressions 'ms' and 'w and their respective functions:
    m_sa_strich =rho*pi*(r0**2)*(r_strich)*N_strich*NA #rho*(4/3)*pi*(r_strich**3)*N_strich*NA #Unit: kg*molecule/mol
    #m_sa_strich =rho*(4/3)*pi*(r_strich**3)*N_strich*NA
    m_sa_fun = cs.Function("ms_a_fun",[N_strich,r_strich],[m_sa_strich])
    m_sa = m_sa_fun(N_exp,r_exp)
    alpha = m_sa/m_ges_0


    w_la_w=np.asarray([0])
    if Dw>0:
        #####
        #Wenn die Masse steigt
        durchmesser = 0.0145
        A=pi*(durchmesser/2)**2
        l=m_ges_0/A/rho #in Mikrometer-Range - realistisch? DB: Ja 15-40 micrometerfilme haben wir

        #mw_infty = 2E-7
        cw_infty=0.017/(1-0.017)#mw_infty/m_ges_0



        beta = crank.crank(T,Dw,l)

        # D=D*cs.exp((beta-1)*4)
        # kt=kt*cs.exp((beta-1)*4)
        D=D*beta
        kt=kt*beta

        if len(excelstr)>2:


            dlvec_exp,w_la_w=ExtractIsotherm()
            #w_la_wmittel=(w_la_w[-1]+w_la_w[0])/2
            #dl_la_mittel=(dlvec_exp[-1]+dlvec_exp[0])/2
            cw_infty_exp=w_la_w/(1-w_la_w)
            #dlvec_exp = np.asarray([0.67,0.4,0.2,0])
            #cw_infty_exp = (np.asarray([1.83E-6,1.87E-6,1.9E-6,1.94E-6])-m_ges_0)/m_ges_0
            from scipy.interpolate import InterpolatedUnivariateSpline,UnivariateSpline
            dlvec_vec=np.linspace(0,1,100)
            cw_infty_vec=InterpolatedUnivariateSpline(dlvec_exp[::-1],cw_infty_exp[::-1],k=3)(dlvec_vec)

            #cw_infty_fun = cs.interpolant("cw_infty_fun",'linear',[dlvec_exp[::-1]],cw_infty_exp[::-1])
            cw_infty_fun = cs.interpolant("cw_infty_fun",'linear',[dlvec_vec],cw_infty_vec)

            dl_la = (m_la_0-m_sa)/(m_la_0-m_sa+mp)
            cw_la = cw_infty_fun(dl_la)

            # Faktor=cw_la/cw_infty
            # gamma=Faktor*(1-alpha)
        ######
            #mw_infty = m_ges_0*cw_la #Wenn die Masse steigt
            cw=beta*cw_la*(1-alpha)
        else:
            #mw_infty = m_ges_0*cw_infty #Wenn die MAsse nicht steigt
            cw=beta*cw_infty*(1-alpha)


        #ausdruck = (mw_infty*beta)/(1+mw_infty*beta) #Irgendwas stimmt da noch nicht ganz
        #ww = (ausdruck*(1-alpha))/(1+ausdruck*alpha)
        ww=cw/(1+cw)
        ww_fun = cs.Function("ww_fun",[N,r,T],[ww])
        mw =cw*m_ges_0 #
        #mw=(ww*m_ges_0)/(1-ww)
        #mw=0
    else:
        mw=0
        ww_fun=cs.Function("ww_fun",[N,r,T],[0])
    m_la=m_la_0 - m_sa
    w_la = m_la/(m_la + mp + mw) #DB Muss sich auf die flüssige Phase beziehen. Sonst ist die Übersättigung im 1 Komponenten Fall nicht konstant
    dmu_sla =dmu_int(cs.vertcat(dl_la,ww)) if len(excelstr)>2 else dmu_int(w_la) #Variable Übersättigung "springt" rauf und runter




    #w_la = (m_la_0 - m_sa)/(m_la_0 + mp + mw)
    n_la =  m_la/M
    #c_la = n_la/Vs
    #m_l=(m_ges_0+ mw - m_sa)
    #V_l=m_l/rho
    c_la = n_la # DB:An sich kürzt sich das Vs in der ODE sowieso aus. Deswegen brauchen wir das gar nicht mehr

    # PC-SAFT is used to calculate the equilibrium chemical potemtial as well as a variety of chemical potentials belonging
    # to a specific mass fraction.
    # First, it is validated whether the data of the desired componen= has already been calculated by checking the
    # existence of the file.



    # DB: Übersättigung sinkt durch Wassersorption. Dies entspricht dem durch SAFT berechneten Wasserwert bei RH=1, da p>pH2OLV ist. Danach steigt diese wieder an bis zum LLEs
    #Eine sinkende Übersättigung dürfte kein Problem machen. Dennoch kann die Übersättigung im momentanten
    #Modell nicht abgebaut werden, da diese maximal zum Sorptionsgleichgewicht abbgebaut wird.
    # Dass heißt der Fall API/Wasser entspricht dem selben Problem wie beim Einkomponenten Fall
    # Die Übersättigung baut sich also nie ab. Sprich es gibt keine Bremse



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
    # TGGW=4500
    # D0=1E-100
    # D=cs.fmin((D-D0)/TGGW*T+D0,D)
    # D=cs.fmin((D-D0)*beta+D0,D)

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

    # kt0=1E-100
    # #kt=cs.fmin((kt-kt0)/TGGW*T+kt0,kt)
    # kt=cs.fmin((kt-kt0)beta+kt0,kt)
    c_r = kt*M/rho
    #dNdt = 1/scale_N*cs.fmax(z*fstar_dif*c_a_0*exp(dmu_sla)*exp(-Wstar/(R*temp))*Vs,0)
    dNdt = cs.fmax(1/scale_N*z*fstar_dif*c_a_0*exp(dmu_sla)*exp(-Wstar/(R*temp)),0)
    drdt = cs.fmax(1/scale_r*c_r*(dmu_sla)**g,0)
    ode=cs.vertcat(dNdt,drdt)

    #Schaltet eine künstliche Bremse für API/Wasser und API ein. Sprich Alpha kein nicht größer als die Drugload werden
    ode = cs.conditional(alpha<(dl),[cs.vertcat(0,0)],cs.vertcat(dNdt,drdt))#-w_la_SAFTeq
    X = cs.vertcat(N,r)
    dummy = cs.MX.sym("dummy")


    #alg_def = {'x':dummy,'z':X,'alg':cs.vertcat(ode[0],0),'ode':1}


    dNdt_fun = cs.Function("dNdt_fun",[N,r,T],[dNdt])
    drdt_fun = cs.Function("drdt_fun",[N,r,T],[drdt])

    # Intitial Conditions Definition
    N0=1/NA/scale_N#ms0/M
    r0 = ((2*sigma*v_a_0)/(R*temp*dmu_sla_0))/scale_r
    N0 = cs.DM(N0)
    r0 = cs.DM(r0)

    # DGL Solving Properties
    ts = np.linspace(0,tend,tsteps+1)
    # x0 = [N0,r0]
    x0 = [N0.full().T[0], r0.full().T[0]]
    opts = {}
    #opts['fsens_err_con'] = True
    #opts['quad_err_con'] = True
    opts ['abstol'] = 1e-6
    opts['reltol'] = 1e-6
    opts['grid'] = ts
    # opts['max_num_steps'] = 10
    #opts['Timeout']=200
    opts['output_t0'] = True
    # F = cs.integrator('sim', 'collocation', alg_def, {'grid':ts})
    # sol_r_scale = F(x0=0,z0=x0)
    # sol_r_scale = sol_r_scale['zf'].full().T
    xvec=X
    fvec=ode
    funs=[alpha_fun,S_fun,dmu_fun,dNdt_fun,drdt_fun,fstar_dif_fun2,fstar_inter_fun2,z_fun,drdt_fun,nstar_fun,ww_fun]
    return xvec,fvec,x0,opts,funs,T
def crystsolve(xvec,fvec,x0,opts,funs,T):#
    ts=opts['grid']
    alpha_fun=funs[0]
    S_fun=funs[1]
    dmu_fun=funs[2]
    dNdt_fun=funs[3]
    drdt_fun=funs[4]
    fstar_dif_fun2=funs[5]
    fstar_inter_fun2=funs[6]
    z_fun=funs[7]
    drdt_fun=funs[8]
    nstar_fun=funs[9]
    ww_fun=funs[10]
    X=xvec
    ode=fvec

    kB = sp.constants.Boltzmann
    NA = sp.constants.Avogadro
    R = kB*NA
    scale_N=1
    scale_r=1

    ode_def = {'x':X,'t':T,'ode':ode}
    F = cs.integrator('sim', 'idas', ode_def, opts)
    sol = F(x0=x0)
    sol = sol['xf'].full().T
    N_vec = sol[:,0]
    r_vec = sol[:,1]
    r_vec_m = r_vec
    r_vec_mum = r_vec*10**6
    r_vec_nm = r_vec*10**9

    # Defining solution vectors of all desired outputs
    alpha_vec = alpha_fun(N_vec,r_vec,ts).full().T[0]
    S_vec = S_fun(N_vec,r_vec,ts).full().T[0]
    #mu_vec = mu_fun(N_vec,r_vec,ts).full().T[0]
    dmu_vec = dmu_fun(N_vec,r_vec,ts).full().T[0]
    dNdt_vec = dNdt_fun(N_vec,r_vec,ts).full().T[0]*NA
    drdt_vec = drdt_fun(N_vec,r_vec,ts).full().T[0]*10**9
    nstar_vec = nstar_fun(N_vec,r_vec,ts).full().T[0]*NA
    z_vec = z_fun(N_vec,r_vec,ts).full().T[0]*NA**(-1/2)
    fstar_dif_vec = fstar_dif_fun2(N_vec,r_vec,ts).full().T[0]*NA#fstar_dif_fun(S_vec,nstar_vec).full().T
    fstar_inter_vec = fstar_inter_fun2(N_vec,r_vec,ts).full().T[0]*NA
    ww_vec=ww_fun(N_vec,r_vec,ts).full().T[0]
    N_vec=N_vec*scale_N*NA
    r_vec=r_vec*scale_r
    t_vec = ts


    cryst_output_dic = {}
    pa_cryst_output_keys = ['t_vec','N_vec','r_vec_nm','alpha_vec','S_vec','dmu_vec','dNdt_vec','drdt_vec','nstar_vec','z_vec','fstar_dif_vec','fstar_inter_vec']
    aw_cryst_output_keys = ['t_vec','N_vec','r_vec_nm','alpha_vec','S_vec','dmu_vec','dNdt_vec','drdt_vec','nstar_vec','z_vec','fstar_dif_vec','fstar_inter_vec','ww_vec']
    Dw=0
    if Dw == 0:
        for i in pa_cryst_output_keys:
            cryst_output_dic.update({i:eval(i)})
    elif Dw > 0:
        for i in aw_cryst_output_keys:
            cryst_output_dic.update({i:eval(i)})

    return cryst_output_dic

def crystfit(texp,sigma,kt,g,D):
    from scipy.interpolate import InterpolatedUnivariateSpline
    rho=1200
    M=0.23
    temp=298.15
    out=cryst(rho, M, temp, sigma, kt, g, D, Dw=0)
    tsim,alphasim=out['t_vec'],out['alpha_vec']
    alpha_fun=InterpolatedUnivariateSpline(tsim,alphasim)
    return alpha_fun(texp)

if __name__ == "__main__":
    N = cs.MX.sym("N")
    r = cs.MX.sym("r")
    T = cs.MX.sym('T',1)
    xvec,fvec,x0,opts,funs,T=cryst(1200, 0.23, 298.15, 5e-03, 8.2e-09, 4, 8.8e-17,N,r,T,1.45E-13)
    cryst_output_dic=crystsolve(xvec,fvec,x0,opts,funs,T)
    #cryst_output_dic=cryst(1200, 0.23, 298.15, 5e-03, 1E-9, 2, 4e-15,1.45E-13)
    print(cryst_output_dic["alpha_vec"])
    print(cryst_output_dic["ww_vec"])
    fig,ax=plt.subplots()
    fig1,ax1=plt.subplots()
    ax.plot(cryst_output_dic["t_vec"]/60,cryst_output_dic["ww_vec"])
    ax1.plot(cryst_output_dic["t_vec"]/60,cryst_output_dic["alpha_vec"])
    plt.show()
