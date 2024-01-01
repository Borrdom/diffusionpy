import os 
from os.path import join
cwd=os.getcwd()
import pandas as pd
pathsaft=join("C:\BöttcherPCSAFT","e-Dipol_PC-SAFT_DoppelT_kij(T)_KOMI mit JC_copol_D")
import subprocess
from subprocess import PIPE,STDOUT
from lmfit import Minimizer,Parameters,report_fit
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
import Thermoplot 
from Thermoplot import Symbol,Colorcode
import time
from PureComponentData import DenMol


def writeParametertxt(p,strp):
    """Writes Value and Name of Parameter in a textfile"""
    p1=join(cwd,"input_file")
    p2=join(pathsaft,"input_file")
    np.savetxt(join(p1,strp+".dat"),[p],fmt='%.6f')
    try:
        np.savetxt(join(p2,strp+".dat"),[p],fmt='%.6f')
    except:
        pass
        
def SetupPCSAFTVector(M1,T,p,*args):
    """Returns a PCSAFTvector of given Inputs for the PCSAFT communication"""
    if  hasattr(T, "__len__"):
        Tvec=T
    else:
        Tvec=np.ones_like(args[0])*T
    M1vec=np.ones_like(args[0])*M1*1000
    if  hasattr(p, "__len__"):
        pvec=p
    else:
        pvec=np.ones_like(args[0])*p
    x=[]
    x.extend((Tvec,M1vec,*args,pvec))
    return x


def PCSAFTCall(pathsaft,inputstr):
    """Communicates with Fortran PCSAFT over executable"""
    pathcall=join(pathsaft,"Debug","PCSAFT.exe")
    comm = subprocess.Popen(pathcall, stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    #inputstr='9\n2\nwater\npvp\nindometacin\nx\n''    
    try:
        grep_stdout = comm.communicate(input=inputstr.encode("utf-8") + b" ")[0]
    except:
        grep_stdout="PCSAFT.exe not found"
    return grep_stdout

def PCSAFT(params,x,*args):
    """Calls PCSAFT and returns gammas,ideal Solubility und Density for a given PCSAFT vector"""
    ncomp=len(args)
    uij = params['uij']
    kij = params['kij']
    mij = params['mij']
    sij = params['sij']
    kijS1 = params['kijS1']
    kijS2 = params['kijS2']
    matinp=np.vstack(x).transpose()
    inputstr='9\n2\n'
    txtstr=''
    for i,val in enumerate(args):
        inputstr+=val+'\n'
    for i,val in enumerate(args[1:]):
        txtstr+=val+'_' 
    inputstr+='x\n'
    if  len(args)==2:
        txtstr+="x.dat"
    else:
        txtstr=txtstr[:-1]+".dat"
    pathinp=join(cwd,"VLE_ATPS",txtstr)
    np.savetxt(pathinp,matinp,fmt='%.6f')
    writeParametertxt(uij,"uij")
    writeParametertxt(kij,"kij")
    writeParametertxt(mij,"mij")
    writeParametertxt(sij,"sij")
    writeParametertxt(kijS1,"kijS1")
    writeParametertxt(kijS2,"kijS2")
    #Öffne PCSAFT.exe
    report=PCSAFTCall(pathsaft,inputstr)
    #print(report)
    #Lese output
    pathout=join(cwd,"output_file","VLE_ATPS.dat")
    matout= np.loadtxt(pathout)
    rho=matout[:,ncomp-1+1]
    gamma2=matout[:,ncomp-1+3]
    gamma1=matout[:,ncomp+3]
    SLE=matout[:,-1]
    #print(gamma2)    
    return gamma2,SLE,rho,gamma1


def PCSAFTSLE(params,Tdata,*args):
    """Creates Value Table gamma w for given Temperatures and iterates over the SLE-Equation to return real Solubilities"""
    wL=[]
    for i,val in enumerate(Tdata):
        w2=np.linspace(0,1,100)
        w1=1-w2
        p=1.013
        _,M1=DenMol(args[1])
        SLEx=SetupPCSAFTVector(M1,val,p,w2,w1)
        SLE=PCSAFT(params,SLEx,*args)
        gamma2=SLE[0]
        va=SLE[1]    
        try:
            #plt.plot(w2,gamma2-va)
            wL.append(InterpolatedUnivariateSpline(w2[10:],(gamma2-va)[10:]).roots()[0])
            
        except:
            wL.append(1)
    wL=np.asarray(wL)  
    print(wL)
    
    return wL

def PCSAFTPVT(params,*args):
    w2=np.linspace(0,1,100)
    T=298.15
    _,M1=DenMol(args[1])
    p=1.013
    w1=1-w2
    PVTx=SetupPCSAFTVector(M1,T,p,w2,w1)
    rho=PCSAFT(params,PVTx,*args)[2]
    return w2,rho

def PVTFIT(params,*args):
    params['uij'].set(vary=1)
    params['kij'].set(vary=0)
    params['mij'].set(vary=1)
    params['sij'].set(vary=1)
    params['kijS1'].set(vary=0)
    params['kijS2'].set(vary=0) 
    
    cwd = os.getcwd()
    path1=join(cwd,"PVTbin.csv")
    Data=pd.read_csv(path1)
    Tdata=Data.values[0:,0]    
    pdata=Data.values[0:,1]
    rhodata=Data.values[0:,2]
    w2=np.zeros(rhodata.shape[0])
    T=Tdata
    _,M1=DenMol(args[1])
    p=pdata
    w1=1-w2
    PVTx=SetupPCSAFTVector(M1,T,p,w2,w1)
    

    def PVTbin(params,x,rhodata):
        rho=PCSAFT(params,PVTx,*args)[2]
        #print(rho)
        print(params)
        return rho-rhodata
    minner = Minimizer(PVTbin, params ,fcn_args=(PVTx, rhodata))
    result = minner.minimize(method='cobyla')
    report_fit(result)
    uijres=result.params["uij"]
    mijres=result.params["mij"]
    sijres=result.params["sij"]
    params['uij'].set(value=uijres,vary=0)
    params['mij'].set(value=mijres,vary=0) 
    params['sij'].set(value=sijres,vary=0)
    rho=PCSAFT(params,PVTx,*args)[2]
    fig,ax=plt.subplots()
    ax.plot(pdata,rhodata,'kx')
    ax.plot(pdata,rho,'rx')
    return params

def PCSAFTVLE(params,*args):
    """Creates Value Table gamma w for a given Temperature and calculates activity of component 2"""

    w2=np.linspace(0,1,100)
    T=298.15
    _,M1=DenMol(args[1])
    p=params['p']
    w1=1-w2
    VLEx=SetupPCSAFTVector(M1,T,p,w2,w1)
    _,M1=DenMol(args[1])
    _,M2=DenMol(args[0])
    def mass2mole(w1,w2):
        x2=w2/M2/(w1/M1+w2/M2)
        return x2
    x2=mass2mole(w1,w2)
    R=8.3145
    T=298.15
    rho02,M2=DenMol(args[0])
    pH2OLV=0.03166
    PI0i=np.exp(M2/rho02*(p-pH2OLV)*1E5/(R*T))
    print(PI0i)
    RH=PCSAFT(params,VLEx,*args)[0]*x2
    return w2,RH

def PCSAFTLLE(params,*args):
    """Creates Value Table gamma w for given Temperatures and iterates over the LLE-Equation to return Solubilities""" 
    
    w2=np.linspace(0,1,1000)
    w1=1-w2
    w3=np.zeros_like(w2)
    p=1.013
    M1=25700/1000
    T=298.15
    LLEx=SetupPCSAFTVector(M1,T,p,w2,w1)
    gamma2=PCSAFT(params,LLEx,*args)[0]
    gamma1=np.nan_to_num(PCSAFT(params,LLEx,*args)[3])
    
    M1=25700/1000
    #M1=357.88/1000
    #M2=18.015/1000
    M2=58.015/1000
    M3=25700/1000
    def mass2mole(w1,w2,w3):
        x2=w2/M2/(w1/M1+w2/M2+w3/M3)
        x3=w3/M3/(w1/M1+w2/M2+w3/M3)
        return x2,x3
    x2,x3=mass2mole(w1,w2,w3)
    x1=(1-x3-x2)
    a1_fun=InterpolatedUnivariateSpline(w2,gamma1*x1)
    a2_fun=InterpolatedUnivariateSpline(w2,gamma2*x2)
    
    TH22=np.gradient(np.log(gamma2),x2)*x2+1
    TH12=np.gradient(np.log(gamma2),x1)*x1
    
    wherenan=np.where(np.isnan(TH22))[0]
    THf22=np.delete(TH22,wherenan)
    w2ff=np.delete(x2,wherenan)
    TH_fun=InterpolatedUnivariateSpline(w2ff,THf22)
    
    TH_fun12=InterpolatedUnivariateSpline(w2,TH12)
    
    TH_roots=TH_fun.roots()
    TH_roots12=TH_fun12.roots()
    fig1,ax1=Thermoplot.Plot("wAceton","$\Gamma$22")
    fig2,ax2=Thermoplot.Plot("wAceton","$\Gamma$12")
    ax1.plot(w2ff,THf22)
    ax2.plot(w2,TH12)
    ax1.plot(TH_roots,TH_fun(TH_roots),'rx')
    ax2.plot(TH_roots12,TH_fun12(TH_roots12),'rx')
    #plt.plot(w2,TH21)
    #plt.plot(w2,TH11)
    def obj(p):
        w2O1,w2O2=p
        return np.abs(np.asarray([a1_fun(w2O1)-a1_fun(w2O2),a2_fun(w2O1)-a2_fun(w2O2)]))
    from scipy.optimize import least_squares
    x0=[TH_roots12[0],TH_roots[1]]
    lb=[0,TH_roots[1]]
    ub=[TH_roots12[0],1]
    res=least_squares(obj,x0,bounds=(lb,ub))
    print(res.x-x0)
    w2x1=res.x[0]
    w2x2=res.x[1]

    return res.x







# create a set of Parameters
params = Parameters()

#Bounds and Initials

lbui=200
ubui=450
lbkij=-0.8
ubkij=0
lbmi=0.005
ubmi=0.15
lbsi=2
ubsi=4
lbkijS1=-1
ubkijS1=1
lbkijS2=-0.1
ubkijS2=0.1
ui0=lbui
kij0=lbkij
mi0=lbmi
si0=lbsi
kijS10=lbkijS1
kijS20=lbkijS2
uistep=50
kijstep=0.05
mistep=0.01
sistep=0.1
kij2step=0.05
params.add('uij', value=ui0, min=lbui,max=ubui,brute_step=uistep)
params.add('kij', value=kij0,min=lbkij,max=ubkij,brute_step=kijstep)
params.add('mij', value=mi0, min=lbmi, max=ubmi,brute_step=mistep)
params.add('sij', value=si0,min=lbsi, max=ubsi,brute_step=sistep)
params.add('kijS1', value=kijS10,min=lbkijS1, max=ubkijS1,brute_step=kij2step)
params.add('kijS2', value=kijS20,min=lbkijS2, max=ubkijS2,brute_step=kij2step)
params['uij'].set(value=205,vary=1)
params['kij'].set(value=-0.148,vary=0)
params['mij'].set(value=0.04,vary=1)
params['sij'].set(value=2.7,vary=1)
params.add('p', value=100)
params['p'].set(value=1.013,vary=0)

#PVPVA64
#params['kijS1'].set(value=0.0922,vary=1)
#params['kijS2'].set(value=-0.000633,vary=1)
#params['p'].set(value=100,vary=1)
#params['uij'].set(value=205,vary=1)
#params['kij'].set(value=-0.146,vary=0)
#params['mij'].set(value=0.0415,vary=0)
#params['sij'].set(value=3.02,vary=1)
#params['uij'].set(value=400.28,vary=1)
#params['kij'].set(value=-0.208,vary=0)
#params['mij'].set(value=0.036032,vary=0)
#params['sij'].set(value=2.169,vary=1)
#params['p'].set(value=100,vary=1)
#params['kijS1'].set(value=-2,vary=0)
#params['kijS2'].set(value=-0.1,vary=0)
#kijS1=0.092180
#kijS2=-0.000633

#params['uij'].set(value=396.69,vary=1)
#params['kij'].set(value=-0.3248,vary=0)
#params['mij'].set(value=0.0159,vary=0)
#params['sij'].set(value=2.19,vary=1)
#Differential Evolution. Without SLE. [4:]


#params['uij'].set(value=389,vary=1)
#params['kij'].set(value=-0.145,vary=0)
#params['mij'].set(value=0.0402,vary=0)
#params['sij'].set(value=2.343,vary=1)
#Pauss Indo water kij [4:]

params['kijS1'].set(value=0.092180,vary=0)
params['kijS2'].set(value=-0.000633,vary=0)
#params['kijS2'].set(value=0,vary=0)


Pol="pvp"
API="indometacin"
Sol="water"
#LLE=PCSAFTLLE(params,Sol,API)
def VLEFIT(params,*args):
    start_time = time.time()
    params['uij'].set(vary=0)
    params['kij'].set(vary=1,value=-0.15)
    params['mij'].set(vary=1,value=0.05)
    params['sij'].set(vary=0)
    params['kijS1'].set(vary=0)
    params['kijS2'].set(vary=0) 
    #cwd = os.getcwd()
    path1=join(cwd,"VLEbin.csv")
    Data=pd.read_csv(path1)    
    RH=Data.values[0:,0]
    ww=Data.values[0:,1]
    w2=ww
    T=298.15
    _,M1vec=DenMol(args[1])
    p=1.013
    w1=1-w2
    w3=np.zeros_like(w1)
    VLEdata=RH
    VLEx=SetupPCSAFTVector(M1vec,T,p,w2,w1)
    _,M1=DenMol(args[1])
    _,M2=DenMol(args[0])
    def mass2mole(w1,w2):
        x2=w2/M2/(w1/M1+w2/M2)
        return x2
    x2=mass2mole(w1,w2)
    
    def VLEbin(params,x,data):
        gamma2=PCSAFT(params,x,*args)[0]
        return gamma2*x2-data
    minner = Minimizer(VLEbin, params ,fcn_args=(VLEx, VLEdata))
    result = minner.minimize(method='nelder')
    report_fit(result)
    kijres=result.params["kij"]
    mijres=result.params["mij"]
    params['uij'].set(vary=1)
    params['kij'].set(value=kijres,vary=0)
    params['mij'].set(value=mijres,vary=0)
    params['sij'].set(vary=1)
    params['kijS1'].set(vary=1)
    params['kijS2'].set(vary=0) 
    print(result.residual)
    report_fit(params)
    print("--- %s seconds ---" % (time.time() - start_time))
    return params

def SLEFIT(params,*args):
    params['uij'].set(vary=0)
    params['kij'].set(vary=0)
    params['mij'].set(vary=0)
    params['sij'].set(vary=0)
    params['kijS1'].set(vary=1,value=0.092180)
    params['kijS2'].set(vary=1,value=-0.000633) 
    #cwd = os.getcwd()
    path1=join(cwd,"SLEbin.csv")
    Data=pd.read_csv(path1)    
    Tdata=Data.values[:,0]
    wLexp=Data.values[:,1]
    print(wLexp)
    def SLEbin(params,x,data):    
        return PCSAFTSLE(params,x,API,Pol)-wLexp
    minner = Minimizer(SLEbin, params ,fcn_args=(Tdata, wLexp))
    result = minner.minimize(method='cobyla')
    report_fit(result)
    kijS1res=result.params["kijS1"]
    kijS2res=result.params["kijS2"]
    
    params['uij'].set(vary=1)
    params['kij'].set(vary=0)
    params['mij'].set(vary=0)
    params['sij'].set(vary=1)
    params['kijS1'].set(value=kijS1res,vary=1)
    params['kijS2'].set(value=kijS2res,vary=1) 
    print(result.residual)
    return params







def TotalFit(params,x,data):
    params['uij'].set(vary=0)
    params['kij'].set(vary=1)
    params['mij'].set(vary=1)
    params['sij'].set(vary=0)
    params['kijS1'].set(vary=0)
    params['kijS2'].set(vary=0)
    #Pol="pvpva64"
    Pol="pvp"
    API="indometacin"
    Sol="water"
    #cwd = os.getcwd()
#    path1=join(cwd,"VLEbin.csv")
#    Data=pd.read_csv(path1)
#    
#    RH=Data.values[:,0]
#    ww=Data.values[:,1]
#    w2=ww
#    T=298.15
#    M1vec=25700
#    p=1.013
#    w1=1-w2
#    w3=np.zeros_like(w1)
#    VLEdata=RH
#    VLEx=SetupPCSAFTVector(M1vec,T,p,w2,w1)
    _,M1=DenMol(Pol)
    _,M2=DenMol(Sol)
    _,M3=DenMol(API)
    def mass2mole(w1,w2,w3):
        x2=w2/M2/(w1/M1+w2/M2+w3/M3)
        x3=w3/M3/(w1/M1+w2/M2+w3/M3)
        return x2,x3
#    x2,x3=mass2mole(w1,w2,w3)
#    
#    VLEbin=lambda params,x,data:    PCSAFT(params,x,Sol,Pol)[0]*x2-data
#    minner = Minimizer(VLEbin, params ,fcn_args=(VLEx, VLEdata))
#    result = minner.minimize(method='nelder')
#    report_fit(result)
#
#    kijres=result.params["kij"]
#    mijres=result.params["mij"]
#    params['uij'].set(vary=1)
#    params['kij'].set(value=kijres,vary=0)
#    params['mij'].set(value=mijres,vary=0)
#    #params['kij'].set(value=-0.146,vary=0)
#    #params['mij'].set(value=0.04,vary=0)
#    params['sij'].set(vary=1)
#    report_fit(result)
    #params['kijS1'].set(vary=1)
    #params['kijS2'].set(vary=1) 
    #params=VLEFIT(params,Sol,Pol)
    print(params)
    params=SLEFIT(params,API,Pol)
    
    x2,x3=mass2mole(x[3],x[2],x[4])
    residual=(PCSAFT(params,x,Sol,Pol,API)[0]*x2-data)
    return residual
# do fit, here with the default least1sq algorithm

#params=PVTFIT(params,Sol,Pol)
path1=cwd+"\\ASDZografi.csv"
Data=pd.read_csv(path1)
#RH=Data.dropna().values[:,0][22:31]
RH=Data.dropna().values[:,0]
Zog=[]
ZogP=[]
lis=[1,2,3,4,5]
nASD=len(lis)
#nData=9
nData=31


for i in lis:
    W=Data.dropna().values[:,i+1]
    P=float(Data.columns[i+1])
    Zog.append(W)
    ZogP.append(P)
ww2=np.array([])
ww3=np.array([])
data=np.array([])

for i in range(nASD):
    #wreal=Zog[i][22:31]
    wreal=Zog[i]
    wPolyASD=ZogP[i]
    w2=wreal
    w3=((1-wPolyASD)*(1-wreal))
    ww2=np.hstack((ww2,w2))
    ww3=np.hstack((ww3,w3))
    data=np.hstack((data,RH))
ww1=1-ww2-ww3
ww2=np.asarray(ww2)
ww3=np.asarray(ww3)
T=np.ones_like(ww2)*298.15
M1vec=np.ones_like(ww2)*25700
params['p'].set(value=1.013)  
p=np.ones_like(ww2)*params['p']
x=[]
x.extend((T,M1vec,ww2,ww1,ww3,p))
VLEterndata=data



#Sol="acetone"
#wL=PCSAFTLLE(params,Sol,Pol)


#w2ac,RHac=PCSAFTVLE(params,Sol,Pol,API)
#fig22,ax22=Thermoplot.Plot(xla,yla)
#ax22.plot(RHac,w2ac,'k-')
#TotalFit(params,x,VLEterndata)
#PCSAFT SLE via Python

#minner = Minimizer(TotalFit, params ,fcn_args=(x, VLEterndata))
#result = minner.minimize(method='nelder')
#kijS2=-0.000633
#params['uij'].set(value=356,vary=1)
#params['kij'].set(value=-0.309,vary=0)
#params['mij'].set(value=0.0152,vary=0)
#params['sij'].set(value=2.22,vary=1)
#params['uij'].set(value=356,vary=1)
#params['kij'].set(value=-0.309,vary=0)
#params['mij'].set(value=0.0152,vary=0)
#params['sij'].set(value=2.22,vary=1)

#params['uij'].set(value=400,vary=1)
#params['kij'].set(value=-0.208,vary=0)
#params['mij'].set(value=0.036,vary=0)
#params['sij'].set(value=2.169,vary=1)
#params['kijS1'].set(value=kijS1,vary=1)
#params['kijS2'].set(value=kijS2,vary=1)
#params['kijS2'].set(value=-0.03,vary=0)
#params['kijS1'].set(value=0,vary=0)
#params['kijS2'].set(value=0,vary=0)
#grid_x = [np.unique(par.ravel()) for par in result.brute_grid]
#params=result.params
#grid=result.brute_grid
#Obj=result.brute_Jout
#plt.plot(grid_x[0],Obj[:,0,0,0])
#plt.plot(grid_x[1],Obj[0,:,0,0])
#plt.plot(grid_x[2],Obj[0,0,:,0])
#plt.plot(grid_x[3],Obj[0,0,0,:])

# calculate final result

#final = data + result.residual
final = data + TotalFit(params,x,VLEterndata)
#final = data + PCSAFT(params,x,data)
# write error report
#report_fit(result)
#params.pretty_print()


# try to plot results

yla="$w_w$[-]"
xla="$p_{w}/p_{0w}^{LV}$ [-]"
fig12,ax12=Thermoplot.Plot(xla,yla)
try:
    for i in range(nASD):
        ax12.plot(final[i*nData:(i+1)*nData-1],x[2][i*nData:(i+1)*nData-1], 'k-',Color=Colorcode[i+1],label="$w_{IND}^F=$"+str(round(1-ZogP[i],1)))
        ax12.plot(data[i*nData:(i+1)*nData-1],x[2][i*nData:(i+1)*nData-1],'x',Color=Colorcode[i+1])
except ImportError:
    pass
ax12.legend()

cwd = os.getcwd()
path1=cwd+"\\SLEbin.csv"
Data=pd.read_csv(path1)
Tdata=Data.values[:,0]
ww=Data.values[:,1]

#Pol="pvpva64"
wL=PCSAFTSLE(params,Tdata,API,Pol)
T_fun=InterpolatedUnivariateSpline(wL[:-1],Tdata[:-1])
w_fun=np.linspace(wL[0],wL[-1],100)
xla="$w_{IND}$[-]"
yla="$T[K]$"

fig1,ax1=Thermoplot.Plot(xla,yla)

ax1.plot(w_fun,T_fun(w_fun),'k')
ax1.plot(ww,Tdata,'kx')
fig1.savefig(join(cwd,"Plots",'SLE.jpeg'))

#paramlist=[]
#sijvec=np.linspace(2,3,3)
#uijvec=np.linspace(200,350,3)
#fig, axs = plt.subplots(3,3)
#for i,vali in enumerate(sijvec):
#    for j, valj in enumerate(uijvec):
#        params['sij'].set(value=sijvec[i],vary=0)
#        params['uij'].set(value=uijvec[j],vary=0)
#        params=VLEFIT(params,Sol,Pol)
#        owndata=[0,1,2,3,4,8]
#        otherdata=[4,5,6,7,9,10,11,12,13,14,15,16,17,18,19]
#        w2,RHSAFT=PCSAFTVLE(params,Sol,Pol)
#        cwd = os.getcwd()
#        path1=join(cwd,"VLEbin.csv")
#        Data=pd.read_csv(path1)
#
#        VLEdata=Data.values[:,0]
#        ww=Data.values[:,1]
#        yla="$w_w$[-]"
#        xla="$p_{w}/p_{0w}^{LV}$ [-]"
#
#        fig2,ax2=Thermoplot.Plot(xla,yla)
#        ax2.plot(RHSAFT,w2,'k')
#        ax2.plot(VLEdata[otherdata],ww[otherdata],'kx')
#        axs[i,j].plot(VLEdata[otherdata],ww[otherdata],'kx')
#        axs[i,j].plot(RHSAFT,w2,'k')
#        for k in range(1,7):
#            ax2.plot(VLEdata[owndata[k-1]],ww[owndata[k-1]],'k'+Symbol[k],MarkerFaceColor=Colorcode[k])
#            axs[i,j].plot(VLEdata[owndata[k-1]],ww[owndata[k-1]],'k'+Symbol[k],MarkerFaceColor=Colorcode[k])
#            
#            fig2.savefig(join(cwd,"Plots",'VLE.jpeg'))
#        fig2.savefig(join(cwd,"Plots",'VLE'+str(i)+str(j)+'.jpeg'))
#        paramlist.append(params)



#params['uij'].set(value=406,vary=1)
#params['kij'].set(value=-0.323,vary=0)
#params['mij'].set(value=0.017,vary=0)
#params['sij'].set(value=2.169,vary=1)
owndata=[0,1,2,3,4,7]
otherdata=[5,6,8,9,10,11,12,13,14,15,16,17,18]
w2,RHSAFT=PCSAFTVLE(params,Sol,Pol)
yla="$w_w$[-]"
xla="$p_{w}/p_{0w}^{LV}$ [-]"
cwd = os.getcwd()
path1=join(cwd,"VLEbin.csv")
Data=pd.read_csv(path1)

VLEdata=Data.values[:,0]
ww=Data.values[:,1]
fig2,ax2=Thermoplot.Plot(xla,yla)
ax2.plot(RHSAFT,w2,'k')
ax2.plot(VLEdata[otherdata],ww[otherdata],'kx')
for k in range(1,7):
    ax2.plot(VLEdata[owndata[k-1]],ww[owndata[k-1]],'k'+Symbol[k],MarkerFaceColor=Colorcode[k])
fig2.savefig(join(cwd,"Plots",'VLE.jpeg'))


w2,RHSAFT=PCSAFTVLE(params,Sol,API)
yla="$w_w$[-]"
xla="$p_{w}/p_{0w}^{LV}$ [-]"
fig22,ax22=Thermoplot.Plot(xla,yla)
cwd = os.getcwd()
path1=join(cwd,"VLEbinIND.csv")
Data=pd.read_csv(path1)
VLEdata=Data.values[:,0]
ww=Data.values[:,1]
ax22.plot(RHSAFT,w2,'k')
for k in range(1,7):
    ax22.plot(VLEdata[k-1],ww[k-1],'k'+Symbol[k],MarkerFaceColor=Colorcode[k])





cwd = os.getcwd()

path1=join(cwd,"PVT_H2OPVP.csv")
Data=pd.read_csv(path1)
rhodata=Data.values[:,1]
ww=Data.values[:,0]
w2,rho=PCSAFTPVT(params,Sol,Pol)

xla="$w_{w}$[-]"
yla="$\\rho$\n [kg/$m^3$]"
fig3,ax3=Thermoplot.Plot(xla,yla)
ax3.plot(w2,rho,'k')
ax3.plot(ww,rhodata,'kx')
fig3.savefig(join(cwd,"Plots",'PVT.jpeg'))





#uijsol=res_1['x']
#print(uijsol)