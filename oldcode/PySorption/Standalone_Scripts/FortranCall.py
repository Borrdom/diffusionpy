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
    T=np.ones_like(args[0])*298.15
    M1vec=np.ones_like(args[0])*25700
    p=np.ones_like(args[0])*1.013
    x=[]
    x.extend((T,M1vec,*args,p))
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
    print(report)
    #Lese output
    pathout=join(cwd,"output_file","VLE_ATPS.dat")
    matout= np.loadtxt(pathout)
    rho=matout[:,ncomp-1+1]
    gamma2=matout[:,ncomp-1+3]
    SLE=matout[:,-1]    
    return gamma2,SLE,rho


def PCSAFTSLE(params,Tdata,*args):
    """Creates Value Table gamma w for given Temperatures and iterates over the SLE-Equation to return real Solubilities"""
    wL=[]
    for i,val in enumerate(Tdata):
        w2=np.linspace(0,1,100)
        w1=1-w2
        p=1.013
        M1=25700
        SLEx=SetupPCSAFTVector(M1,val,p,w2,w1)
        SLE=PCSAFT(params,SLEx,*args)
        gamma2=SLE[0]
        va=SLE[1]    
        try:
            wL.append(np.max(InterpolatedUnivariateSpline(w2[10:],(gamma2-va)[10:]).roots()))
        except:
            wL.append(1)
    wL=np.asarray(wL)  
    return wL

def PCSAFTVLE(params,*args):
    """Creates Value Table gamma w for a given Temperature and calculates activity of component 2"""
    w2=np.linspace(0,1,100)
    T=np.ones_like(w2)*298.15
    M1vec=np.ones_like(w2)*25700
    p=np.ones_like(w2)*1.013
    w1=1-w2
    w3=np.zeros_like(w1)
    VLEx=SetupPCSAFTVector(M1,T,p,*args)
    M1=25700/1000
    M2=18.015/1000
    M3=357.88/1000
    def mass2mole(w1,w2,w3):
        x2=w2/M2/(w1/M1+w2/M2+w3/M3)
        x3=w3/M3/(w1/M1+w2/M2+w3/M3)
        return x2,x3
    x2,x3=mass2mole(w1,w2,w3)
    return w2,PCSAFT(params,VLEx,*args)[0]*x2


cwd = os.getcwd()
path1=join(cwd,"VLEbin.csv")
Data=pd.read_csv(path1)

RH=Data.values[:,0]
ww=Data.values[:,1]



    #params['kij2'].set(value=kijres,vary=0)   
#def CreateXvector(M1,T=298.15,w2,w1,w3,p=1.013):

# create a set of Parameters
params = Parameters()

#Bounds and Initials

lbui=200
ubui=450
lbkij=-0.4
ubkij=0
lbmi=0.01
ubmi=0.15
lbsi=2
ubsi=3
lbkij2=-0.4
ubkij2=0
ui0=lbui
kij0=lbkij
mi0=lbmi
si0=lbsi
kij20=lbkij2
uistep=50
kijstep=0.05
mistep=0.01
sistep=0.1
kij2step=0.05
params.add('uij', value=ui0, min=lbui,max=ubui,brute_step=uistep)
params.add('kij', value=kij0,min=lbkij,max=ubkij,brute_step=kijstep)
params.add('mij', value=mi0, min=lbmi, max=ubmi,brute_step=mistep)
params.add('sij', value=si0,min=lbsi, max=ubsi,brute_step=sistep)
#params.add('kij2', value=kij20,min=lbkij2, max=ubkij2,brute_step=kij2step)
params['uij'].set(value=205,vary=1)
params['kij'].set(value=-0.146,vary=0)
params['mij'].set(value=0.04,vary=0)
params['sij'].set(value=2.7,vary=1)

#params['uij'].set(value=400,vary=1)
#params['kij'].set(value=-0.205,vary=0)
#params['mij'].set(value=0.035,vary=0)
#params['sij'].set(value=2.159,vary=1)

#params['kij2'].set(value=2.7) 

#params['uij'].set(value=350,vary=0)
#params['kij'].set(vary=)
#params['mij'].set(value=0.07,vary=0)
#params['sij'].set(value=2,vary=0)


Pol="pvp"
API="indometacin"
Sol="water"

def VLEFIT(params):
    params['uij'].set(vary=0)
    params['kij'].set(vary=1)
    params['mij'].set(vary=1)
    params['sij'].set(vary=0)
    #params['kijS1'].set(vary=0)
    #params['kijS2'].set(vary=0) 
    Pol="pvp"
    API="indometacin"
    Sol="water"
    #cwd = os.getcwd()
    path1=join(cwd,"VLEbin.csv")
    Data=pd.read_csv(path1)
    
    RH=Data.values[:,0]
    ww=Data.values[:,1]
    w2=ww
    T=np.ones_like(w2)*298.15
    M1vec=np.ones_like(w2)*25700
    p=np.ones_like(w2)*1.013
    w1=1-w2
    w3=np.zeros_like(w1)
    VLEdata=RH
    VLEx=[]
    VLEx.extend((T,M1vec,w2,w1,p))
    M1=25700/1000
    M2=18.015/1000
    M3=357.88/1000
    def mass2mole(w1,w2,w3):
        x2=w2/M2/(w1/M1+w2/M2+w3/M3)
        x3=w3/M3/(w1/M1+w2/M2+w3/M3)
        return x2,x3
    x2,x3=mass2mole(w1,w2,w3)
    
    VLEbin=lambda params,x,data:    PCSAFT(params,x,data,Sol,Pol)[0]*x2-data
    minner = Minimizer(VLEbin, params ,fcn_args=(VLEx, VLEdata))
    result = minner.minimize(method='nelder')
    report_fit(result)
    kijres=result.params["kij"]
    mijres=result.params["mij"]
    params['uij'].set(vary=1)
    params['kij'].set(value=kijres,vary=0)
    params['mij'].set(value=mijres,vary=0)
    params['sij'].set(vary=1)
    
    return params







def TotalFit(params,x,data):
    params['uij'].set(vary=0)
    params['kij'].set(vary=1)
    params['mij'].set(vary=1)
    params['sij'].set(vary=0)

    #params['kijS1'].set(vary=0)
    #params['kijS2'].set(vary=0) 
    Pol="pvp"
    API="indometacin"
    Sol="water"
    #cwd = os.getcwd()
    path1=cwd+"\\VLEbin.csv"
    Data=pd.read_csv(path1)
    
    RH=Data.values[:,0]
    ww=Data.values[:,1]
    w2=ww
    T=np.ones_like(w2)*298.15
    M1vec=np.ones_like(w2)*25700
    p=np.ones_like(w2)*1.013
    w1=1-w2
    w3=np.zeros_like(w1)
    VLEdata=RH
    VLEx=[]
    VLEx.extend((T,M1vec,w2,w1,p))
    M1=25700/1000
    M2=18.015/1000
    M3=357.88/1000
    def mass2mole(w1,w2,w3):
        x2=w2/M2/(w1/M1+w2/M2+w3/M3)
        x3=w3/M3/(w1/M1+w2/M2+w3/M3)
        return x2,x3
    x2,x3=mass2mole(w1,w2,w3)
    #uij350
    #kij-0.15
    #mij0.07
    #sij2
    #uij205
    #kij-0.146
    #mij0.04
    #sij2

    VLEbin=lambda params,x,data:    PCSAFT(params,x,data,Sol,Pol)[0]*x2-data
    minner = Minimizer(VLEbin, params ,fcn_args=(VLEx, VLEdata))
    result = minner.minimize(method='nelder')
    report_fit(result)

    kijres=result.params["kij"]
    mijres=result.params["mij"]
    params['uij'].set(vary=1)
    params['kij'].set(value=kijres,vary=0)
    params['mij'].set(value=mijres,vary=0)
    #params['kij'].set(value=-0.146,vary=0)
    #params['mij'].set(value=0.04,vary=0)
    params['sij'].set(vary=1)
    report_fit(result)
    #params['kijS1'].set(vary=1)
    #params['kijS2'].set(vary=1) 
    #params=VLEFIT(params)
    
    x2,x3=mass2mole(x[3],x[2],x[4])
    
    return (PCSAFT(params,x,data,Sol,Pol,API)[0]*x2-data)
# do fit, here with the default least1sq algorithm


path1=cwd+"\\ASDZografi.csv"
Data=pd.read_csv(path1)
RH=Data.dropna().values[:,0]
Zog=[]
ZogP=[]
lis=[1,2,3,4,5]
nASD=len(lis)
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
p=np.ones_like(ww2)*1.013
x=[]
x.extend((T,M1vec,ww2,ww1,ww3,p))
VLEterndata=data

#SLE



#TotalFit(params,x,VLEterndata)
#PCSAFT SLE via Python

minner = Minimizer(TotalFit, params ,fcn_args=(x, VLEterndata))
result = minner.minimize(method='nelder')
#grid_x = [np.unique(par.ravel()) for par in result.brute_grid]
params=result.params
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


wL=PCSAFTSLE(params,Tdata,API,Pol)
T_fun=InterpolatedUnivariateSpline(wL,Tdata)
w_fun=np.linspace(wL[0],wL[-1],100)
xla="$w_{IND}$[-]"
yla="$T[K]$"

fig1,ax1=Thermoplot.Plot(xla,yla)

ax1.plot(w_fun,T_fun(w_fun),'k')
ax1.plot(ww,Tdata,'kx')
fig1.savefig(join(cwd,"Plots",'SLE.jpeg'))

paramlist=[]
sijvec=np.linspace(2,3,3)
uijvec=np.linspace(200,350,3)
fig, axs = plt.subplots(3,3)
for i,vali in enumerate(sijvec):
    for j, valj in enumerate(uijvec):
        params['sij'].set(value=sijvec[i],vary=0)
        params['uij'].set(value=uijvec[j],vary=0)
        params=VLEFIT(params)
        owndata=[0,1,2,3,4,8]
        otherdata=[4,5,6,7,9,10,11,12,13,14,15,16,17,18,19]
        w2,RHSAFT,ww,VLEdata=PCSAFTVLE(params)
        yla="$w_w$[-]"
        xla="$p_{w}/p_{0w}^{LV}$ [-]"

        fig2,ax2=Thermoplot.Plot(xla,yla)
        ax2.plot(RHSAFT,w2,'k')
        ax2.plot(VLEdata[otherdata],ww[otherdata],'kx')
        axs[i,j].plot(VLEdata[otherdata],ww[otherdata],'kx')
        axs[i,j].plot(RHSAFT,w2,'k')
        for k in range(1,7):
            ax2.plot(VLEdata[owndata[k-1]],ww[owndata[k-1]],'k'+Symbol[k],MarkerFaceColor=Colorcode[k])
            axs[i,j].plot(VLEdata[owndata[k-1]],ww[owndata[k-1]],'k'+Symbol[k],MarkerFaceColor=Colorcode[k])
            
            fig2.savefig(join(cwd,"Plots",'VLE.jpeg'))
        fig2.savefig(join(cwd,"Plots",'VLE'+str(i)+str(j)+'.jpeg'))
        paramlist.append(params)



params['uij'].set(value=406,vary=1)
params['kij'].set(value=-0.208,vary=0)
params['mij'].set(value=0.036,vary=0)
params['sij'].set(value=2.169,vary=1)
owndata=[0,1,2,3,4,8]
otherdata=[4,5,6,7,9,10,11,12,13,14,15,16,17,18,19]
w2,RHSAFT,ww,VLEdata=PCSAFTVLE(params)
yla="$w_w$[-]"
xla="$p_{w}/p_{0w}^{LV}$ [-]"
#
fig2,ax2=Thermoplot.Plot(xla,yla)
ax2.plot(RHSAFT,w2,'k')
ax2.plot(VLEdata[otherdata],ww[otherdata],'kx')
for k in range(1,7):
    ax2.plot(VLEdata[owndata[k-1]],ww[owndata[k-1]],'k'+Symbol[k],MarkerFaceColor=Colorcode[k])
fig2.savefig(join(cwd,"Plots",'VLE.jpeg'))



cwd = os.getcwd()

path1=cwd+"\\PVT_H2OPVP.csv"
Data=pd.read_csv(path1)
rhodata=Data.values[:,1]
ww=Data.values[:,0]

w2=np.linspace(0,1,100)
T=np.ones_like(w2)*298.15
M1vec=np.ones_like(w2)*25700
p=np.ones_like(w2)*1.013
w1=1-w2
w3=np.zeros_like(w1)
VLEdata=RH
x=[]
x.extend((T,M1vec,w2,w1,p))
PCSAFT(params,x,RH,Sol,Pol)


xla="$w_{w}$[-]"
yla="$\\rho$\n [kg/$m^3$]"
fig3,ax3=Thermoplot.Plot(xla,yla)
ax3.plot(w2,rho,'k')
ax3.plot(ww,rhodata,'kx')
fig3.savefig(join(cwd,"Plots",'PVT.jpeg'))





#uijsol=res_1['x']
#print(uijsol)