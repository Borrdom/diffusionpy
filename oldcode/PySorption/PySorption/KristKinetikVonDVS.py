from scipy.interpolate import InterpolatedUnivariateSpline,UnivariateSpline
import pandas as pd

from os.path import join
from os import getcwd
import tkinter
from tkinter import *
from tkinter.filedialog import askopenfilename
import matplotlib.pyplot as plt
import numpy as np
import ntpath
from lmfit import Parameters, minimize
from lmfit.printfuncs import report_fit
def SorptionAndKris(dvspath,RH=0.9,DL=0.7,T=25,n=1):
    cwd=getcwd()
    molecule1="pvpva64"
    molecule2="naproxen"
    rho1=1180
    rho2=1250
    
    RHstr=str(max(int(RH*100),RH*100))+"%"
    Tstr=str(max(int(T),T))+"°C"
    DLstr=str(max(int(DL*100),DL*100))+"%"
    wPolyASD=1-DL
    wAPIASD=DL
    rhoASD=(wPolyASD/rho1+wAPIASD/rho2)**-1

    Outputpath=join(cwd,"KristundSorptionskinetik",RHstr+"_"+Tstr,"DL"+DLstr)
    simpath=join(cwd,"FeedWasser",RHstr+"_"+Tstr,"_".join(["EXP",RHstr,molecule1,molecule2,Tstr]))

    INQ=pd.read_excel(dvspath)
    TSorption=INQ["t[min]"].dropna().values.astype(float)
    MSorption=INQ["w[-]"].dropna().values.astype(float)
    TKrist=INQ["tK[min]"].dropna().values.astype(float)[::50]
    MKrist=INQ["wK[-]"].dropna().values.astype(float)[::50]
    M0=INQ["M0[mg]"].dropna().values.astype(float)
    D0=0.0145
    A0=np.pi/4*D0**2
    L=M0/rhoASD/A0

    #####
    
    SAFT=pd.read_excel(simpath+".xlsx")
    #EXP=pd.read_excel(exppath+".xlsx")
    wAPISIMraw=SAFT["wAPI[-]"].values
    wwSIMraw=SAFT["ww[-]"].values
    wwSLE=SAFT["wwSLE[-]"].dropna().values
    wAPISLE=SAFT["wAPISLE[-]"].dropna().values

    #wAPISIM=np.hstack((wAPISIMraw[wAPISIMraw<wAPIASD-0.09],wAPIASD)) #Precision problem 0.09
    #wwSIM=np.hstack((wwSIMraw[wAPISIMraw<wAPIASD-0.09],MKrist[0])) #A great problem arises when the water mass surpasses the crystalline ASD
    wAPISIM=wAPISIMraw[wAPISIMraw<=wAPIASD] #Precision problem 0.09
    wwSIM=wwSIMraw[wAPISIMraw<=wAPIASD]
    XwSIM=wwSIM/(1-wwSIM)
    
    def MakeItBijective(ww,wAPI):
        hasmaximum=np.argmax(ww) if wAPI[np.argmax(ww)]!=0 else None
        hasminimum=np.argmin(ww) if wAPI[np.argmin(ww)]!=1 else None
        hasmaximum=None
        hasminimum=None
        if hasmaximum!=None and hasminimum!=None:
            raise NameError("Crazy it has both maximum and minimum. Better abort")
    
        elif hasmaximum!=None:
            wwsortleft,wAPIsortleft=np.sort(ww[:hasmaximum+1]),wAPI[:hasmaximum+1][np.argsort(ww[:hasmaximum+1])]
            wwsortright,wAPIsortright=np.sort(ww[1+hasmaximum:]),wAPI[hasmaximum+1:][np.argsort(ww[1+hasmaximum:])]
            wAPIleft_fun=InterpolatedUnivariateSpline(wwsortleft,wAPIsortleft,k=1,ext=3)
            wAPIright_fun=InterpolatedUnivariateSpline(wwsortright,wAPIsortright,k=1,ext=3)
            wamax=wAPI[hasmaximum]
            output=lambda x,y: wAPIleft_fun(x) if y<wamax else wAPIright_fun(x)
        elif hasminimum!=None:
            wwsortleft,wAPIsortleft=np.sort(ww[:hasminimum+1]),wAPI[:hasminimum+1][np.argsort(ww[:hasminimum+1])]
            wwsortright,wAPIsortright=np.sort(ww[1+hasminimum:]),wAPI[hasminimum+1:][np.argsort(ww[1+hasminimum:])]
            wAPIleft_fun=InterpolatedUnivariateSpline(wwsortleft,wAPIsortleft,k=1,ext=3)
            wAPIright_fun=InterpolatedUnivariateSpline(wwsortright,wAPIsortright,k=1,ext=3)
            wamin=wAPI[hasminimum]
            output=lambda x,y: wAPIleft_fun(x) if y<wamin else wAPIright_fun(x)
        else:
            wwsort,wAPIsort=np.sort(ww),wAPI[np.argsort(ww)]
            output=lambda x,y=1:InterpolatedUnivariateSpline(wwsort,wAPIsort,k=2,ext=0)(x) #the higher the order the more problematic
        return output
    wtest=np.linspace(0,1,100)
    wAPISIM_fun=MakeItBijective(wwSIM,wAPISIM)
    wwtest=np.linspace(0,np.max(wwSIM),100)
    wwSIM_fun=InterpolatedUnivariateSpline(np.sort(wAPISIM_fun(wwtest)),wwtest[np.argsort(wAPISIM_fun(wwtest))])
    wwPhiL0=wwSIM_fun(wAPIASD)
    wwEXP=MKrist-(MKrist[0]-wwPhiL0) #shift the water fraction of the experiment so that the predicted water amount for that feedconcentration is zero
    XwEXP=wwEXP/(1-wwEXP)
    PhiLeq=((1-wAPISIM_fun(wwSLE))/wPolyASD*(1-wwSLE)+wwSLE)**-1
    PhiSeq=1-PhiLeq #calculate the equilibrium solid fraction#calculate the zero krist point of the current sorption isotherm
     # when this is smaller than zero the experiment wont start at phi=0 but the reverse is everytime 0 so better less than more
    wwMax=wwSLE*PhiLeq
    XwSIM_fun=InterpolatedUnivariateSpline(wAPISIM,XwSIM)
    fit_params = Parameters() 
    fit_params.add('phiL', value=1, min=0, max=1,vary=1)
    fit_params.add('Alpha', value=0, min=0, max=1,vary=1)

    def PHIL(params,wwEX,wPolyASD):
        phiL=params["phiL"].value
        phiLL=wwEX+wPolyASD/(1-wAPISIM_fun(wwEX*1/phiL))*(1-wwEX)#safe
        return np.nan_to_num(abs(phiL-phiLL))
    def ALPHA(params,XwEX,DL):
        Alpha=params["Alpha"].value
        XaEX=(DL-Alpha)/(1-Alpha)
        AAlpha=1-XwEX/XwSIM_fun(XaEX)
        print(Alpha)
        return np.nan_to_num(abs(Alpha-AAlpha))
    alpha = np.asarray([minimize(ALPHA,args=(XwEXP[i],DL), params=fit_params,method="least-squares",tol=1E-6).params["Alpha"] for i,val in enumerate(XwEXP)])
    PhiS=alpha*(1-wwEXP)
    PhiL=1-PhiS
    # PhiL = np.asarray([minimize(PHIL,args=(wwEXP[i],wPolyASD), params=fit_params,method="least-squares",tol=1E-6).params["phiL"] for i,val in enumerate(wwEXP)])
    # PhiS=1-PhiL
    # alpha=PhiS/PhiSeq
    # alpha=PhiS/(1-wwEXP)
    fig11,ax11=plt.subplots()
    ax11.plot(wAPISIMraw,wwSIMraw,'co',label="Isotherm_points")
    ax11.plot(wAPISIM,wwSIM,'go',label="Isotherm_points_used")
    ax11.plot(wtest,wwSIM_fun(wtest),'g-',label="Isotherm")
    ax11.plot(wAPIASD,MKrist[0],'kx',label="Starting_water_fraction")
    ax11.plot(wAPIASD,wwEXP[0],'ro',label="Shifted_water_fraction")
    ax11.plot(wAPISIM_fun(wwEXP*1/PhiL),wwEXP*1/PhiL,'r-',label="Shifted_Trajectory")
    ax11.plot(wAPISIM_fun((wwEXP+(MKrist[0]-wwPhiL0))*1/PhiL),(wwEXP+(MKrist[0]-wwPhiL0))*1/PhiL+MKrist[0]-wwPhiL0,'k-',label="Trajectory")
    ax11.plot([wAPISLE,1],[wwSLE,0],'mo-',label="Crystalline_ASD")
    ax11.set_ylabel("$w_w$[-]")
    ax11.set_xlabel("$w_{API}^F$[-]")
    ax11.legend()
    fig13,ax13=plt.subplots()
    ax13.plot(TKrist,wAPISIM_fun(wwEXP*1/PhiL),'kx')
    ax13.set_ylabel("$w_{API}^F$[-]")
    ax13.set_xlabel("t[min]")

    
    def FitKrist(t,alpha):
        fit_params = Parameters() 
        
        fit_params.add('n', value=n, min=1, max=8,vary=0)
        fit_params.add('k', value=-4*n, min=-50, max=0,vary=1)
        fit_params.add('tau', value=0, min=0, max=100000,vary=1)
        def JMAK(params,t2,alpha2):
            k=10**np.asarray(params["k"].value)
            n=int(np.asarray(params["n"].value))
            tau=np.asarray(params["tau"].value)
            t2=t2*60
            alpha=1-np.exp(-k*(t2-tau)**n)
            #print(alpha-alpha2)
            #print(tau)
            return alpha-alpha2
        out = minimize(JMAK,args=(t[alpha>0.15]-t[0],alpha[alpha>0.15]), params=fit_params,method="least-squares")
        kopt=10**out.params['k'].value
        nopt=out.params['n'].value
        tauopt=out.params['tau'].value/60
        #print(kopt)
        #print(tauopt)
        #print(out.params['n'].value)
        t3=np.linspace(t[0],t[-1]*10,1000)
        opt=JMAK(out.params,t3-t[0],np.zeros_like(t3))    
    
        fig,ax=plt.subplots()
        ax.plot(t-t[0]-tauopt,alpha,'kx')
        ax.plot(t3-t[0]-tauopt,opt,'k-')
        ax.set_xlabel("t [min]")
        ax.set_ylabel("alpha[-]")
        ax.set_xlim([-tauopt,6000])
        ax.set_ylim([0,1])
        twrite=t-t[0]-tauopt
        tsimwrite=t3-t[0]-tauopt
        tstart=t
        tseperate=t-t[0]
        datawrite=pd.concat((pd.DataFrame(tstart),pd.DataFrame(tseperate),pd.DataFrame(twrite),pd.DataFrame(alpha),pd.DataFrame(wAPISIM_fun(wwEXP*1/PhiL)),pd.DataFrame(tsimwrite),pd.DataFrame(opt),pd.DataFrame(np.asarray([kopt])),pd.DataFrame(np.asarray([nopt])),pd.DataFrame(np.asarray([tauopt]))),ignore_index=True, axis=1)
        datawrite.columns=["tStart[min]","tSeperate[min]","tKrist[min]","alphaExp[-]","wAPIFEXP[-]","tSim[min]","alphaSim[-]","kopt[1/s]","n[-]","tauopt[min]"]
        
        filename=ntpath.basename(dvspath).replace(".xlsx","Out_Krist.xlsx")
        writename=join(Outputpath,filename)
        datawrite.to_excel(writename)
        
    def FitSorption(TSorption,MSorption,L):
        fit_params = Parameters() 
        fit_params.add('D', value=1, min=.001, max=10000,vary=1)
        fit_params.add('tau', value=1, min=0, max=100000,vary=1)
        t=TSorption
        deltat=np.round(np.average(t[1:]-t[:-1]),2)
        #print(deltat)
        t=t-t[0]
        m=MSorption
        mp=(m-m[0])/(m[-1]-m[0])
        def Cranc1(params,t2,mp2):
            n=12
            
            LCranc=2*L*1E-6
            D=np.asarray(params["D"].value)*1E-14
            tau=np.asarray(params["tau"].value)
            Forier=0
            t2=t2*60
            for cc in range(n):
                Forier=Forier+1/(2*cc+1)**2*np.exp(-((2*cc+1)**2*np.pi**2*D*(t2-tau))/(LCranc**2))
                mperc=1-8/np.pi**2*Forier
            return np.fmin(np.fmax(mperc,0),1)-mp2
        
        out = minimize(Cranc1,args=(t[mp>0],mp[mp>0]), params=fit_params,method="powell")
        Dopt=out.params['D'].value*1E-14
        tauopt=out.params['tau'].value/60
        print(Dopt)
        t3=np.linspace(t[0],t[-1],1000)
        opt=Cranc1(out.params,t3,np.zeros_like(t3))    
        mplotexp=mp*(-m[0]+m[-1])+m[0]
        mplot=opt*(-m[0]+m[-1])+m[0]
        t99=t3[next(i for i,v in enumerate(opt-0.99>0) if v)]
        t98=t3[next(i for i,v in enumerate(opt-0.98>0) if v)]
        t95=t3[next(i for i,v in enumerate(opt-0.95>0) if v)]
        fig,ax=plt.subplots()
        ax.plot(t-tauopt,mplotexp,'kx')
        ax.plot(t3-tauopt,mplot,'k-')
        ax.set_xlabel("t [min]")
        ax.set_ylabel("m [mg]")
        ax.plot([t99,t99],[m[0],mplot[t3==t99][0]],'b--')
        ax.plot([0,t99],[mplot[t3==t99][0],mplot[t3==t99][0]],'b--',label='99% of GGW')
        ax.plot([t98,t98],[m[0],mplot[t3==t98][0]],'g--')
        ax.plot([0,t98],[mplot[t3==t98][0],mplot[t3==t98][0]],'g--',label='98% of GGW')
        ax.plot([t95,t95],[m[0],mplot[t3==t95][0]],'r--')
        ax.plot([0,t95],[mplot[t3==t95][0],mplot[t3==t95][0]],'r--',label='95% of GGW')
        ax.legend()
        twrite=t-tauopt
        tsimwrite=t3-tauopt
        datawrite=pd.concat((pd.DataFrame(twrite),pd.DataFrame(mplotexp),pd.DataFrame(tsimwrite),pd.DataFrame(mplot),pd.DataFrame(np.asarray([Dopt])),pd.DataFrame(np.asarray([t99])),pd.DataFrame(np.asarray([t98])),pd.DataFrame(np.asarray([t95])),pd.DataFrame(np.asarray([tauopt]))), ignore_index=True, axis=1)
        datawrite.columns=["t[min]","wExp[-]","tSim[min]","wSim[-]","Dopt[m^2/s]","t99[min]","t98[min]","t95[min]","tauopt[min]"]
        filename=ntpath.basename(dvspath).replace(".xlsx","Out_Diffusion.xlsx")
        writename=join(Outputpath,filename)
        datawrite.to_excel(writename)
    FitKrist(TKrist,alpha)
    FitSorption(TSorption,MSorption,L)




class window():
    def __init__(self):
        self.main()
        self.path=askopenfilename()
        self.n=1          
    
    def EnterBtn(self):
        RH=float(self.Box0.get())
        DL=float(self.Box1.get())
        T=float(self.Box2.get())
        n=float(self.Box3.get())
        #plt.close("all")
        SorptionAndKris(self.path,RH=RH,DL=DL,T=T,n=n)
        
        
    def searchdirectory(self):
        self.root.withdraw()
        self.path= askopenfilename()
        self.root.deiconify()
        self.label2.config(text='Current Path is '+ self.path)
    
    def main(self):
     """ Main Window for GUI"""
     self.root = Tk()
     self.root.attributes("-topmost", True)

     #root.geometry("350x150+450+280")  # "width x height + position right + position down"
     self.root.title("Predict Diffusion")
     #root.iconbitmap("TU.ico")#here the TU Dortmund Icon has be in the same folder as the python script and 
                              #the image as to be an icon type.
     var0 = StringVar(self.root, value='0.9')
     self.Box0=Entry(self.root, textvariable=var0)
     self.Box0.pack(side="top", fill="x")
     Labeltext0 = StringVar()
     Labeltext0.set("RH[-]")
     labelDir0=Label(self.root, textvariable=Labeltext0, height=4)
     labelDir0.pack(side="left")
     
     var1 = StringVar(self.root, value='0.7')
     self.Box1=Entry(self.root, textvariable=var1)
     self.Box1.pack(side="top", fill="x")
     Labeltext1 = StringVar()
     Labeltext1.set(r'DL[-]')
     labelDir1=Label(self.root, textvariable=Labeltext1, height=4)
     labelDir1.pack(side="left")
     
     var2 = StringVar(self.root, value='25')
     self.Box2=Entry(self.root, textvariable=var2)
     self.Box2.pack(side="top", fill="x")
     Labeltext2 = StringVar()
     Labeltext2.set(r'T[°C]')
     labelDir2=Label(self.root, textvariable=Labeltext2, height=4)
     labelDir2.pack(side="left")
     
     var3 = StringVar(self.root, value='1')
     self.Box3=Entry(self.root, textvariable=var3)
     self.Box3.pack(side="top", fill="x")
     Labeltext3 = StringVar()
     Labeltext3.set(r'n[-]')
     labelDir3=Label(self.root, textvariable=Labeltext3, height=4)
     labelDir3.pack(side="left")
     
     self.Btn0 = Button(self.root, text="enter", padx=50, pady=2, font='Helvetica 8 italic', command=self.EnterBtn)
     self.Btn0.pack(side="left")
     self.Btn1 = Button(self.root, text="browse", padx=50, pady=2, font='Helvetica 8 italic', command=self.searchdirectory)
     self.Btn1.pack(side="left")
     self.Btn2 = Button(self.root, text="close_plots", padx=50, pady=2, font='Helvetica 8 italic', command=lambda:plt.close("all"))
     self.Btn2.pack(side="left")
     
     self.label2=Label(self.root,text='Current Path is not chosen click browse',bg='gold',fg='blue')
     self.label2.pack(side="top")
     
     self.root.mainloop()  
if __name__=="__main__":
    #EndpointOptimization()
    window1=window()
    
#plt.plot(wwSIM,wAPISIM_fun(wwEXP,1))
#plt.plot(wwEXP_fun.derivative(n=1)(wAPIEXP))
#plt.plot(wwSIM_fun.derivative(n=1)(wAPISIM))
#pd.read_xlsx("")