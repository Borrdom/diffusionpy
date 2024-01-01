from PySorption.DiffusionFinal import Film
from PySorption.dcc import window
import os
import tkinter as tk
from functools import partial
import matplotlib.pyplot as plt
import ast
import numpy as np
from PyCSAFT.list_manipulations import averagelisttimeseries
#import LeastFitDiffusionFunktionen as LDF
from PySorption.LeastFitDiffusion import FitDiffusionData
from PyCSAFT.get_pcpar import get_par
from  PySorption.DiffusionsmodellPlots import Plotlist
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

#from LeastFitDiffusion import FitDiffusionData
#je höher konz dest niederiger muss dehnung sein
#Größe der Federkonstanten 1E7/1E7 schneller im Relax GGW vs 1E6/1E6 langsamer im relax GGW
#sigma/ rho2 steigt linear an und je nach relax zeit sinkt es auf GGWab
#E1 den GGW wert steuert
#E0 steuert die steigung von sigma/rho2
File=True
class empty:
    result=""
    components=({"Mi":65,"rho0": 1190,"name":"pvpva64"},{"Mi":0.0018015,"rho0": 997,"name":"water"},{"Mi":65,"rho0": 1190,"name":"pvpva64"})
    wPolyASD=1
    pass
originalpath=os.getcwd()
window1=window() if File else empty()
Data=window1.result
Film1=Film(window1.result,window1.components,window1.wPolyASD)# if len(window1.components)==2 else ASDFilm(window1.result,*window1.components)
Film1.File=File if not File else Film1.File
Film1.exp_boundary()

class Checkbar(tk.Frame):
   def __init__(self, parent=None, picks=[],commands= lambda : None, side="top", anchor="w"):
      tk.Frame.__init__(self, parent)
      self.vars = []
      for pick,command in zip(picks,commands):
         var = tk.IntVar()
         chk = tk.Checkbutton(self, text=pick, variable=var,command=command)
         chk.pack(side=side, anchor=anchor, expand="yes")
         self.vars.append(var)
   def state(self):
      return map((lambda var: var.get()), self.vars)


class window3:
    def __init__(self,Film1,Data):
        self.Film=Film1
        self.Data=Data
        self.main()

    def enter(self):
        self.Film.indJ=int(self.Box0.get())
    def indexjump(self,indJ):
        indJ=int(self.spinbox.get())-1 if self.spinbox.get()!='' else 0
        self.Film.indJ=indJ
        self.label3.config(text="Current Jump is "+str(indJ))
        self.PlotRHSteps()
        plt.close('all')
    def PlotRHSteps(self):
        matplotlib.use("TkAgg")
        x=np.arange(0,self.Film.nHum+1,1)
        y=np.hstack((0,self.Film.Feuchtereal)) 
        self.line1.set_data(x,y)
        self.line2.set_data(x[self.Film.indJ:self.Film.indJ+2],y[self.Film.indJ:self.Film.indJ+2])
        [val.remove() for i,val in enumerate(self.line3)]
        self.line3=[]
        for index in range(len(x)-1):
            self.line3.append(self.ax.text(x[index]+0.1, y[index]-0.02, x[index+1], size=7))
        self.canvas.draw()

    def Desorption(self):
        self.Film.des=True if self.Film.des==False else False
        self.Film.exp_boundary()
        #[val.destroy() for i,val in enumerate(self.BtnX)]
        #self.BtnX =[tk.Button(self.frame, text="Jump "+str(i), padx=0, pady=2, font='Helvetica 7 bold', command=partial(self.indexjump, i)) for i in range(self.Film.nHum)]
        #[val.pack(side="top", fill="x") for i,val in enumerate(self.BtnX)]
        self.spinbox.destroy()
        self.spinbox = tk.Spinbox(self.frame, from_ = 1, to = self.Film.nHum,increment = 1,command=partial(self.indexjump, 1))
        self.spinbox.pack(side="top", fill="x")
        self.label4.config(text="Current Mode: \tDesorption") if self.Film.des==True else self.label4.config(text="Mode is Absorption")
        self.Film.PCSAFT(self.Film.p,self.Film.T)
        self.PlotRHSteps()
    def Ideal(self):
        self.Film.ideal=1 if self.Film.ideal==0 else 0
        strideal="Not-Ideal" if self.Film.ideal==0 else "Ideal"
        self.label5.config(text="Current Mode: \t"+strideal)
    def NewFilm(self):

        self.Film.solvent=self.BoxB.get()
        pure,kij=get_par([self.Film.molecule,self.Film.solvent,self.Film.molecule2]) if hasattr(self,"molecule2") else get_par([self.Film.molecule,self.Film.solvent,self.Film.molecule])
        self.Film.__init__(self.Film.DATA_FILE,pure,wPolyASD=self.Film.wPolyASD)
        self.Film.A0=float(self.BoxA.get())
        self.Film.Durchmesser=(self.Film.A0*4/np.pi)**(1/2)
        self.Film.exp_boundary()
        self.Film.PCSAFT(self.Film.p,self.Film.T)
    def Origin(self):
        self.Film.Origin=1 if self.Film.Origin==0 else 0
        strorigin="Origin Off" if self.Film.Origin==0 else "Origin On"
        self.label8.config(text="Current Mode: \t"+strorigin)

    def ActivityFromIsotherm(self):
        self.Film.SAFT=True if self.Film.SAFT==False else False
        strSAFT="PC-SAFT" if self.Film.SAFT==True else "comes from the Isotherm"
        self.label6.config(text="Current Activity: \t"+strSAFT)
        self.Film.PCSAFT(self.Film.p,self.Film.T)
    def NET(self):
        self.Film.NET=True if self.Film.NET==False else False
        strNET="On" if self.Film.NET==True else "Off"
        self.label7.config(text="NET-GP is "+strNET)
        self.Film.PCSAFT(self.Film.p,self.Film.T)
    def Actual(self):
        self.Film.Actual=True if self.Film.Actual==False else False
        strActual="On" if self.Film.Actual==True else "Off"
        self.label9.config(text="Actual RH is "+strActual)
        self.Film.exp_boundary()
        self.Film.PCSAFT(self.Film.p,self.Film.T)
    def SetConditions(self):
        self.Film.File=False
        self.Film.RHd=ast.literal_eval(self.BoxRHd.get())
        self.Film.wd=ast.literal_eval(self.Boxwd.get())
        self.Film.L0d=float(self.BoxL0d.get())
        self.Film.td=ast.literal_eval(self.Boxtd.get())
        self.Film.exp_boundary()
        self.Film.PCSAFT(self.Film.p,self.Film.T)

    def Diffusion1D(self):
        self.Film.Mode="NEAT"
        
        self.Film.D12=float(self.Box0.get())
        self.Film.Diffusion1D()
        #from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        #canvas = FigureCanvasTkAgg(self.Film.FigSorption, master=self.root3)
        #canvas.get_tk_widget().pack()
        #canvas.draw()
        #plt.close("all")
    def Dissolution1D(self):
        self.Film.Mode="Diss"
        self.Film.D12=float(self.Box0.get())
        self.Film.Diffusion1D()
    def Diffusion1DMEOS(self):
        #self.Film.Diffusion1DMEOS(D12=self.Film.D12,E0=self.Film.E0,tau0=self.Film.tau0,nJ=self.Film.nJ)
        #self.Film.Diffusion1DMEOS()
        self.Film.Mode="MEOS"
        self.Film.Diffusion1D()
        #self.Film.SortpionPlot()
    def RecalcIsotherm(self):
        self.Film.Recalculate_Isotherm(self.Film.p, self.Film.T)
        self.Film.PCSAFT(self.Film.p,self.Film.T)
    def Least(self):
        FitDiffusionData(self.Film,self.Data,Fit=True)

    def Prediction(self):
        FitDiffusionData(self.Film,self.Data,Fit=False)
    def Diffusion1DMEOSAll(self):
        self.Film.Mode="MEOS"
        self.Film.Diffusion1D()
        #self.Film.Diffusion1DMEOS()


        tveclist,wtreflist,wtlist,wtsimlist,tsimlist,wtficklist=[],[],[],[],[],[]
        for i,val in enumerate(Data.values()):
            self.Film.DATA_FILE=os.path.join("Pure",self.Film.molecule,val)
            self.Film.exp_boundary()
            print("Jump Number is "+str(self.Film.indJ))
            #self.Film.Diffusion1DMEOS()

            self.Film.Diffusion1D()
            tveclist.append([self.Film.tlist[self.Film.indJ]])
            tsimlist.append([self.Film.tsim])
            wtreflist.append([self.Film.wtref])
            wtlist.append([self.Film.wt])
            wtsimlist.append([self.Film.wtsim])
            #psi2IItlist.append([self.Film.psi2IIt])
            #psitlist.append([self.Film.psit])
        #self.Film.TriebkraftPlot()
        #self.Film.ContourSpatial()
        #self.Film.Mode="NEAT"
        #self.Film.ideal=0
        for i,val in enumerate(Data.values()):
            self.Film.DATA_FILE=os.path.join("Pure",self.Film.molecule,val)

            self.Film.exp_boundary()
            self.Film.E0=1E6
            self.Film.tau0=1
            self.Film.nJ=1
            #self.Film.Diffusion1DMEOS()
            self.Film.Diffusion1D()
            wtficklist.append([self.Film.wtsim])

        wtave, wtstd=averagelisttimeseries(tveclist,wtlist)
        wtsimave, wtsimstd=averagelisttimeseries(tsimlist,wtsimlist)
        wtrefave, wtrefstd=averagelisttimeseries(tveclist,wtreflist)
        wtfickave, wtfickstd=averagelisttimeseries(tsimlist,wtficklist)
        #psitave, psitstd=LDF.averagelisttimeseries(tsimlist,psitlist)
        #psi2IItave, psi2IItstd=LDF.averagelisttimeseries(tsimlist,psi2IItlist)

        tarr=np.array(tveclist)
        tarrsim=np.array(tsimlist)
        nt=tarr.shape[2]
        for i in range(nt):
            tarr[:,:,i]=tarr[:,:,i]-tarr[:,:,0]
            tarrsim[:,:,i]=tarrsim[:,:,i]-tarrsim[:,:,0]
        tendn=tarr[:,:,-1]
        tendnsim=tarrsim[:,:,-1]
        amin=np.argmin(tendn,axis=0)
        aminsim=np.argmin(tendnsim,axis=0)
        tveca=np.asarray([tarr[val,i,:] for i,val in enumerate(amin)])
        tvecasim=np.asarray([tarrsim[val,i,:] for i,val in enumerate(aminsim)])

        #fig,ax=plt.subplots()

        #ax.plot(np.sqrt(tveca[0,:]/60),wtave[0,:],'k-')
        #ax.plot(np.sqrt(tveca[0,:]/60),wtrefave[0,:],'kx')
        #[ax.plot(np.sqrt(tveclist[i][0]/60),wtlist[i][0],'r-') for i,val in enumerate(wtlist)]
        #[ax.plot(np.sqrt(tveclist[i][0]/60),wtreflist[i][0],'ro') for i,val in enumerate(wtreflist)]
        #ax.errorbar(x=np.sqrt(tveca[0,:]/60),y=wtave[0,:],yerr=wtstd[0,:],fmt='k-')
        #ax.errorbar(x=np.sqrt(tveca[0,:]/60),y=wtrefave[0,:],yerr=wtrefstd[0,:],fmt='kx')
        #DiffusionsmodellPlots.WriteHistorycsv(tveca[0,:],wtave[0,:],wtrefave[0,:],wtstd[0,:],wtrefstd[0,:])

        #DiffusionsmodellPlots.Plotlist([tveca[0,:],wtrefave[0,:],wtrefstd[0,:],tvecasim[0,:],wtsimave[0,:],tvecasim[0,:],wtfickave[0,:]])
        print("ALL")
        Plotlist([tveca[0,:]/60,tveca[0,:]/60],[wtrefave[0,:],wtrefave[0,:]],[tvecasim[0,:]/60,tvecasim[0,:]/60],[wtsimave[0,:],wtfickave[0,:]],ye=[wtrefstd[0,:],wtrefstd[0,:]],Origin=True,xla="t/min",yla="ww/-",filename="test")
        print("ALL")
        #DiffusionsmodellPlots.TriebkraftPlot(tveca[0,:],psi2IItave[0,:],psitave[0,:],name="Triebkraft")
        #DiffusionsmodellPlots.SingleSorptionPlot(tveca[0,:],wtrefave[0,:],wtrefstd[0,:],tveca[0,:],wtave[0,:],tveca[0,:],wtave[0,:])

    def LeastMEOS(self):
        #self.Film.Diffusion1DMEOS(D12=self.Film.D12,E0=self.Film.E0,E1=self.Film.E1,eta0=self.Film.eta0)
        D12True=input("Vary D12? <<y>> or <<n>>\n")=="y"
        exTrue=input("Vary exponent? <<y>> or <<n>>\n")=="y"
        ETrue=input("Vary E? <<y>> or <<n>>\n")=="y"
        tauTrue=input("Vary tau? <<y>> or <<n>>\n")=="y"
        TswitchTrue=input("Vary Tswitch? <<y>> or <<n>>\n")=="y"
        C3True=input("Vary C3? <<y>> or <<n>>\n")=="y"
        Jumps=np.asarray(ast.literal_eval(self.BoxFit.get()))
        self.Film.LeastMEOS(D12True,exTrue,ETrue,tauTrue,TswitchTrue,C3True,Jumps=Jumps)

#        self.Film.SortpionPlot()
    def LeastMEOSGES(self):
        #self.Film.Diffusion1DMEOS(D12=self.Film.D12,E0=self.Film.E0,E1=self.Film.E1,eta0=self.Film.eta0)
        D12True=input("Vary D12? <<y>> or <<n>>\n")=="y"
        exTrue=input("Vary exponent? <<y>> or <<n>>\n")=="y"
        ETrue=input("Vary E? <<y>> or <<n>>\n")=="y"
        tauTrue=input("Vary tau? <<y>> or <<n>>\n")=="y"
        C3True=input("Vary C3? <<y>> or <<n>>\n")=="y"
        TswitchTrue=input("Vary Tswitch? <<y>> or <<n>>\n")=="y"
        Jumps=np.asarray(ast.literal_eval(self.BoxFit.get()))
        for i in range(self.Film.nHum):
            self.Film.indJ=i
            self.Film.LeastMEOS(D12True,exTrue,ETrue,tauTrue,TswitchTrue,C3True,Jumps=Jumps)
        #self.Film.SortpionPlot()

    def enterMEOS(self):
        self.Film.D12=float(self.Box0.get())

        self.Film.E0=np.asarray(ast.literal_eval(self.Box1.get()))
        self.Film.exponent=np.asarray(ast.literal_eval(self.Box2.get()))
        self.Film.tau0=np.asarray(ast.literal_eval(self.Box3.get()))
        self.Film.nJ=int(self.Box5.get())
        self.Film.Xhi=float(self.Box6.get())
        self.Film.Tswitch=float(self.BoxTswitch.get())
        self.Film.C3=float(self.BoxC3.get())


        Error=False
        if np.size(self.Film.E0)<self.Film.nJ or np.size(self.Film.tau0)<self.Film.nJ:
            Error=True
        #self.Diffusion1DMEOSAll() if not Error else None #self.Diffusion1DMEOSAll() if more than one 04.09.20
        if not Error:
            self.Diffusion1DMEOSAll() if isinstance(Data,dict) else self.Diffusion1DMEOS()

            #plt.pause(1)
            #pyautogui.PAUSE =20
            #myScreenshot = pyautogui.screenshot()
            #path="C:\MEOSEXPERIMENTS"

            #stringoptions="D12"+self.Box0.get()+"E0"+self.Box1.get()+"tau0"+self.Box3.get()+"nJump"+str(self.Film.indJ)
            #myScreenshot.save(path+"\\"+stringoptions+".png")
        else:
            print("Size of E or tau is smaller than the number of elements")
        
        Tgspan=self.Film.T-self.Film.Tgspan
        Tgreal=self.Film.T-self.Film.Tgreal
        WLF=self.Film.WLF
        AR=self.Film.AR
        WLFreal=self.Film.WLFreal
        ARreal=self.Film.ARreal
        CO=np.fmin(WLF,AR)
        COreal=np.fmin(WLFreal,ARreal)
        #Plotlist([Tgspan,Tgspan,Tgspan],[WLF,AR,CO],[Tgspan,Tgspan,Tgspan],[WLF,AR,CO],Origin=self.Film.Origin,xla="T-Tg/K",yla="eta/Pas",filename="ViscoPlot")
        Plotlist([Tgspan,Tgspan,Tgspan,Tgreal,Tgreal,Tgreal],[WLF,AR,CO,WLFreal,ARreal,COreal],[Tgspan,Tgspan,Tgspan,Tgreal,Tgreal,Tgreal],[WLF,AR,CO,WLFreal,ARreal,COreal],Origin=self.Film.Origin,xla="T-Tg/K",yla="eta/Pas",filename="ViscoPlot")
        Plotlist([Tgspan,Tgspan,Tgspan,Tgreal,Tgreal,Tgreal],[np.log10(WLF),np.log10(AR),np.log10(CO),np.log10(WLFreal),np.log10(ARreal),np.log10(COreal)],[Tgspan,Tgspan,Tgspan,Tgreal,Tgreal,Tgreal],[np.log10(WLF),np.log10(AR),np.log10(CO),np.log10(WLFreal),np.log10(ARreal),np.log10(COreal)],Origin=False,xla="T-Tg/K",yla="eta/Pas",filename="ViscoPlot")
    def enterKris(self):
        self.Film.D12=float(self.Box0.get())

        self.Film.DAPI=np.asarray(ast.literal_eval(self.Box1.get()))
        self.Film.sigma=np.asarray(ast.literal_eval(self.Box2.get()))
        self.Film.kt=np.asarray(ast.literal_eval(self.Box3.get()))
        self.Film.g=float(self.Box5.get())
        Film1.Diffusion1D()


    def DisplayMEOS(self):
        #self.Btnb.destroy()
        self.Film.Mode="MEOS"
        #self.Btnc.destroy()
        self.BtnC1.destroy()
        self.Box0.destroy()
        self.labelDir0.destroy()
        self.BtnEOS =tk.Button(self.frameRout, text="Diffusionsmodell mit MEOS", padx=0, pady=2, font='Helvetica 7 bold', command=self.enterMEOS)
        self.BtnEOS.pack(side="top", fill="x")
        self.BtnMEOSFIT=tk.Button(self.frameRout, text="LeatFitONE", padx=0, pady=2, font='Helvetica 7 bold', command=self.LeastMEOS)
        self.BtnMEOSFIT.pack(side="top", fill="x")
        #varC3 = tk.StringVar(self.frameParam, value='2156.18581')
        varC3 = tk.StringVar(self.frameParam, value='41290')
        varFitindJ = tk.StringVar(self.frameParam, value='[3,4]')
        varTswitch = tk.StringVar(self.frameParam, value='310.464113')
        #self.Btn4 =tk.Button(self.frameRout, text="Least MEOS", padx=0, pady=2, font='Helvetica 7 bold', command=self.LeastMEOS)
        #self.Btn4.pack(side="top", fill="x")
        #self.Btn5 =tk.Button(self.frameRout, text="Least MEOSGES", padx=0, pady=2, font='Helvetica 7 bold', command=self.LeastMEOSGES)
        #self.Btn5.pack(side="top", fill="x")
        #var0 = tk.StringVar(self.frameParam, value='1.5E-13')
        #var0 = tk.StringVar(self.frameParam, value='0.8E-15')
        #var2 = tk.StringVar(self.frameParam, value='0.039')
        var0 = tk.StringVar(self.frameParam, value='3.9382E-13')
        #var0 = tk.StringVar(self.frameParam, value='6E-16')
        
        var2 = tk.StringVar(self.frameParam, value='14.034')
        #var2 = tk.StringVar(self.frameParam, value='0.12')
        #var3 = tk.StringVar(self.frameParam, value='1')
        #var5 = StringVar(self.root3, value='4')
        var5 = tk.StringVar(self.frameParam, value='1')
        #var1 = StringVar(self.root3, value="[1444940000,611018000,361484000,72075800]")
        #var1 = tk.StringVar(self.frameParam, value="5E9")
        #var1 = tk.StringVar(self.frameParam, value="5E9")
        var1 = tk.StringVar(self.frameParam, value="4E9")
        #var1 = tk.StringVar(self.frameParam, value="1E13")
        #"[456704000,321794000,122470000,26876500]" old
        #"[0.0774479630649733,0.36881875713636,1.37361808470257,10.5053372740917]"
        #"[0.00152538685660807,0.0769520156410277,0.50008174346004,4.26619792363875]" tau
        #"[1444940000,611018000,361484000,72075800]"E slightly differen
        #var3=StringVar(self.root3, value="[0.00152538685660807,0.0769520156410277,0.50008174346004,4.26619792363875]")
        #var3 = tk.StringVar(self.frameParam, value="0.361")
        var3 = tk.StringVar(self.frameParam, value="0.451")
        #var3 = tk.StringVar(self.frameParam, value="1")
        var6 = tk.StringVar(self.frameParam, value='0')
        #var6 = tk.StringVar(self.frameParam, value='0.06')
        self.Box0=tk.Entry(self.frameParam, textvariable=var0)

        Labeltext0 = tk.StringVar()
        Labeltext0.set("D12")
        labelDir0=tk.Label(self.frameParam, textvariable=Labeltext0, height=1)
        
        labelDir0.grid(row=0,column=0)
        self.Box0.grid(row=0,column=1)
        self.Box1=tk.Entry(self.frameParam, textvariable=var1)

        Labeltext1 = tk.StringVar()
        Labeltext1.set("E0")
        labelDir1=tk.Label(self.frameParam, textvariable=Labeltext1, height=1)
        labelDir1.grid(row=0,column=2)
        self.Box1.grid(row=0,column=3)
        self.Box2=tk.Entry(self.frameParam, textvariable=var2)

        Labeltext2 = tk.StringVar()
        Labeltext2.set("exponent")
        labelDir2=tk.Label(self.frameParam, textvariable=Labeltext2, height=1)
        labelDir2.grid(row=0,column=4)
        self.Box2.grid(row=0,column=5)
        self.Box3=tk.Entry(self.frameParam, textvariable=var3)

        Labeltext3 = tk.StringVar()
        Labeltext3.set("tau0")
        labelDir3=tk.Label(self.frameParam, textvariable=Labeltext3, height=1)
        labelDir3.grid(row=1,column=0)
        self.Box5=tk.Entry(self.frameParam, textvariable=var5)
        self.Box3.grid(row=1,column=1)


        Labeltext5 = tk.StringVar()
        Labeltext5.set("nJ")
        labelDir5=tk.Label(self.frameParam, textvariable=Labeltext5, height=1)
        labelDir5.grid(row=1,column=2)
        self.Box5.grid(row=1,column=3)


        Labeltext6 = tk.StringVar()
        Labeltext6.set("Xhi")
        labelDir6=tk.Label(self.frameParam, textvariable=Labeltext6, height=1)
        labelDir6.grid(row=1,column=4)
        self.Box6=tk.Entry(self.frameParam, textvariable=var6)
        self.Box6.grid(row=1,column=5)

        LabeltextFit = tk.StringVar()
        LabeltextFit.set("ListOfJumpsIndecestoFit")
        labelDirFit=tk.Label(self.frameParam, textvariable=LabeltextFit, height=1)
        labelDirFit.grid(row=0,column=8)
        self.BoxFit=tk.Entry(self.frameParam, textvariable=varFitindJ)
        self.BoxFit.grid(row=1,column=8)

        LabeltextTswitch = tk.StringVar()
        LabeltextTswitch.set("Tswitch")
        labelDirTswitch=tk.Label(self.frameParam, textvariable=LabeltextTswitch, height=1)
        labelDirTswitch.grid(row=0,column=6)
        self.BoxTswitch=tk.Entry(self.frameParam, textvariable=varTswitch)
        self.BoxTswitch.grid(row=0,column=7)

        LabeltextC3 = tk.StringVar()
        LabeltextC3.set("C3")
        labelDirC3=tk.Label(self.frameParam, textvariable=LabeltextC3, height=1)
        labelDirC3.grid(row=1,column=6)
        self.BoxC3=tk.Entry(self.frameParam, textvariable=varC3)
        self.BoxC3.grid(row=1,column=7)

    def DisplayKris(self):
        self.Film.endpoint=False
        self.Film.Mode="Kris"
        self.Film.ideal=0
        self.Film.SAFT=True
        self.Film.NET=False
        self.Film.exp_boundary()
        cwd=os.getcwd()
        molecule1=self.Film.molecule
        molecule2=self.Film.molecule2
        RHstr="90%"
        Tstr="25°C"
        DL=1-self.Film.wPolyASD
        simpath=os.path.join(cwd,"LeastFit","FeedWasser",RHstr+"_"+Tstr,"_".join(["EXP",RHstr,molecule1,molecule2,Tstr]))
        import pandas as pd
        SAFT=pd.read_excel(simpath+".xlsx")
            
        wAPISIMraw=SAFT["wAPI[-]"].values[::-1]
        wwSIMraw=SAFT["ww[-]"].values[::-1]
        alpharaw=(DL-wAPISIMraw)/(1-wAPISIMraw)
        rhoASD=self.Film.rho0ASD
        
        self.Film.alphavec=alpharaw
        self.Film.rho2GGWvec=wwSIMraw/(1-wwSIMraw)*rhoASD
        self.Film.RHraw=0.9

        self.Film.Mode="Kris"
        self.Btnc.destroy()
        self.BtnCh.destroy()
        self.Box0.destroy()
        self.labelDir0.destroy()
        self.BtnEOS =tk.Button(self.frameRout, text="Diffusionsmodell mit Kris", padx=0, pady=2, font='Helvetica 7 bold', command=self.enterKris)
        self.BtnEOS.pack(side="top", fill="x")
        self.BtnMEOSFIT=tk.Button(self.frameRout, text="LeatFitONE", padx=0, pady=2, font='Helvetica 7 bold', command=lambda:print("Not Implemented, YET!"))
        self.BtnMEOSFIT.pack(side="top", fill="x")

        var0 = tk.StringVar(self.frameParam, value='1E-13')
        var1 = tk.StringVar(self.frameParam, value="6.5E-19")
        var2 = tk.StringVar(self.frameParam, value='2.97E-2')
        var3 = tk.StringVar(self.frameParam, value="2E-12")
        var5 = tk.StringVar(self.frameParam, value='3.92')

        self.Box0=tk.Entry(self.frameParam, textvariable=var0)

        Labeltext0 = tk.StringVar()
        Labeltext0.set("D12")
        labelDir0=tk.Label(self.frameParam, textvariable=Labeltext0, height=1)
        labelDir0.grid(row=0,column=0)
        self.Box0.grid(row=0,column=1)
        self.Box1=tk.Entry(self.frameParam, textvariable=var1)

        Labeltext1 = tk.StringVar()
        Labeltext1.set("DAPI")
        labelDir1=tk.Label(self.frameParam, textvariable=Labeltext1, height=1)
        labelDir1.grid(row=0,column=2)
        self.Box1.grid(row=0,column=3)
        self.Box2=tk.Entry(self.frameParam, textvariable=var2)

        Labeltext2 = tk.StringVar()
        Labeltext2.set("sigma")
        labelDir2=tk.Label(self.frameParam, textvariable=Labeltext2, height=1)
        labelDir2.grid(row=0,column=4)
        self.Box2.grid(row=0,column=5)
        self.Box3=tk.Entry(self.frameParam, textvariable=var3)

        Labeltext3 = tk.StringVar()
        Labeltext3.set("kt")
        labelDir3=tk.Label(self.frameParam, textvariable=Labeltext3, height=1)
        labelDir3.grid(row=1,column=0)
        self.Box5=tk.Entry(self.frameParam, textvariable=var5)
        self.Box3.grid(row=1,column=1)


        Labeltext5 = tk.StringVar()
        Labeltext5.set("g")
        labelDir5=tk.Label(self.frameParam, textvariable=Labeltext5, height=1)
        labelDir5.grid(row=1,column=2)
        self.Box5.grid(row=1,column=3)

        


    def main(self):
     """The main function initializes the GUI-window."""
     self.root3 = tk.Tk()


     self.root3.attributes("-topmost", True)
     #root.geometry("350x150+450+280")  # "width x height + position right + position down"
     self.root3.title("Diffusion Model")
     #root.iconbitmap("TU.ico")#here the TU Dortmund Icon has be in the same folder as the python script and
                              #the image as to be an icon type.
     self.root3.grid_columnconfigure((0,1), weight=1)
     Filmismadeof1=self.Film.molecule
     Filmismadeof2=" and "+self.Film.molecule2+" with an API Feed of "+str(round((1-self.Film.wPolyASD)*100))+"%"  if self.Film.wPolyASD!=1 else ""
     Filmismadeof=Filmismadeof1+Filmismadeof2
     self.label1=tk.Label(self.root3,text='Your Film consist of :'+Filmismadeof,bg='gold',fg='black',anchor="w", justify="left")
     self.label3=tk.Label(self.root3,text='Default Jump:  0',bg='gold',fg='black',anchor="w", justify="left")

     #self.label2=tk.Label(self.root3,text='Current Path is '+self.Film.DATA_FILE,bg='gold',fg='black')
     self.label4=tk.Label(self.root3,text="Default Mode: Absorption",bg='red',fg='black',anchor="w", justify="left")
     self.label5=tk.Label(self.root3,text="Default Mode: Ideal",bg='green',fg='black',anchor="w", justify="left")
     self.label6=tk.Label(self.root3,text="Default Mode: Activity from Isotherm",bg='green',fg='black',anchor="w", justify="left")
     self.label7=tk.Label(self.root3,text="Default Mode: NET-GP Off",bg='green',fg='black',anchor="w", justify="left")
     self.label8=tk.Label(self.root3,text="Default Mode: Origin Off",bg='green',fg='black',anchor="w", justify="left")
     self.label9=tk.Label(self.root3,text="Default Mode: Actual RH Off",bg='green',fg='black',anchor="w", justify="left")

     self.label1.pack(side="top", fill="x")
     self.label3.pack(side="top", fill="x")
     #self.label2.pack(side="top", fill="x")
     self.label4.pack(side="top", fill="x")
     self.label5.pack(side="top", fill="x")
     self.label6.pack(side="top", fill="x")
     self.label7.pack(side="top", fill="x")
     self.label8.pack(side="top", fill="x")
     self.label9.pack(side="top", fill="x")

     #Btn0 = tk.Button(self.root3, text="enter", padx=0, pady=2, font='Helvetica 7 bold', command=self.enter)




     #Btn3 = tk.Button(self.root3, text="Diffusionsmodell mit mechanischer Zustandsgleichung", padx=0, pady=2, font='Helvetica 7 bold', command=self.Diffusion1DMEOS)
     self.frame=tk.LabelFrame(self.root3,text="Sorption Steps",font='Helvetica 7 bold')
     self.framefit=tk.LabelFrame(self.root3,text="Fits & Predictions",font='Helvetica 7 bold')
     self.frameModes=tk.LabelFrame(self.root3,text="Mode Changer",font='Helvetica 7 bold')
     self.frameRout=tk.LabelFrame(self.root3,text="Routines",font='Helvetica 7 bold')
     self.frameParam=tk.LabelFrame(self.root3,text="Parameter",font='Helvetica 7 bold')
     self.frameFilm=tk.LabelFrame(self.root3,text="Custom film",font='Helvetica 7 bold')
     self.frame.pack(side="top", fill="x")
     self.framefit.pack(side="top", fill="x")
     self.frameModes.pack(side="top", fill="x")
     self.frameRout.pack(side="top", fill="x")
     self.frameParam.pack(side="top", fill="x")
     self.frameFilm.pack(side="top", fill="x")

     Btna = tk.Button(self.frameRout, text="Model current step", padx=0, pady=2, font='Helvetica 7 bold', command=self.Diffusion1D)
     self.Btnb = tk.Button(self.frameRout, text="Least-squares", padx=0, pady=2, font='Helvetica 7 bold', command=self.Least)
     self.Btnc = tk.Button(self.frameRout, text="Prediction", padx=0, pady=2, font='Helvetica 7 bold', command=self.Prediction)
     #Btnd = tk.Button(self.frameModes, text="Desorption/Absorption", padx=0, pady=2, font='Helvetica 7 bold', command=self.Desorption)
     #Btnd = tk.Checkbutton(self.frameModes, text="Desorption on /Absorption off", padx=0, pady=2, font='Helvetica 7 bold', command=self.Desorption)
     #Btne = tk.Button(self.frameModes, text="Ideal/Non-Ideal", padx=0, pady=2, font='Helvetica 7 bold', command=self.Ideal)
     #Btne = tk.Checkbutton(self.frameModes, text="Ideal on / Non-Ideal on", padx=0, pady=2, font='Helvetica 7 bold', command=self.Ideal)
     #Btnf = tk.Button(self.frameModes, text="Activity from PCSAFT/Isotherm", padx=0, pady=2, font='Helvetica 7 bold', command=self.ActivityFromIsotherm)
     #Btnf = tk.Checkbutton(self.frameModes, text="Activity from PC-SAFT on / Isotherm Off", padx=0, pady=2, font='Helvetica 7 bold', command=self.ActivityFromIsotherm)
     #Btng = tk.Button(self.frameModes, text="NET-GP", padx=0, pady=2, font='Helvetica 7 bold', command=self.NET)
     #Btng = tk.Checkbutton(self.frameModes, text="NET-GP on / NET-GP off", padx=0, pady=2, font='Helvetica 7 bold', command=self.NET)
     Btnh = tk.Button(self.frameRout, text="Recaluclate the PC-SAFT modeling", padx=0, pady=2, font='Helvetica 7 bold', command=self.RecalcIsotherm)
     #Btni = tk.Button(self.frameModes, text="Origin On/Off", padx=0, pady=2, font='Helvetica 7 bold', command=self.Origin)
     #Btni = tk.Checkbutton(self.frameModes, text="Origin files on / Origin files off", padx=0, pady=2, font='Helvetica 7 bold', command=self.Origin)
     #Btnj = tk.Button(self.frameModes, text="Actual On/Off", padx=0, pady=2, font='Helvetica 7 bold', command=self.Actual)
     #Btnj = tk.Checkbutton(self.frameModes, text="Actual RH on / Target RH off", padx=0, pady=2, font='Helvetica 7 bold', command=self.Actual)
     BtnSetC = tk.Button(self.frameFilm, text="Set conditions", padx=0, pady=2, font='Helvetica 7 bold', command=self.SetConditions)
     picks=['Desorption on /Absorption off', 'Ideal on / Non-Ideal on', "Activity from PC-SAFT on / Isotherm Off", 'NET-GP on / NET-GP off',"Origin files on / Origin files off","Actual RH on / Target RH off"]
     commands=[self.Desorption,self.Ideal,self.ActivityFromIsotherm,self.NET,self.Origin,self.Actual]

     #self.BtnX =[tk.Button(self.frame, text="Jump "+str(i), padx=0, pady=2, font='Helvetica 7 bold', command=partial(self.indexjump, i)) for i in range(self.Film.nHum)]
     #[val.pack(side="top", fill="x") for i,val in enumerate(self.BtnX)]
     #current_value = tk.StringVar(value=0)
     self.spinbox = tk.Spinbox(self.frame, from_ = 1, to = self.Film.nHum,increment = 1,command=partial(self.indexjump, 1))
     self.spinbox.pack(side="top", fill="x")
    #  listNodes = tk.Listbox(self.frame, font="Helvetica 7 bold")
    #  listNodes.pack(side="top", fill="x")

    #  scrollbar = tk.Scrollbar(self.frame, orient="vertical")
    #  scrollbar.config(command=listNodes.yview)
    #  scrollbar.pack(side="top", fill="x")
    #  listNodes.config(yscrollcommand=scrollbar.set)

    #  for i in range(self.Film.nHum):
    #     listNodes.insert(tk.END, "Jump "+str(i))

    
     self.BtnCh=tk.Button(self.frameModes, text="Display MEOS", padx=0, pady=2, font='Helvetica 7 bold', command=self.DisplayMEOS)
     self.BtnC1=tk.Button(self.frameModes, text="Display Kris", padx=0, pady=2, font='Helvetica 7 bold', command=self.DisplayKris)
     self.BtnAuf=tk.Button(self.frameModes, text="Auflösungsmodell", padx=0, pady=2, font='Helvetica 7 bold', command=self.Dissolution1D)
     self.BtnCh.pack(side="top", fill="x")
     self.BtnC1.pack(side="top", fill="x")
     self.BtnAuf.pack(side="top", fill="x")
     #Btn0.pack(side="top", fill="x")
     Btna.pack(side="top", fill="x")
     self.Btnb.pack(side="top", fill="x")
     self.Btnc.pack(side="top", fill="x")
    #  Btnd.pack(side="top", fill="x")
    #  Btnd["anchor"]="w"
    #  Btne.pack(side="top", fill="x")
    #  Btne["anchor"]="w"
    #  Btnf.pack(side="top", fill="x")
    #  Btnf["anchor"]="w"
    #  Btng.pack(side="top", fill="x")
    #  Btng["anchor"]="w"
    #  Btni.pack(side="top", fill="x")
    #  Btni["anchor"]="w"
    #  Btnj.pack(side="top", fill="x")
    #  Btnj["anchor"]="w"
     lng = Checkbar(self.frameModes,picks=picks,commands=commands)
     lng.pack(side="top", fill="x")


     var0 = tk.StringVar(self.frameParam, value='1E-14')
     self.Box0=tk.Entry(self.frameParam, textvariable=var0)
     self.Box0.grid(row=0,column=1)
     Labeltext0 = tk.StringVar()
     Labeltext0.set("D12")
     self.labelDir0=tk.Label(self.frameParam, textvariable=Labeltext0, height=1)
     self.labelDir0.grid(row=0,column=0)
     varA = tk.StringVar(self.frameParam, value=self.Film.A0)
     varB = tk.StringVar(self.frameParam, value=self.Film.solvent)
     self.BoxA=tk.Entry(self.frameFilm, textvariable=varA)
     self.BoxA.grid(row=0,column=1)

     LabeltextA = tk.StringVar()
     LabeltextA.set("Area")

     self.labelDirA=tk.Label(self.frameFilm, textvariable=LabeltextA, height=1)
     self.labelDirA.grid(row=0,column=0)

     LabeltextB = tk.StringVar()
     LabeltextB.set("Solvent")
     self.BoxB=tk.Entry(self.frameFilm, textvariable=varB)
     self.BoxB.grid(row=0,column=1)
     self.labelDirB=tk.Label(self.frameFilm, textvariable=LabeltextB, height=1)
     self.labelDirB.grid(row=0,column=0)

     varRHd = tk.StringVar(self.frameParam, value="[25,50,75,50,25,0]")
     varwd = tk.StringVar(self.frameParam, value="[0.1,0.2,0.3,0.2,0.1,0]")
     varL0d = tk.StringVar(self.frameParam, value="1E-5")
     vartd = tk.StringVar(self.frameParam, value="[100,100,100,100,100,100]")

     #varRHd = tk.StringVar(self.frameParam, value="[100]")
     #varwd = tk.StringVar(self.frameParam, value="[0.122930896]")
     #varL0d = tk.StringVar(self.frameParam, value="1E-3")
     #vartd = tk.StringVar(self.frameParam, value="[68400]")

     self.BoxRHd=tk.Entry(self.frameFilm,textvariable=varRHd)
     self.BoxRHd.grid(row=0,column=0)
     self.Boxwd=tk.Entry(self.frameFilm,textvariable=varwd)
     self.Boxwd.grid(row=1,column=0)
     self.Boxtd=tk.Entry(self.frameFilm,textvariable=vartd)
     self.Boxtd.grid(row=2,column=0)
     self.BoxL0d=tk.Entry(self.frameFilm,textvariable=varL0d)
     self.BoxL0d.grid(row=3,column=0)






     BtnA = tk.Button(self.frameFilm, text="New Film", padx=0, pady=2, font='Helvetica 7 bold', command=self.NewFilm)
     BtnSetC.grid(row=2,column=1)


     BtnA.grid(row=3,column=1)
     self.fig,self.ax=plt.subplots(figsize=(3, 1))
     self.ax.set_ylabel("RH")
     self.ax.set_ylim([0,1])
     self.ax.set_yticks([0,0.2,0.4,0.6,0.8,1])
     self.ax.grid()
     self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame)
     
     matplotlib.use("TkAgg")
     x=np.arange(0,self.Film.nHum+1,1)
     y=np.hstack((0,self.Film.Feuchtereal))   
     self.line1, =self.ax.step(x,y,'b',where='pre')
     self.line2, =self.ax.step(x[self.Film.indJ:self.Film.indJ+2],y[self.Film.indJ:self.Film.indJ+2],'r',where='pre')
     self.line3=[]
     for index in range(len(x)-1):
        self.line3.append(self.ax.text(x[index]+0.1, y[index]-0.02, x[index+1], size=7))
     self.canvas.get_tk_widget().pack(side="bottom")
     self.root3.mainloop()

window33=window3(Film1,Data)
