import numpy as np
import matplotlib.pyplot as plt
import tkinter
from tkinter import *
from tkinter import filedialog

def EndpointOptimization(mend=22.5,L=300,path="InQuestion.xlsx"):
    """Experimental Data must be supplied in the form of an xls or xlsx file with columns:
        t[min]   m[mg]
        300      1
        ...      ...
        600      5
        The Mass mGGW that sould be known and is supplied via the GUI.
        The Thickness is only important for the Diffusion coefficient"""
    from lmfit import Parameters, minimize
    from lmfit.printfuncs import report_fit

    fit_params = Parameters() 
    fit_params.add('D', value=1, min=.001, max=10000,vary=1)
    fit_params.add('tend', value=8, min=0, max=100,vary=1)
    #fit_params.add('mend', value=22.5, min=0, max=30,vary=0)
    #mend=fit_params['mend'].value
    import pandas as pd
    INQ=pd.read_excel(path)
    t=INQ["t[min]"].values
    deltat=np.round(np.average(t[1:]-t[:-1]),2)
    print(deltat)

    t=t-t[0]
    m=INQ["m[mg]"].values
    mp=(m-m[0])/(mend-m[0])
    mpend=1
    def Cranc1(params,t,mp):
        tend=np.asarray(params["tend"].value)*100
        t2=np.hstack((t,tend))*60
        mp2=np.hstack((mp,1))

        n=12
        Forier=0
        LCranc=2*L*1E-6
        D=np.asarray(params["D"].value)*1E-14
        for cc in range(n):
            Forier=Forier+1/(2*cc+1)**2*np.exp(-((2*cc+1)**2*np.pi**2*D*t2)/(LCranc**2))
            mperc=1-8/np.pi**2*Forier
        weight=np.ones_like(mp2)
        weight[-1]=weight[-1]
        return (mperc-mp2)*weight

    #Cranc2=lambda params,t,mp:np.sum(Cranc1(params,t,mp))**2
    out = minimize(Cranc1,args=(t,mp), params=fit_params,method="powell")
    Dopt=out.params['D'].value*1E-14
    tend=out.params['tend'].value*100

    
    print(Dopt)
    t2=np.hstack((t,tend))
    mp2=np.hstack((mp,mpend))
    
    
    
    t3=np.linspace(t[0],tend,1000)
    weight=np.ones_like(t3)
    weight[-1]=weight[-1]
    opt=Cranc1(out.params,t3[:-1],np.zeros_like(t3[:-1]))/weight+np.hstack((np.zeros_like(t3[:-1]),1))
    
    weight2=np.ones_like(t2)
    weight2[-1]=weight2[-1]
    optsametime=Cranc1(out.params,t2[:-1],mp)/weight2+mp2
    
    mplotexp=mp2*(-m[0]+mend)+m[0]
    mplot=opt*(-m[0]+mend)+m[0]
    
    msametimesim=optsametime*(-m[0]+mend)+m[0]
    
    texpsim=np.arange(t[0],tend,deltat)
    optexpsim=(Cranc1(out.params,texpsim,np.zeros_like(texpsim))+np.hstack((np.zeros_like(texpsim),1)))*(-m[0]+mend)+m[0]
    
    t99=t3[next(i for i,v in enumerate(opt-0.99>0) if v)]
    t98=t3[next(i for i,v in enumerate(opt-0.98>0) if v)]
    t95=t3[next(i for i,v in enumerate(opt-0.95>0) if v)]
    
    fig,ax=plt.subplots()
    ax.plot(t2,mplotexp,'kx')
    ax.plot(t3,mplot,'k-')
    ax.plot(texpsim,optexpsim[:-1],'b-')
    ax.set_xlabel("t [min]")
    ax.set_ylabel("m [mg]")
    ax.plot([t99,t99],[m[0],mplot[t3==t99][0]],'b--')
    ax.plot([0,t99],[mplot[t3==t99][0],mplot[t3==t99][0]],'b--',label='99% of GGW')
    ax.plot([t98,t98],[m[0],mplot[t3==t98][0]],'g--')
    ax.plot([0,t98],[mplot[t3==t98][0],mplot[t3==t98][0]],'g--',label='98% of GGW')
    ax.plot([t95,t95],[m[0],mplot[t3==t95][0]],'r--')
    ax.plot([0,t95],[mplot[t3==t95][0],mplot[t3==t95][0]],'r--',label='95% of GGW')
    ax.legend()
    datawrite=pd.concat((pd.DataFrame(t2),pd.DataFrame(mplotexp),pd.DataFrame(msametimesim),pd.DataFrame(texpsim),pd.DataFrame(optexpsim),pd.DataFrame(np.asarray([Dopt])),pd.DataFrame(np.asarray([t99])),pd.DataFrame(np.asarray([t98])),pd.DataFrame(np.asarray([t95]))), ignore_index=True, axis=1)
    datawrite.columns=["t[min]","mExp[mg]","mSim[Exp]","tpred[min]","mpred[mg]","Dopt[m^2/s]","t99[min]","t98[min]","t95[min]"]
    datawrite.to_excel(path.replace(".xlsx","Out.xlsx"))
class window():
    def __init__(self):
        self.path="InQuestion.xlsx"
        self.main()
        
        
    
    def EnterBtn(self):
        mend=float(self.Box0.get())
        L=float(self.Box1.get())
        plt.close("all")
        EndpointOptimization(mend,L,self.path)
        
        
    def searchdirectory(self):
        self.root.withdraw()
        self.path= filedialog.askopenfilename()
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
     var0 = StringVar(self.root, value='22.5')
     self.Box0=Entry(self.root, textvariable=var0)
     self.Box0.pack(side="top", fill="x")
     Labeltext0 = StringVar()
     Labeltext0.set("mGGW[mg]")
     labelDir0=Label(self.root, textvariable=Labeltext0, height=4)
     labelDir0.pack(side="left")
     
     var1 = StringVar(self.root, value='300')
     self.Box1=Entry(self.root, textvariable=var1)
     self.Box1.pack(side="top", fill="x")
     Labeltext1 = StringVar()
     Labeltext1.set(r'Dicke[mikrometer]')
     labelDir1=Label(self.root, textvariable=Labeltext1, height=4)
     labelDir1.pack(side="left")
     
     self.Btn0 = Button(self.root, text="enter", padx=50, pady=2, font='Helvetica 8 italic', command=self.EnterBtn)
     self.Btn0.pack(side="left")
     self.Btn1 = Button(self.root, text="browse", padx=50, pady=2, font='Helvetica 8 italic', command=self.searchdirectory)
     self.Btn1.pack(side="left")
     
     self.label2=Label(self.root,text='Current Path is not chosen click browse',bg='gold',fg='blue')
     self.label2.pack(side="top")
     
     self.root.mainloop()  
if __name__=="__main__":
    #EndpointOptimization()
    window1=window()
    