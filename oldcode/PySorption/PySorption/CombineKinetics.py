import pandas as pd
from os.path import join
from os import getcwd
import tkinter
from tkinter import *
from tkinter.filedialog import askopenfilenames
import matplotlib.pyplot as plt
from ReadXandYdatafromExcel import CombineData,CombineKinetikData
import itertools
import numpy as np
import os
class window():
    def __init__(self):
        self.pathlist=[]
        self.strt="tges[min]"
        self.strx="wges[-]"
        self.strtsim="tgessim[min]"
        self.strxsim="wgessim[-]"
        self.BoxData=[]
        self.kinetics=1
        self.main()
        
    
    def EnterBtn(self):
        self.strt=self.Box0.get()
        self.strx=self.Box1.get()
        self.strtsim=self.Box2.get()
        self.strxsim=self.Box3.get()
        self.kinetics=int(self.Box4.get())

        dirname=os.path.dirname(self.pathlist[0])
        if self.kinetics:
            t,m,mave,mstd=CombineKinetikData(self.strt,self.strx,self.pathlist)
            tsim,msim,msimave,msimstd=CombineKinetikData(self.strtsim,self.strxsim,self.pathlist)
        else:
            t,m,mave,mstd=CombineData(self.strt,self.strx,self.pathlist)
            tsim,msim,msimave,msimstd=CombineData(self.strtsim,self.strxsim,self.pathlist)


        filenames=[val.split("/")[-1] for i,val in enumerate(self.pathlist)]
        try:
            Filmnumbers=[val[val.find("Film"):val.find("Film")+7].replace("_","") for i,val in enumerate(filenames)] if filenames[0].find("Film")!=-1 else "unknown" 
        except:
            Filmnumbers=[ val.get() for i,val in enumerate(self.BoxData)]
        
        head=list(np.asarray([[self.strt+str(i),self.strx+str(i)] for i in range(len(m))]).flatten())
        headsim=list(np.asarray([[self.strtsim+str(i),self.strxsim+str(i)] for i in range(len(msim))]).flatten())
        tm=tuple(zip(t,m))
        tmsim=tuple(zip(tsim,msim))
        arr=[]
        for val in tm:
            arr+=val
        arr+=[mave]+[mstd]
        for val in tmsim:
            arr+=val
        arr+=[msimave]+[msimstd]
        datasave=pd.DataFrame(tuple(arr))
        header=head+[self.strx+"ave"]+[self.strx+"std"]+headsim+[self.strxsim+"ave"]+[self.strxsim+"std"]
    
        datawrite=pd.concat(tuple([pd.DataFrame(val) for i,val in enumerate(arr)]),ignore_index=True, axis=1)
        datawrite.columns=header
        name="Combined_"+self.BoxData[0].get()+"_"+"_".join(Filmnumbers)
        parentname=os.path.dirname(dirname)
        datawrite.to_excel(os.path.join(parentname,name+".xlsx"))
        self.t,self.m=t+[t[0]],m+[mave]
        self.tsim,self.msim=tsim+[tsim[0]],msim+[msimave]
        self.terr,self.merr=[np.zeros_like(mstd)]*len(m+[mave]),[mstd]*len(m+[mave])
        self.Btn3.config(state="normal")
        self.name=name
    def OriginBtn(self):
        
        from PyPlorigin import PyOrigin
        x,y,xs,ys,xe,ye=PyOrigin.MusterDictOrigin()
        x["values"],y["values"]=self.t,self.m
        xs["values"],ys["values"]=self.tsim,self.msim
        xe["values"],ye["values"]=self.terr,self.merr
        xlabel="t/min"
        ylabel="w/-"
        filename=self.name
        PyOrigin.Plot(x,y,xs,ys,xe,ye,xlabel,ylabel,filename,Vorlage=self.BoxVorlage.get()+".opj")
    def searchdirectory(self):
        self.root.withdraw()
        self.pathlist= list(askopenfilenames())
        self.root.deiconify()
        #self.pathlist.append(path)
        self.label2.config(text='Current Paths are \n'+"\n".join([str(i+1)+". "+ val for i,val in enumerate(self.pathlist)]))
        #temp=path.split("/")[-1]
        #filmname=temp[temp.find("Film"):temp.find("Film")+7].replace("_","") if temp.find("Film")!=-1 else "unknown"
        filmname=os.path.basename(os.path.normpath(os.path.dirname(self.pathlist[0])))
        varData = StringVar(self.root, value=filmname)
        BoxData=Entry(self.root, textvariable=varData)
        self.BoxData.append(BoxData)
        #self.BoxData[len(self.pathlist)-1].pack(side="top", fill="x")
        BoxData.pack(side="top", fill="x")

    
    def reset(self):
        self.pathlist=[]
        [val.destroy() for i,val in enumerate(self.BoxData)]
        self.BoxData=[]
        self.label2.config(text='Current Path is not chosen click browse')
    
    def main(self):
     """ Main Window for GUI"""
     self.root = Tk()
     self.root.attributes("-topmost", True)

     #root.geometry("350x150+450+280")  # "width x height + position right + position down"
     self.root.title("Combine Files")
     #root.iconbitmap("TU.ico")#here the TU Dortmund Icon has be in the same folder as the python script and 

     var0 = StringVar(self.root, value=self.strt)
     self.Box0=Entry(self.root, textvariable=var0)
     self.Box0.pack(side="top", fill="x")
     Labeltext0 = StringVar()
     Labeltext0.set("xHeader")
     labelDir0=Label(self.root, textvariable=Labeltext0, height=4)
     labelDir0.pack(side="left")
     
     var1 = StringVar(self.root, value=self.strx)
     self.Box1=Entry(self.root, textvariable=var1)
     self.Box1.pack(side="top", fill="x")
     Labeltext1 = StringVar()
     Labeltext1.set("yHeader")
     labelDir1=Label(self.root, textvariable=Labeltext1, height=4)
     labelDir1.pack(side="left")
     
     var2 = StringVar(self.root, value=self.strtsim)
     self.Box2=Entry(self.root, textvariable=var2)
     self.Box2.pack(side="top", fill="x")
     Labeltext2 = StringVar()
     Labeltext2.set("xsimHeader")
     labelDir2=Label(self.root, textvariable=Labeltext2, height=4)
     labelDir2.pack(side="left")
     
     var3 = StringVar(self.root, value=self.strxsim)
     self.Box3=Entry(self.root, textvariable=var3)
     self.Box3.pack(side="top", fill="x")
     Labeltext3 = StringVar()
     Labeltext3.set("ysimHeader")
     labelDir3=Label(self.root, textvariable=Labeltext3, height=4)
     labelDir3.pack(side="left")

     var4 = StringVar(self.root, value=self.kinetics)
     self.Box4=Entry(self.root, textvariable=var4)
     self.Box4.pack(side="top", fill="x")
     Labeltext4 = StringVar()
     Labeltext4.set("Kinetics?")
     labelDir4=Label(self.root, textvariable=Labeltext4, height=4)
     labelDir4.pack(side="left")
     
     self.Btn0 = Button(self.root, text="enter", padx=50, pady=2, font='Helvetica 8 italic', command=self.EnterBtn)
     self.Btn0.pack(side="left")
     self.Btn1 = Button(self.root, text="browse", padx=50, pady=2, font='Helvetica 8 italic', command=self.searchdirectory)
     self.Btn1.pack(side="left")
     self.Btn2 = Button(self.root, text="reset", padx=50, pady=2, font='Helvetica 8 italic', command=self.reset)
     self.Btn2.pack(side="left")
     self.Btn3 = Button(self.root, text="origin", padx=50, pady=2, font='Helvetica 8 italic', command=self.OriginBtn,state="disabled")
     self.Btn3.pack(side="left")

     varVorlage = StringVar(self.root, value="Origin")
     self.BoxVorlage=Entry(self.root, textvariable=varVorlage)
     self.BoxVorlage.pack(side="bottom", fill="x")
     LabeltextVorlage = StringVar()
     LabeltextVorlage.set("Name der Vorlage")
     labelDirVorlage=Label(self.root, textvariable=LabeltextVorlage, height=4)
     labelDirVorlage.pack(side="bottom")


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


