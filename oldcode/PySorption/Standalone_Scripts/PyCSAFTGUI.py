from PyCSAFTFinal2 import Mixture
import PureDataBase
import os
from tkinter import *
import tkinter.font as font 
from functools import partial
import matplotlib.pyplot as plt
originalpath=os.getcwd()
from os.path import join
#originalpath='C:/Users/Moritz/Documents/Test'

def DataBaseImport(*args):
    """Enter components as tuple, get tuple of dicts"""
    return tuple(getattr(PureDataBase,val) for i,val in enumerate(args))


class PCSAFTGUI:
    def __init__(self):
        self.complist=[]
        self.main()

    def addcomponents(self,x):
        self.complist.append(x) if x not in self.complist and len(self.complist)<3 else None
        self.Btn2["state"]="disabled" if  len(self.complist)<3 else "active"
        self.Btn4["state"]="disabled" if  len(self.complist)<3 else "active"
        self.Btn5["state"]="disabled" if  len(self.complist)<3 else "active"
        self.label2.config(text="_".join(self.complist))
    def VLE(self):
        self.Tuples()
        mix=Mixture(*self.comp)
        mix.Isotherm(float(self.Box0.get()))
    def MultIsotherm(self):
        self.Tuples()
        mix=Mixture(*self.comp)
        mix.MultIsothermOrigin(T=298.15)
        
    def Gibbs(self):
        self.Tuples()
        mix=Mixture(*self.comp)
        mix.nASD=int(self.Box1.get())
        mix.nVLE=int(self.Box2.get())
        mix.nLLE=int(self.Box3.get())
        mix.nGplot=int(self.Box4.get())
        mix.GibbsDiagramm(1)

    def Tuples(self):
        self.comp=tuple([getattr(PureDataBase,val) for i,val in enumerate(self.complist)])#
    def Reset(self):
        self.complist=[]
        self.label2.config(text="_".join(self.complist))
    def main(self):
     """The main function initializes the GUI-window."""
     self.root2 = Tk()
     self.root2.attributes("-topmost", True)
     #root.geometry("350x150+450+280")  # "width x height + position right + position down"
     self.root2.title("Diffusion Model")
     #root.iconbitmap("TU.ico")#here the TU Dortmund Icon has be in the same folder as the python script and 
                              #the image as to be an icon type.
     self.label2=Label(self.root2,text=' ',bg='gold',fg='blue')
     
     Database=[val for i,val in enumerate(dir(PureDataBase)) if isinstance(getattr(PureDataBase,val),dict)]
     Database.remove("__builtins__")
     Btn1 = [Button(self.root2, text=val,font='Helvetica 8 italic', command=partial(self.addcomponents,val)) for i,val in enumerate(Database)]
     self.Btn2 = Button(self.root2,state="disabled", text="VLE", font='Helvetica 8 italic', command=self.VLE)   
     self.Btn4 = Button(self.root2,state="disabled", text="Gibbs", font='Helvetica 8 italic', command=self.Gibbs)   
     self.Btn5 = Button(self.root2,state="disabled", text="MultIso", font='Helvetica 8 italic', command=self.MultIsotherm)   

     [val.grid(column=i%8,row=i//8) for i ,val in enumerate(Btn1)]
     self.Btn2.grid()
     self.Btn4.grid()
     self.Btn5.grid()
     self.Btn3=Button(self.root2, text="Reset", font='Helvetica 8 italic', command=self.Reset)
     self.Btn3.grid()
     self.label2.grid()
     
     var0 = StringVar(self.root2, value="1")
     self.Box0=Entry(self.root2, textvariable=var0)
     self.Box0.grid()
     var1 = StringVar(self.root2, value="21")
     self.Box1=Entry(self.root2, textvariable=var1)
     self.Box1.grid()
     var2 = StringVar(self.root2, value="20")
     self.Box2=Entry(self.root2, textvariable=var2)
     self.Box2.grid()
     var3 = StringVar(self.root2, value="20")
     self.Box3=Entry(self.root2, textvariable=var3)
     self.Box3.grid()
     var4 = StringVar(self.root2, value="10")
     self.Box4=Entry(self.root2, textvariable=var4)
     self.Box4.grid()
     Labeltext0 = StringVar()
     Labeltext0.set("wFAPI")
     labelDir0=Label(self.root2, textvariable=Labeltext0, height=4)
     labelDir0.grid()
     self.root2.mainloop()

window1=PCSAFTGUI()

