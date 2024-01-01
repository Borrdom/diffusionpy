import SchemDraw
import SchemDraw.elements as elm
from SchemDraw import Segment
from matplotlib.cbook import get_sample_data
import matplotlib.pyplot as plt



import os
def AddSpring(Eeq,i,d):
    d.add(elm.Line('left'))
    d.add(elm.Line('left'))
    #d.add(elm.Resistor(d='left', label="{:.1e}".format(Eeq)+'\n[Pa]',botlabel="$E_{eq}$"))
    
def AddMaxwell(EJ,etaJ,i,d):
    d.add(elm.Line('up'))
    d.push()
    d.add(elm.Resistor(d='right', label="{:.1e}".format(float(EJ))+'\n[Pa]',botlabel="$E_{"+str(i)+"}$"))
    n=d.add(elm.RBox(d='right', label="{:.1e}".format(float(etaJ))+'\n[Pas]',botlabel="$\eta_{"+str(i)+"}$"))
    n.segments.append(Segment([[0, 0], [0.5*1.41, 0]]))
    n.segments.append(Segment([[0.5*1.41, -0.1*1.41], [0.5*1.41, 0.1*1.41]]))
    d.add(elm.Line(d='down'))
    d.pop()
def AddMaxwell2(EJ,etaJ,i,d):
    
    d.push()
    d.add(elm.Resistor(d='right', botlabel="$E_{"+str(i)+"}$"))
    n=d.add(elm.RBox(d='right', botlabel="$\eta_{"+str(i)+"}$"))
    n.segments.append(Segment([[0, 0], [0.5*1.41, 0]]))
    n.segments.append(Segment([[0.5*1.41, -0.1*1.41], [0.5*1.41, 0.1*1.41]]))
    d.add(elm.Line(d='down'))
    d.pop()
def EMPTY(EJ,etaJ,i,d):
    d.push()
    d.add(elm.Dot())
    d.add(elm.Line(d='left'))
    d.add(elm.Line(d='right'))
    d.add(elm.Line(d='right'))
    d.add(elm.Dot())
    d.add(elm.Line(d='right'))
    d.add(elm.Dot())
    d.add(elm.Line(d='right'))
    d.add(elm.Line(d='left'))
    d.add(elm.Line(d='down'))
    d.pop()
def DrawSheet(Eeq,EJvec,etaJvec):
    d = SchemDraw.Drawing()
    #nJ=len(EJvec)
    AddSpring(Eeq,0,d)
    hmm=[AddMaxwell(val,etaJvec[i],(i),d) for i,val in enumerate(EJvec)]
    schemfig = d.draw()
    #schemfig.save(os.getcwd()+"\Model.png",dpi=100)
    fig,ax=plt.subplots(figsize=(4,5))
    #im=plt.imread(os.getcwd()+"\Model.png")
    #ax.imshow(im)
    #ax.axis("off")
    #plt.show()
    #schemfig.show()
    return fig, ax
def DrawSheetG(EJvec,etaJvec):
    d = SchemDraw.Drawing()
    #nJ=len(EJvec)
    #AddSpring(Eeq,0,d)
    AddMaxwell2(1,1,0,d)
    d.add(elm.Line('up'))
    EMPTY(1,1,0,d)
    d.add(elm.Line('up'))
    AddMaxwell2(1,1,0,d)
    #schemfig = d.draw()
    #schemfig.save(os.getcwd()+"\Model.png",dpi=100)
    fig,ax=plt.subplots(figsize=(4,5))
    #im=plt.imread(os.getcwd()+"\Model.png")
    #ax.imshow(im)
    #ax.axis("off")
    #plt.show()
    #schemfig.show()
    return fig, ax

if __name__=="__main__":
    import numpy as np
    nJ=5
    Eeq=1E7
    E1=1E7
    eta1=1E14
    EJ=np.ones(nJ)*E1
    etaJ=np.ones(nJ)*eta1
    DrawSheet(Eeq,EJ,etaJ)
    DrawSheetG(EJ,etaJ)


    
    

