import pandas as pd

from PyCSAFTNil import Mixture

#from PureDataBase import naproxen,acetone,pvp,pvpva
#import PureDataBase
import matplotlib.pyplot as plt
import numpy as np
from lmfit import Parameters,minimize
import os
from move_figure import move_figure
from get_pcpar import get_par

#Solvents=["ethanol","acetonitrile","ethylacetate","toluene","hexane","cyclohexane","carbontatrachloride","chloroform","tetrahydrofuran","dichloromethane","methylisobutylketone","butanol","isopropylalcohol","acetone","methanol","dimethylsulfoxide","octanol"]

#Solvents=["acetone","acetonitrile","ethanol","ethylacetate","water"]


Polymers=["pvpva","hpmcas"]
apiname="naproxen"
mode="Exp_Lit"
path=os.path.join(os.getcwd(),"Projekt_Kruemeltee",mode)
files=os.listdir(path)
Solvents=[val.replace("_celecoxib.xlsx","") for i,val in enumerate(files)]
#for i,val in enumerate(Polymers):
#    Solvents=np.asarray(Solvents)[np.asarray(Solvents)!=val].tolist()

#chloroform was not in sugar
#methanol has a big kij was in sugar but not anymore
#["methanol","toluene","hexane","cyclohexane","carbontetrachloride","acetone","ethylacetate","methylisobutylketone","tetrahydrofuran","dimethylsulfoxide","pvpva","hpmcas"]
#All removing that are associating to get the 3 Parameters right["acetonitrile","carbontetrachloride","tetrahydrofuran","methylisobutylketone","tetrahydrofuran","acetone","dimethylsulfoxide","ethylacetate","methanol","ethanol","butanol","isopropylalcohol,"octanol","water","pvpva","hpmcas"]

#Removing all indued association and association
Ausnahmen=["acetonitrile","carbontetrachloride","tetrahydrofuran","methylisobutylketone","tetrahydrofuran","acetone","dimethylsulfoxide","ethylacetate","methanol","ethanol","butanol","isopropylalcohol","octanol","water","pvpva","hpmcas"]

#Get to know why kijs from hexane toluol and cyclohexane are so big
#Ausnahmen=["toluene","dichloromethane","chloroform","cyclohexane","acetonitrile","carbontetrachloride","tetrahydrofuran","methylisobutylketone","tetrahydrofuran","acetone","dimethylsulfoxide","ethylacetate","methanol","ethanol","butanol","isopropylalcohol","octanol","water","pvpva","hpmcas"]


#Removing all induced-assocs and non assocs

#Ausnahmen=["toluene","hexane","cyclohexane","dichloromethane","chloroform","acetonitrile","carbontetrachloride","tetrahydrofuran","methylisobutylketone","tetrahydrofuran","acetone","dimethylsulfoxide","ethylacetate","pvpva","hpmcas"]

#Ethanol und Isopropyl machen Nan
#Ausnahmen=["isopropylalcohol","water","octanol","dichloromethane","butanol","methanol","toluene","hexane","cyclohexane","dichlormethane","chloroform","acetonitrile","carbontetrachloride","tetrahydrofuran","methylisobutylketone","tetrahydrofuran","acetone","dimethylsulfoxide","ethylacetate","pvpva","hpmcas"]



#Removing all assocs
#Ausnahmen=["toluene","hexane","cyclohexane","dichlormethane","chloroform","methanol","ethanol","butanol","isopropylalcohol","octanol","water","pvpva","hpmcas"]

for i,val in enumerate(Ausnahmen):
    Solvents=np.asarray(Solvents)[np.asarray(Solvents)!=val].tolist()

class SolubilityFitting:
    def __init__(self,apiname,solventnames):
        self.apiname=apiname
        self.solventnames=solventnames
    def BinarySolubilities(self, par,solventname):
        path=os.path.join(os.getcwd(),"Projekt_Kruemeltee",mode,solventname+"_"+self.apiname+".xlsx")
        data=pd.read_excel(path)
        T=data["T[°C]"].dropna().values+273.15
        w=data["w[-]"].dropna().values
        #solvent=getattr(PureDataBase,solventname)
        solvent=get_par([solventname])[0][0]
        api=get_par([self.apiname])[0][0]
        dummy=get_par(["pvp"])[0][0]
        #api=getattr(PureDataBase,self.apiname)
        #dummy=getattr(PureDataBase,"pvp")
        api["mi"]=par["miMi"].value*api["Mi"]*1000/100
        api["sigi"]=par["sigi"].value
        api["ui"]=par["ui"].value*100
        api["epsAiBi"]=par["epsAiBi"].value*1000
        dikij={}
        dikijb={}
        dikij[apiname+solventname]=par["kija"+solventname].value
        dikijb[apiname+solventname]=par["kijb"+solventname].value
        print(par["kija"+solventname].value)
        #dikij[apiname+solventname]=par["kija"].value
        #dikijb[apiname+solventname]=par["kijb"].value
        BinarySLE=Mixture(solvent,dummy,api,wPolyASD=0,dikij=dikij,dikijb=dikijb)
        wsim=[]
        BinarySLE.PCSAFT(psys=1.013E5,T=298.15)
        rhos=BinarySLE.Pure(psys=1.013E5,liq=True,com=2,T=298.15)["rho"]
        self.rhos=rhos
        par["rhos"].value=rhos
        print(rhos)
        print(solventname)
        if rhos>500 and rhos<2000:
            for i,val in enumerate(T):
                dictL=BinarySLE.SLE(psys=1.013E5,wwASD=0,T=val)
                wsim.append(dictL["wi"][2])
                #print(dictL["rho"])
            wsim=np.asarray(wsim)
        else:
            wsim=np.ones_like(w) #better not go for zero because low soluble substance


        #fig,ax=plt.subplots(figsize=(2,2))
        #ax.plot(T,wsim,'ko')
        #ax.plot(T,w,'kx')
        #ax.set_title(solventname)
        #self.j=0
        #self.j+=1
        #move_figure(fig, 0+200*self.j, 200)
        #plt.pause(.001)

        print((wsim-w)/w) # Maybe relative error :D
        return (wsim-w)/w
    def SimForPlot(self, par,solventname):
        path=os.path.join(os.getcwd(),"Projekt_Kruemeltee",mode,solventname+"_"+self.apiname+".xlsx")
        data=pd.read_excel(path)
        T=data["T[°C]"].dropna().values+273.15
        w=data["w[-]"].dropna().values
        # solvent=getattr(PureDataBase,solventname)
        # api=getattr(PureDataBase,self.apiname)
        # dummy=getattr(PureDataBase,"pvp")
        solvent=get_par([solventname])[0][0]
        api=get_par([self.apiname])[0][0]
        dummy=get_par(["pvp"])[0][0]
        api["mi"]=par["miMi"].value*api["Mi"]*1000/100
        api["sigi"]=par["sigi"].value
        api["ui"]=par["ui"].value*100
        api["epsAiBi"]=par["epsAiBi"].value*1000
        dikij={}
        dikijb={}
        dikij[apiname+solventname]=par["kija"+solventname].value
        dikijb[apiname+solventname]=par["kijb"+solventname].value

        #dikij[apiname+solventname]=par["kija"].value
        #dikijb[apiname+solventname]=par["kijb"].value

        BinarySLE=Mixture(solvent,dummy,api,wPolyASD=0,dikij=dikij,dikijb=dikijb)
        wsim=[]
        BinarySLE.PCSAFT(psys=1.013E5,T=298.15)
        #rhos=BinarySLE.Pure(psys=1.013E5,liq=True,com=2,T=298.15)["rho"]
        #print(rhos)
        #print(solventname)

        for i,val in enumerate(T):
            dictL=BinarySLE.SLE(psys=1.013E5,wwASD=0,T=val,i=2)
            wsim.append(dictL["wi"][2])
                #print(dictL["rho"])
        wsim=np.asarray(wsim)

        return wsim,w,T
    def Plot(self, par):
        wsim,w,T=[],[],[]
        for i,val in enumerate(self.solventnames):
            t1,t2,t3=self.SimForPlot(par,val)
            wsim.append(t1)
            w.append(t2)
            T.append(t3)
        import OriginVorlagePlot
        xe,ye=[np.zeros_like(val) for i,val in enumerate(T)],[np.zeros_like(val) for i,val in enumerate(T)]
        x,xs,y,ys=T,T,w,wsim
        xData,yData,xSim,ySim,xErr,yErr=OriginVorlagePlot.MusterDictOrigin()
        xlab="T[°C]"
        ylab="w[-]"
        xu="°C"
        yu="-"
        xc="temperature"
        yc="solubility"
        name="Solubilities_"+mode
        xData["values"],yData["values"],xSim["values"],ySim["values"],xErr["values"],yErr["values"]=x,y,xs,ys,xe,ye
        xData["unit"],yData["unit"],xSim["unit"],ySim["unit"],xErr["unit"],yErr["unit"]=xu,yu,xu,yu,xu,yu
        xData["comment"],yData["comment"],xSim["comment"],ySim["comment"],xErr["comment"],yErr["comment"]=xc,yc,xc,yc,xc,yc

        OriginVorlagePlot.Plot(xData,yData,xSim,ySim,xErr,yErr,ylab,xlab,name)
        return

    def Solubilities(self,par):
        binarysolubilities=np.asarray([])
        for i,val in enumerate(self.solventnames):
            par["kija"+val]=self.par["kija"+val]
            par["kijb"+val]=self.par["kijb"+val]
            opt=self.kijopt(par,val)
            par["kija"+val]=opt.params["kija"+val]
            par["kijb"+val]=opt.params["kijb"+val]
            self.par["kija"+val]=opt.params["kija"+val]
            self.par["kijb"+val]=opt.params["kijb"+val]
            #binarysolubilities=np.hstack((binarysolubilities,self.BinarySolubilities(par,val),opt.params["kija"+val]))

            binarysolubilities=np.hstack((binarysolubilities,self.BinarySolubilities(par,val)))

        #plt.close("all")
        print(par)
        print(binarysolubilities)
        #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4069839/ True Density of CBX 1350 kg/m^3
        binarysolubilities=np.hstack((binarysolubilities,(self.rhos-1350)/10000))
        return binarysolubilities

    def locopt(self,par):
        par["miMi"].vary=0
        par["sigi"].vary=0
        par["ui"].vary=0
        par["epsAiBi"].vary=0
        self.par=par
        opt=minimize(self.Solubilities,par,method="powell",tol=1E-3)
        #fitter=Minimizer()
        #fitter._asteval.symtable['lorentzian']
        return opt

    def locoptpoly(self,par):
        par["miMi"].vary=0
        par["sigi"].vary=0
        par["ui"].vary=0
        par["epsAiBi"].vary=0
        self.par=par
        opt=minimize(self.Solubilities,par,method="cobyla",tol=1E-3)

        return opt
    def kijopt(self,par,solventname):
        for i,val in enumerate(par):
            par[val].vary=0
        par["kija"+solventname].vary=1

        par["kijb"+solventname].vary=0
        opt=minimize(lambda par:self.BinarySolubilities(par,solventname),par,method="least_squares",ftol=1E-3)

        return opt
    def globopt(self,par):
        opt=locopt(par)
        opt["miMi"].vary=0
        opt["sigi"].vary=0
        opt["ui"].vary=1
        return self.BinarySolubilities(opt,self.polyname,self.apiname)

Mixture1=SolubilityFitting(apiname,Solvents)
par=Parameters()
par.add('rhos', value=0, vary=0)
#par.add('bound', expr="rhos-1300")

#par.add('ui', value=2.082, min=1.5, max=4,vary=1)
#par.add('miMi', value=3.3, min=0.5, max=5,vary=1)
#par.add('sigi', value=2.7, min=1, max=4,vary=1)
#par.add('epsAiBi', value=1.596, min=0, max=5,vary=1)

#par.add('ui', value=1.9738722840585259, min=1.5, max=4,vary=1)
#par.add('miMi', value=4.297955891472619, min=0.5, max=5,vary=1) #Dichte 1300
#par.add('sigi', value=2.880674080809237, min=1, max=4,vary=1)
#par.add('epsAiBi', value=0.9802773156864133, min=0, max=5,vary=1)
#par.add('kija', value=0, min=-.1, max=.1,vary=0)
#par.add('kijb', value=0, min=-.1, max=.1,vary=0)

# par.add('ui', value=2.28336852082822, min=1.5, max=4,vary=1)
# par.add('miMi', value=1.9119014569237, min=0.5, max=5,vary=1) #Dichte 2700
# par.add('sigi', value=2.93408444810157, min=1, max=4,vary=1)
# par.add('epsAiBi', value=0.271345957583785, min=0, max=5,vary=1)


# par.add('ui', value=3.05064056, min=1.5, max=4,vary=1)
# par.add('miMi', value=1.625101394, min=0.5, max=5,vary=1)
# par.add('sigi', value=3.305618723, min=1, max=4,vary=1) #Dichte 2400 mit an Polymere angepasst
# par.add('epsAiBi', value=0.64156794, min=0, max=5,vary=1)


# par.add('ui', value=, min=1.5, max=4,vary=1)
# par.add('miMi', value=, min=0.5, max=5,vary=1)
# par.add('sigi', value=, min=1, max=4,vary=1)   #Dichte 3000 nich an Polymere angepasst
# par.add('epsAiBi', value=, min=0, max=5,vary=1)



#par.add('ui', value=3.98983939769083, min=1.5, max=4.5,vary=1)
#par.add('miMi', value=2.60526461959088, min=0.5, max=5,vary=1)
#par.add('sigi', value=3.46389488119008, min=1, max=4,vary=1)   #Dichte 1461 nicht an Polymere angepasst Dichte Penalty
#par.add('epsAiBi', value=0.126399501110963, min=0, max=5,vary=0)

par.add('ui', value=4.05747230664321, min=1.5, max=4.5,vary=1)
par.add('miMi', value=2.89922993078485, min=0.5, max=5,vary=1)
par.add('sigi', value=3.4268214389218, min=1, max=4,vary=1)   #Dichte 1461 nicht an Polymere angepasst Dichte Penalty
par.add('epsAiBi', value=0.824, min=0, max=5,vary=0)


#2.7 3.3
#kijadict=[0.0788,-0.0937,-0.086,-0.0989,0.07834]
#kijbdict=[0.0834,-0.071,-0.01016,-0.0558,-0.012]
#kijadict=[0,0,0,0,0]
#kijadict=[-0.02,-0.044,-0.049,-0.023,-0.0275]
kijadict=[0.012669995071665593,-0.005214594521115445,-0.0004029752058303526]
kijbdict=[0,0,0,0]
kijadict=[0 for i,val in enumerate(Solvents)]
kijbdict=[0 for i,val in enumerate(Solvents)]
[par.add('kija'+val, value=kijadict[i], min=-.2, max=.2,vary=0) for i,val in enumerate(Solvents)]
[par.add('kijb'+val, value=kijbdict[i], min=-.2, max=.2,vary=0) for i,val in enumerate(Solvents)]
#kij to sensitive?
#maybe additional fiting inside for kijs
#acetonitril very good
#now it varies also the values
#ui 200.58
#miMI 0.01
#2.414
#epsAiBi 600

#176
#0.0076
#2.56
#288.7
#[-0.02,-0.0044,-0.0049,-0.023,-0.0275]
#
#With density constraint!
#kijb and kija not working
#Parameters([('ui', <Parameter 'ui', value=1.7653660471622383, bounds=[1.5:4]>), ('miMi', <Parameter 'miMi', value=0.7322715602946933, bounds=[0.5:5]>), ('sigi', <Parameter 'sigi', value=2.609381600187967, bounds=[1:4]>), ('epsAiBi', <Parameter 'epsAiBi', value=2.730611329488338, bounds=[0:15]>), ('kijaacetone', <Parameter 'kijaacetone', value=0.09446964629989707, bounds=[-0.1:0.1]>), ('kijaacetonitrile', <Parameter 'kijaacetonitrile', value=-0.008035947457750242, bounds=[-0.1:0.1]>), ('kijaethanol', <Parameter 'kijaethanol', value=0.004198017745413091, bounds=[-0.1:0.1]>), ('kijaethylacetate', <Parameter 'kijaethylacetate', value=-0.02524816256284415, bounds=[-0.1:0.1]>), ('kijawater', <Parameter 'kijawater', value=-0.041135622514311715, bounds=[-0.1:0.1]>), ('kijbacetone', <Parameter 'kijbacetone', value=0 (fixed), bounds=[-0.1:0.1]>), ('kijbacetonitrile', <Parameter 'kijbacetonitrile', value=0 (fixed), bounds=[-0.1:0.1]>), ('kijbethanol', <Parameter 'kijbethanol', value=0 (fixed), bounds=[-0.1:0.1]>), ('kijbethylacetate', <Parameter 'kijbethylacetate', value=0 (fixed), bounds=[-0.1:0.1]>), ('kijbwater', <Parameter 'kijbwater', value=0 (fixed), bounds=[-0.1:0.1]>)])
#How will i know which kij must be a vector dependent on solvents
opt=Mixture1.locopt(par)

Mixture1.Plot(opt.params)
#Mixture1.Plot(par)
#Mixture1.par["miMi"]=opt.params["miMi"]/100
#Mixture1.par["sigi"]=opt.params["sigi"]
#Mixture1.par["ui"]=opt.params["ui"]*100
#Mixture1.par["epsAiBi"]=opt.params["epsAiBi"]*1000

#Succesfull optimization

#1.9738722840585259
#4.297955891472619
#2.880674080809237
#0.9802773156864133

#0.026961792850123956
#0.012669995071665593
#-0.005214594521115445
#-0.0004029752058303526
#-0.09996967584162647

# Mixture1.par=par
# Mixture1.Solubilities(par)
# opt=type('test', (object,), {})()
# opt.params=Parameters()
# opt.params.add('ui', value=1.9738722840585259, min=1.5, max=4,vary=1)
# opt.params.add('miMi', value=4.297955891472619, min=0.5, max=5,vary=1)
# opt.params.add('sigi', value=2.880674080809237, min=1, max=4,vary=1)
# opt.params.add('epsAiBi', value=0.9802773156864133, min=0, max=5,vary=1)
# [opt.params.add('kija'+val, value=kijadict[i], min=-.1, max=.1,vary=1) for i,val in enumerate(Solvents)]
# [opt.params.add('kijb'+val, value=kijbdict[i], min=-.1, max=.1,vary=0) for i,val in enumerate(Solvents)]
for i,val in enumerate(Solvents):
    opt.params['kija'+val].value=Mixture1.par["kija"+val].value
Mixture2=SolubilityFitting("celecoxib",Polymers)
[opt.params.add('kija'+val, value=0, min=-.2, max=.2,vary=0) for i,val in enumerate(Polymers)]
[opt.params.add('kijb'+val, value=0, min=-.2, max=.2,vary=0) for i,val in enumerate(Polymers)]
opt=Mixture2.locoptpoly(opt.params)
# opt.params.add('kija'+"pvpva", value=-0.052, min=-.1, max=.1,vary=1)
# opt.params.add('kijb'+"pvpva", value=0, min=-.1, max=.1,vary=0)
# opt.params.add('kija'+"hpmcas", value=-0.052, min=-.1, max=.1,vary=1)
# opt.params.add('kijb'+"hpmcas", value=0, min=-.1, max=.1,vary=0)
import pybroom as br
import pandas as pd
dt = br.tidy(opt)
#[dt['kija'+val]=opt.params["kija"+val].value for i,val in enumerate(Solvents)]
dt.to_excel(apiname+mode+"_fitting.xlsx")
#Best results for acetone and acetonnitril not asssociating but ethylacetate induced assocuating
#Good density for acetone and acetonnitril ethylacetat not asssociating


desired=Mixture1.SimForPlot(opt.params,"pvpva")
desired=Mixture1.SimForPlot(opt.params,"hpmcas")
#Mixture2.Plot(opt.params)
#Solubilities(par)
#opt=minimize(globopt,par,method="cobyla")
