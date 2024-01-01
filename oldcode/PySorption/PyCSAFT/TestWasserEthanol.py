from PyCSAFT.PyCSAFTNil import Mixture
from PyCSAFT import get_pcpar
import numpy as np
T=298.15
pure, kij,ksw=get_pcpar.get_par(["water","ethanol"],T=T)
#pure, kij,ksw=get_pcpar.get_par(["hexane","cyclohexane"],T=T)
water,ethanol=pure
#water["sigi"]=2.7927
# ethanol["Mi"]=46.069/1000
# ethanol["mi"]=2.382697894
# ethanol["sigi"]=3.1771
# ethanol["ui"]=198.24
# ethanol["epsAiBi"]=2653.4
# ethanol["kapi"]=0.032384


water["Mi"]=18.015/1000
water["mi"]=1.20469
#water["sigi"]=2.7927
water["ui"]=353.95
water["epsAiBi"]=2425.67
water["kapi"]=0.04509
#water["sigi"]=2.7927
kij['waterethanol']=0
kij['ethanolwater']=0

kij['hexanecyclohexane']=0
kij['cyclohexanehexane']=0
mix=Mixture(water,ethanol,dikij=kij)
xi=np.asarray([[0,1],[0.25,0.75],[0.5,0.5],[0.75,0.25],[1,0]])

#xi=np.asarray([[0,1],[0.1,0.9],[0.2,0.8],[0.3,0.7],[0.4,0.6],[0.5,0.5],[0.6,0.4],[0.7,0.3],[0.8,0.2],[0.9,0.1],[1,0]])
results=[mix.PTFlash(psys=1E5,T=298.15,liq=1,xi=xi[i,:]) for i,val in enumerate (xi[:,0])]
mix.WriteToExcel([results],["PCSAFT"],Datatype="PTFlash")