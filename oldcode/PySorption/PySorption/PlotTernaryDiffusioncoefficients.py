from PySorption.DiffusionFinal import Film
from PySorption.dcc import window
import matplotlib.pyplot as plt
window1=window(DataDir="SorptionData")
Film1=Film(window1.result,window1.components,window1.wPolyASD)
Film1.SAFT=True
Film1.NET=True
Film1.Origin=True
Film1.CalcTerDiffcoeff(Cranc=False,Plot=True)
plt.show()
