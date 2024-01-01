# This is a simple example how to execute the project properly
from PySorption.DiffusionFinal import Film # Film is an object which inherits method for the diffusion modelling
from PySorption.dcc import window # dcc is a graphical interface to browse through to your Sorptiondata and automatically detects the necessary components and drugload
window1=window(DataDir="SorptionData")  # within the dcc object is the get_par class which extracts the pure component parameters for each components
Film1=Film(window1.result,window1.components,window1.wPolyASD) # Film
Film1.SAFT
Film1.CalcTerDiffcoeff()
