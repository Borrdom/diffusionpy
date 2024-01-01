import pandas as pd
#from tkinter.filedialog import askopenfilename
import os
import pandas.io.formats.excel
pandas.io.formats.excel.header_style = None
import numpy
import pkg_resources.py2_warn
import multiprocessing
#original=askopenfilename()

#pathna=os.path.join(os.getcwd(),"ToBeReduced")
pathna=os.getcwd()
#originallist=os.listdir(pathna)
originallist=[f for f in os.listdir(pathna) if (f.endswith('.xlsx') or f.endswith('.xls'))]
print("Folgende Excel-Dateien gefunden:\n")
for i,val in enumerate(originallist):
    print(str(i)+": "+val)
userchoice=input("Welche Datei? Gebe die entsprechende Nummer ein\n")
original=os.path.join(pathna,originallist[int(userchoice)])
Data=pd.read_excel(original)

print("Ich habe folgende Eingaben von des Phades <<"+original+">> erfasst:\n")
Datanumberheader=Data[Data['DVS-Advantage-Plus-Data-File'].str.contains("minutes").fillna(False)]
idx=Datanumberheader.index[0]
spacing=10
Dataupper=Data.loc[:idx-1]
Datared=Data.loc[idx::spacing]
DataNew=pd.concat([Dataupper,Datared])
#DataNew=pd.concat([pd.DataFrame(DataNew.columns),DataNew],ignore_index=True,axis=1)
filename=original.split(".")[-1]
newfilename=original.split(".")[0]+"_red."+filename
DataNew.to_excel(newfilename,index=False)

print("\nErfolg! Dateien im Pfad "+newfilename+" gespeichert.\n")
print("Vielen Dank der Nutzung dieser Software. Es gilt:\n<<THE BEER-WARE LICENSE>> (Revision 42) <<dominik.borrmann@tu-dortmund.de>> schrieb dieses Programm. Solange Sie diesen Vermerk nicht entfernen, können Sie mit der Datei machen, was Sie möchten. Wenn wir uns eines Tages treffen und Sie denken, die Datei ist es wert, können Sie mir dafür ein Bier ausgeben.\nDominik Borrmann")
input("Drücken Sie eine beliebige Taste zum Schließen dieses Fensters.")