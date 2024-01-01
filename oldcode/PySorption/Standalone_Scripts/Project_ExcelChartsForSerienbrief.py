# import pandas as pd
# from tkinter.filedialog  import askopenfilename
# #sheet_name=None

# writer = pd.ExcelWriter('Mappe1.xlsx', engine='xlsxwriter')
# workbook = writer.book
# worksheet = writer.sheets['Tabelle1']
# charts={}
# for i,val in enumerate(sheets):
#     charts[val]=sheets[val].add_chart({'type': 'column'})


import excel2img
import os
cwd=os.getcwd()
from xlrd import open_workbook
import xlsxwriter
#hidden imports
import numpy
import pkg_resources.py2_warn
import multiprocessing
import sys
inputfile="NamenUndZellen.xlsx"
inputsname=os.path.join(cwd,inputfile)
#saveranges=pd.read_excel(inputsname)
saveranges=open_workbook(inputsname).sheet_by_index(0)
names=[val.value for i,val in enumerate(saveranges.col(0))][1:]
ranges=[val.value for i,val in enumerate(saveranges.col(1))][1:]
filenames=[val.value for i,val in enumerate(saveranges.col(2)) if len(val.value)>0][1:]


outname=os.path.join(cwd,"Diagramme")
os.makedirs(outname) if not os.path.exists(outname) else None
auswertung="Auswertung_Diagrammausgabe"
print("Hallo! Bitte beachten Sie, dass Ihre Eingaben aus einer Datei "+inputfile+" ausgelesen werden.\n")
print("Diese muss sich im selben Ordner wie diese .EXE befinden.\n")
print("Die ausgelesen Dateien müssen sich in einen Unterordner dieser .EXE mit den Ordnername <<"+auswertung+">> befinden.\n")
#print("Zudem werden nur .xlsx Dateien unterstützt\n")
print("Die Zellen können entweder einzelnd: sprich <<A1>> oder als Reihe angeben werden <<A1:C3>>\n")
print("Ich habe folgende Eingaben von der Datei <<"+inputsname+">> erfasst:\n")
print("Namen:\t Zellen:\n")
[print(names[i]+"\t"+ranges[i]+"\n") for i,val in enumerate(names)]
print("Und die folgenden Dateien sollen ausgelesen werden:\n")
[print(filenames[i]+"\n") for i,val in enumerate(filenames)]
None if input("Waren das Ihre Eingaben? Bestätigen mit <<ja>>.\t")=="ja" else sys.exit(0)
print("\nBitte warten!")
file_dict={"Dateiname":[],"Tabellenname":[],"Bildname":[],"Gesamtname":[]}
def Serienbrief(filename,names,ranges):
    pathname=os.path.join(cwd,auswertung)
    excelname1=os.path.join(pathname,filename+".xlsx")
    excelname2=os.path.join(pathname,filename+".xls")
    
    if os.path.isfile(excelname1):
        fullname=excelname1
    elif os.path.isfile(excelname2):
        fullname=excelname2
    else:
        print("Warnung: Datei"+excelname1+"nicht gefunden.")
    
    #book=pd.read_excel(fullname,sheet_name=None,header=None)
    book=open_workbook(fullname)
    sheetnames=book.sheet_names()
    def SaveChart(filename,sheetname,saveranges):
        for j,valj in enumerate(names):
            imagename=os.path.join(outname,"_".join([filename,sheetname,names[j]])+".png")
            excel2img.export_img(fullname, imagename, sheetname, ranges[j])
            file_dict["Dateiname"].append(filename)
            file_dict["Tabellenname"].append(sheetname)
            file_dict["Bildname"].append(names[j])
            file_dict["Gesamtname"].append("_".join([filename,sheetname,names[j]]))

    for i,val in enumerate(sheetnames):
        SaveChart(filename,val,saveranges)
        
for k,valk in enumerate(filenames):
    Serienbrief(valk,names,ranges)





workbook = xlsxwriter.Workbook(os.path.join(outname,"Historie.xlsx"))
worksheet = workbook.add_worksheet()

my_dict = file_dict
col_num = 0
for key, value in my_dict.items():
    worksheet.write(0, col_num, key)
    worksheet.write_column(1, col_num, value)
    col_num += 1

workbook.close()
print("\nErfolg! Dateien im Pfad "+outname+" als .png Dateien abgespeichert und Historie.xlsx erstellt.\n")
print("Vielen Dank der Nutzung dieser Software. Es gilt:\n<<THE BEER-WARE LICENSE>> (Revision 42) <<dominik.borrmann@tu-dortmund.de>> schrieb dieses Programm. Solange Sie diesen Vermerk nicht entfernen, können Sie mit der Datei machen, was Sie möchten. Wenn wir uns eines Tages treffen und Sie denken, die Datei ist es wert, können Sie mir dafür ein Bier ausgeben.\nDominik Borrmann")
input("Drücken Sie eine beliebige Taste zum Schließen dieses Fensters.")

