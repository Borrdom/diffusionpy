#loop over all csv and xlsx files in the current directory
import os
import pandas as pd


def changeheader(filename,newheader,id):
    oldheader=pd.read_excel(filename).columns.tolist() if id else pd.read_csv(filename).columns.tolist()
    frame=pd.read_excel(filename) if id else pd.read_csv(filename)
    while len(newheader)<len(oldheader):
        newheader.append("dummy")
    frame.columns=newheader
    frame.to_excel(filename) if id else pd.to_csv(filename)
    print(oldheader)
    print(newheader)
    return

if __name__=="__main__":
    newheader=["w2/- 0",
    "w2/- 1",
    "w2/- 2",
    "w2ave/- 0",
    "w2std/- 0",
    "RH/- 0",
    "RH/- 1",
    "RH/- 2",
    "RHave/- 0",
    "RHstd/- 0",
    "w2bar/- 0",
    "w2bar/- 1",
    "w2bar/- 2",
    "w2barave/- 0",
    "w2barstd/- 0",
    "DCranc/m/s^2 0",
    "DCranc/m/s^2 1",
    "DCranc/m/s^2 2",
    "DCrancave/m/s^2 0",
    "DCrancstd/m/s^2 0",
    "DFick/m/s^2 0",
    "DFick/m/s^2 1",
    "DFick/m/s^2 2",
    "DFickave/m/s^2 0",
    "DFickstd/m/s^2 0",
    "DStefan/m/s^2 0",
    "DStefan/m/s^2 1",
    "DStefan/m/s^2 2",
    "DStefanave/m/s^2 0",
    "DFStefanstd/m/s^2 0"]
    files=os.listdir()
    filexlsx=[val.endswith(".xlsx") for i,val in enumerate(files)]
    filecsv=[val.endswith(".csv") for i,val in enumerate(files)]
    [changeheader(val,newheader,1)for i,val in enumerate(files) if filexlsx[i] ]
    #[changeheader(val,newheader,0)   for i,val in enumerate(files) if filecsv[i]]

#read all headers in these files and replace them with another header
