# Coded and written by Moritz Hallermann
import os
from tkinter import Button,Tk
import tkinter.font as font
from functools import partial
from importlib import import_module
from PyCSAFT.get_pcpar import get_par
dirs = []

class window():
    def __init__(self,Pure=None,DataDir="SorptionData"):
     self.Pure=Pure
     self.recent="default"
     self.wPolyASD=1
     self.DataDir=DataDir
     self.solvent="default"
     self.molecule="default"
     self.initialpathname=os.path.join(os.getcwd(),self.DataDir)
     self.pathname=self.initialpathname
     self.main()

    def return_several_excelfiles(self,excelfile):
     """If the argument "All excel files" is passed to the function, a dictionary is being
     created which connects chronologically the path of each exclfile to a number starting
     form zero to the respective number of the sum of all excelfiles minus one. This dictionary is returned
     and the window forced to close.
     If just the name of a single excelfile (out of several excelfiles) is passed to the function,
     the path of this specific excelfile is returned and the window forced to close."""

     if excelfile == "All files":
      finalDir = self.pathname#str(os.getcwd())
      dict_sev_excelfiles = {}
      for i,n in enumerate(self.sevExcelfiles):
       dict_sev_excelfiles.update({i:finalDir + "\\" + n})
      self.result=dict_sev_excelfiles
      self.root.destroy()
     else:
      finalDir = self.pathname#str(os.getcwd())
      self.result=str(finalDir + '\\' + excelfile)
      self.root.destroy()
    def click_sev_excelfiles(self):
     """If several excelfiles excist, for each excelfile an individual button is being created.
     Each button is conntected to the function called "return_several_excelfiles" and the button
     passes its name as an argument to the just mentioned function.
     Addtionally a button called "All excel files" is being created which displays all excelfiles
     at once. This button is also connected to the function "return_several_excelfiles" but with the
     argument "All excel files"."""
     self.sevExcelfiles = os.listdir(self.pathname)
     self.excelBtns = [None] * len(self.sevExcelfiles)
     for x in range(len(self.excelBtns)):
      action_with_arg = partial(self.return_several_excelfiles, self.sevExcelfiles[x])
      self.excelBtns[x] = Button(self.root, text=str(self.sevExcelfiles[x]), padx=50, pady=2, font='Helvetica 10 bold', command=action_with_arg)
      self.excelBtns[x].pack(side="top", fill="x")
     action_with_arg_all = partial(self.return_several_excelfiles, "All files")
     self.btnsevEf = Button(self.root, text="All files", padx=50, pady=2, font='Helvetica 10 bold', command=action_with_arg_all)
     self.btnsevEf.pack(side="top", fill="x")


    def return_single_excelfile(self):
     """The path of the selected excelfile is being returned and the window is forced to close."""

     singleExcelfileName = str(self.singleExcelfile[0])
     finalDir = self.pathname #str(os.getcwd())
     self.result=str(finalDir + '\\' + singleExcelfileName)
     #self.result=str(singleExcelfileName)
     #os.chdir(originalpath)
     self.root.destroy()
     #print(str(finalDir + '\\' + singleExcelfileName))



    def click_single_excelfile(self):
     """A button with the name of the single excelfile is being created and connected to the function
     "return_single_excelfile" by the command."""


     self.singleExcelfile = os.listdir(self.pathname)
     self.btnsinEf = Button(self.root, text=str(self.singleExcelfile[0]), padx=50, pady=2, font='Helvetica 10 bold', command=self.return_single_excelfile)
     self.btnsinEf.pack(side="top", fill="x")


    def click_btn(self,dirName):
     """The click_btn function first examines if there already have been buttons created for excelfiles.
     If so, they are being destroyed.

     Afterwards the function checks what the current path is. If the current path prolonged with the
     name of the file or folder (passed argument) is a directory, the path is changed to this directory;
     once again all the files and folders of the subdirectory are saved in a list. If one of the
     files in the subdirectory contains the extension of an excelfile, the currently existing buttons
     are destroyed and for the the function "dynBtngen" is excuted for the new directory.

     In case the subdirectory contains excelfiles, the function examines whether it contains one
     excelfile (all current buttons are destroyed and the function "click_single_excelfile" is
     executed) or several excelfiles (all current buttons are destroyed and the function
     "click_sev_excelfiles" is executed)."""
     #cwd=os.getcwd()
     oldname=self.pathname
     path = oldname#str(os.getcwd())
     self.pathname=os.path.join(oldname,dirName)
     parent=os.path.abspath(os.path.join(oldname, os.pardir))
     if "Pure" in os.listdir(parent):
        self.molecule=dirName

     if self.molecule in str(os.listdir(parent)):
        self.solvent=dirName

     if "ASD" in os.listdir(parent):
        try:
         self.molecule, self.molecule2=dirName.split("_")
        except:
         pass

     self.recent=dirName
     if hasattr(self,'btnsinEf'):
      self.btnsinEf.destroy()
     if hasattr(self,'excelBtns'):
      for x in range(0, len(self.excelBtns)):
       self.excelBtns[x].destroy()
     if hasattr(self,'btnsevEf'):
       self.btnsevEf.destroy()
     if os.path.join(path,dirName):#os.path.isdir(path + '/' + dirName):
      current_subDir = os.listdir(self.pathname)
      extent = ".xlsx"
      extentOld = ".xls"
      extentPXRD=".ras"
      extentDek=".TXT"
      extentDek2=".txt"
      extentSotax=".csv"
      if self.recent.replace(".","").isdigit():
        self.wPolyASD=float(self.recent)
      if (not any(extent in s for s in current_subDir) and not any(extentOld in s for s in current_subDir) and not any(extentPXRD in s for s in current_subDir) and not any(extentDek in s for s in current_subDir) and not any(extentDek2 in s for s in current_subDir)):
       for x in range(0, len(self.btn)):
        self.btn[x].destroy()
       self.currentDir.clear()
       self.dynBtngen()
      elif (len(current_subDir) == 1):
       for x in range(0, len(self.btn)):
        self.btn[x].destroy()
       self.currentDir.clear()

       if self.recent.replace(".","").isdigit():
        self.wPolyASD=float(self.recent)

       self.click_single_excelfile()


      else:
       for x in range(0, len(self.btn)):
        self.btn[x].destroy()
       self.currentDir.clear()

       self.click_sev_excelfiles()


    def dynBtngen(self):
     """The dynBtngen function extracts the current directory and creates a list called "currentDir" with
     all containing folders and files. Then a list called "btn" is created which equales
     the length of previous list "currentDir". In a for loop each entry of the list "btn" is overwritten
     with widget of a button; each with the name of the respective file or folder and the command which
     passes this name as an argument to the function "click_btn"."""


     if self.pathname==self.initialpathname:#os.getcwd()==originalpath+"\\"+self.DataDir:
         self.backBtn.configure(state="disabled")
     else:
         self.backBtn.configure(state="normal")
     self.currentDir = os.listdir(self.pathname)


     self.btn = [None] * len(self.currentDir)

     # print(self.pathname)
     for x in range(0, len(self.btn)):
      action_with_arg = partial(self.click_btn, self.currentDir[x])
      self.btn[x] = Button(self.root, text=str(self.currentDir[x]), padx=50, pady=2, font='Helvetica 10 bold', command=action_with_arg)
      self.btn[x].pack(side="top", fill="x")
      if "Pure" in self.currentDir and self.Pure==True:
          if str(self.currentDir[x])=="ASD":
           self.btn[x].configure(state="disabled")
      elif "Pure" in self.currentDir and self.Pure==False:
          if str(self.currentDir[x])=="Pure":
           self.btn[x].configure(state="disabled")
      elif str(self.currentDir[x])=="Notselected":
           self.btn[x].destroy()



    def main(self):
     """The main function initializes the GUI-window, is running all the time and contains basic
     design settings. The intial path of the starting directory has to be typed in. Since
     the "Go back"-Button should always appear, it is part of the main function and is connected via the
     command to the partial method. If the "Go back"-Button is clicked, it passes the argument ".." to
     function "click_btn" and executes it. Eventually the function "dynBtngen" is executed."""

     self.root = Tk()
     self.root.attributes("-topmost", True)

     self.root.title("Project Maxwell v0.3")
     #root.iconbitmap("TU.ico")#here the TU Dortmund Icon has be in the same folder as the python script and
                              #the image as to be an icon type.

     #os.chdir(originalpath+'\\'+self.DataDir)

     action_with_arg = partial(self.click_btn, "..")
     self.backBtn = Button(self.root, text="Go back", padx=10, pady=3, font='Helvetica 8 italic', command=action_with_arg)
     self.backBtn.pack(anchor = "w", side = "bottom")

     self.dynBtngen()

     self.root.mainloop()
     #os.chdir(originalpath)
     if self.DataDir!="PXRD":
         pure,kij,ksw=get_par([self.molecule,self.solvent,self.molecule2]) if hasattr(self,"molecule2") else get_par([self.molecule,self.solvent,self.molecule])

         pol,water,api=pure

         self.components=(pol,water,api) # if hasattr(self,"molecule2") else (pol,water)# Import must be this way! 29.06.20

if __name__ == '__main__':

 window1=window(Pure=False,DataDir="SorptionData")
 #window1=window(Pure=False,DataDir="DissolutionData")
