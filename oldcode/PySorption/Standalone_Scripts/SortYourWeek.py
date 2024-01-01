import time, datetime, os
from os import listdir
from os.path import isfile, join
import os
import shutil


def sort_files_in_a_folder(mypath,weekpath):
    '''
    A function to sort the files in a download folder
    into their respective categories
    '''
    files = [f for f in listdir(mypath) if isfile(join(mypath, f))]   
    file_type_variation_list=[]
    filetype_folder_dict={}
    for file in files:
        filetype=file.split('.')[1]
        if filetype not in file_type_variation_list:
            file_type_variation_list.append(filetype)
            new_folder_name=mypath+'/'+ filetype + '_folder'
            filetype_folder_dict[str(filetype)]=str(new_folder_name)
            if os.path.isdir(new_folder_name)==True:  #folder exists
                continue
            else:
                os.mkdir(new_folder_name)
    for file in files:
        src_path=mypath+'/'+file
        filetype=file.split('.')[1]
        if filetype in filetype_folder_dict.keys():
            dest_path=filetype_folder_dict[str(filetype)]
            shutil.move(src_path,dest_path)
        print(src_path + '>>>' + dest_path)
        
        #newloc=dest_path+'/'+file
        #desktop = weekpath
        target = os.path.join(dest_path,file)
        path =target+'.lnk'
        shell = win32com.client.Dispatch("WScript.Shell")
        shortcut = shell.CreateShortCut(path)
        shortcut.Targetpath = target
        #shortcut.IconLocation = icon
        shortcut.WindowStyle = 7 # 7 - Minimized, 3 - Maximized, 1 - Normal
        shortcut.save()
        shutil.move(path,weekpath)

#path=r'C:\Users\domin\OneDrive\Promotion\Wissen'
from os import getcwd
path=getcwd().replace("ProgrammCode\\Diffusionsmodell\\Python_Boettcher","Wissen")
    
today = datetime.date.today()  
year,month,day=tuple(str(today).split("-"))

week=str(today.isocalendar()[1])

stri="\Woche"+week+"_"+year

try:
    os.mkdir(path+stri)
except:
    print(path+stri+" is already existing")
pathweek=path+stri

import win32com.client
import pythoncom
import os


#path2=r'C:\Users\domin\OneDrive\Promotion'

path2=getcwd().replace("ProgrammCode\\Diffusionsmodell\\Python_Boettcher","")
    

print(os.listdir(path2))
fold=os.listdir(path2)
for i, val in enumerate(fold):
    dest_path=path2
    try:
        os.mkdir(pathweek+'\\'+val)
    except:
        pass
    mypath=path2+'\\'+val
    sort_files_in_a_folder(mypath,pathweek+'\\'+val)

# pythoncom.CoInitialize() # remove the '#' at the beginning of the line if running in a thread.
#desktop = r'C:\Users\domin\Desktop' # path to where you want to put the .lnk
#path = os.path.join(desktop, 'Neues Textdokument.lnk')
#target = r'C:\Users\domin\Desktop\Neues Textdokument.txt'
#icon = r'C:\path\to\icon\resource.ico' # not needed, but nice

#shell = win32com.client.Dispatch("WScript.Shell")
#shortcut = shell.CreateShortCut(path)
#shortcut.Targetpath = target
#shortcut.IconLocation = icon
#shortcut.WindowStyle = 7 # 7 - Minimized, 3 - Maximized, 1 - Normal
#shortcut.save()


