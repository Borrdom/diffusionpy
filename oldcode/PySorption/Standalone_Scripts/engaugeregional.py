from tkinter.filedialog import askopenfilename

filename = askopenfilename()
with open (filename, "r") as myfile:
    file=myfile.readlines()
file=[val.replace(',','.') for i,val in enumerate(file)]
#file=[val.replace(';',',') for i,val in enumerate(file)]
myfile.close()
with open (filename, "w") as myfile:
    [myfile.write(val) for i,val in enumerate(file)]
myfile.close()