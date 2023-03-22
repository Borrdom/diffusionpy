import tkinter as tk
from tkinter.filedialog import askopenfilename
from ttkwidgets.autocomplete import AutocompleteEntryListbox
import pyodbc
import pandas as pd
conn_str = (r'DRIVER={Microsoft Access Driver (*.mdb, *.accdb)};'r'DBQ='+askopenfilename())
conn = pyodbc.connect(conn_str)
sheets=[row.table_name for row in conn.cursor().tables()]
unary=sheets[-2]
binary=sheets[-1]
query_unary = "SELECT * from {}".format(unary)
query_binary = "SELECT * from {}".format(binary)
df1 = pd.read_sql(query_unary, conn)
df2 = pd.read_sql(query_binary, conn)
window = tk.Tk()
Label1=tk.Label(window, text="Nothing")
Label1.pack()
entries=[]
def okButton_cmd():
    entries.append(entry.get())
    Label1.config(text=','.join(entries))
entry = AutocompleteEntryListbox(window, width=20, completevalues=df1["Substance"].tolist())
quitButton = tk.Button(window, text="Quit", command=window.quit)
okButton = tk.Button(window, text="OK", command=okButton_cmd)
entry.pack()
quitButton.pack()
okButton.pack()

window.mainloop()
print(entries)