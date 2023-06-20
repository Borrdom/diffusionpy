
import pyodbc
import pandas as pd
from tkinter.filedialog import askopenfilename
from textdistance import levenshtein
from itertools import combinations
import numpy as np
def levensthein_pd(df,key1,key2):
   """calculate the lvensthein distance for a pandas dataframe"""   
   df["Levensthein distance"]=df.apply(lambda x: levenshtein.distance(x[key1].lower(),key2), axis=1)#.replace('\d+', '')
   choice=df["Levensthein distance"].idxmin()
   choice2=np.where(df["Substance"].str.replace('\d+', '').str.contains(r'\b{}\b'.format(key2),case=False,regex=True))[0]
   choice2= choice if choice2.size==0 else choice2[0]
   choice= choice2 if df["Levensthein distance"].min()>1 else choice
   return df.iloc[[choice],:]


def get_par(subst_input,path,header_1=None,header_2=None):
   """get parameters from a database"""  
   conn_str = (r'DRIVER={Microsoft Access Driver (*.mdb, *.accdb)};'r'DBQ='+path)
   conn = pyodbc.connect(conn_str)
   sheets=[row.table_name for row in conn.cursor().tables()]
   unary=sheets[-2]
   binary=sheets[-1]
   query_unary = "SELECT * from {}".format(unary)
   query_binary = "SELECT * from {}".format(binary)
   df1 = pd.read_sql(query_unary, conn)
   df2 = pd.read_sql(query_binary, conn)
   unaryname=df1.columns[0]
   binaryname=df2.columns[0]
   header_1=df1.columns if header_1 is None else header_1
   header_2=df2.columns if header_2 is None else header_2
   appended_data = []
   appended_data2 = []
   for i,val in enumerate(subst_input):
      #appended_data.append(df1[df1[unaryname].str.contains(r'\b{}\b'.format(i),case=False,regex=True)])
      df1_search=levensthein_pd(df1,unaryname,val)
      appended_data.append(df1_search)
      #subst_input[i]=df1_search[unaryname].iloc[0]
   combi=list(combinations(subst_input,2))
   for i,j in combi:
      #appended_data2.append(levensthein_pd(df2,binaryname,i+"+"+j))
      df2_search=df2[df2[binaryname].str.replace('\d+', '').str.contains(r'(?=.*\b{}\b)(?=.*\b{}\b)'.format(i,j),case=False,regex=True)]
      df2_search=df2_search.append({"Binary Mixture": i+"+"+j, "kij":0},ignore_index=True) if df2_search.empty else df2_search
      df2_search[unaryname+" I"]=i
      df2_search[unaryname+" II"]=j
      appended_data2.append(df2_search.iloc[[0],:])
   appended_data = pd.concat(appended_data)
   appended_data2 = pd.concat(appended_data2)
   return appended_data[header_1], appended_data2[header_2]

if __name__=="__main__":
   subst_input = ["water","pvpva","hexane"]
   header_1 =['Substance','Segment Number','Segment Diameter','Energy Parameter','Associating Volume','Associating Energy','H_melt','T_melt',"Molar Mass","Delta cp","Segment Number/Molar Mass"]
   header_2 = ['Binary Mixture','kij','Substance Category I','Substance Category II','Publication']
   a,b=get_par(subst_input,askopenfilename(),header_1, header_2)