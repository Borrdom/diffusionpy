# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 18:26:18 2020

@author: Moritz Hallermann
"""
import pyodbc
import pandas as pd
import regex as re
from os.path import join
import numpy as np
import os

# subst_input = ['1-butanol','water','ethanol','acetone','Methyl Propanoate','ethyl acetate']

def get_pcpp(subst_input,T=298.15):
 pd.options.mode.chained_assignment = None
 dict_pcpp_all_raw = dict()
 dict_pcpp_all = dict()
 dict_pciap_all_raw = dict()
 dict_pciap_all_organized = dict()
 dict_pciap_all = dict()
 NoneType = type(None)

 del_info = []
 del_info_0 = []
 del_info_1 = []
 del_info_2 = []
 del_info_3 = []
 del_info_4 = []
 del_info_5 = []
 del_info_6 = []
 check_for_dup_keys = []
 dup_keys = []
 dup_pubs = []

 del_info_key_names_1 = ['Substance','Segment Number','Segment Diameter','Energy Parameter',
                        'Associating Volume', 'Associating Energy','H_melt','T_melt','Substance Category',
                        'Publication 1st']
 del_info_key_names_2 = ['kij','Substance Category I','Substance Category II','Publication']
 func_kij_temp = "lambda T: "

 path=join(os.getcwd(),"Nil","parameter.accdb")
 # Input information of the desired Mircosoft Access File (MAF).
 conn_str = (
    r'DRIVER={Microsoft Access Driver (*.mdb, *.accdb)};'
    r'DBQ='+path
    )

 # A connection to the MAF is being generated.
 conn = pyodbc.connect(conn_str)

 # The MAF is opened.
 cur = conn.cursor()

 # A specifc table of the MAF is selected
 pp = 'Reinstoffparameter'
 iap = 'Wechselwirkungsparameter'
 query_pp = "SELECT * from {}".format(pp)
 query_iap = "SELECT * from {}".format(iap)

 # Pandas is accessing the MAF and converts it to a dataframe.
 df_pp = pd.read_sql(query_pp, conn)
 df_iap = pd.read_sql(query_iap, conn)

# First, it is checked whether the input has at least one entry (substance). If so, a for loop iterates through
# the input and searches the input substance in the previously generated dataframe (in the column called
# "Sustance"). While doing so, both the input substance and the substance in the dataframe (originally from
# the MAF) is converted to lower case. This prevents errors and makes the search case-insensitive.
 if len(subst_input) >= 1:
  try:
   for i in subst_input:
    #df_pp_subst = df_pp[df_pp.Substance.str.lower().isin(['{}'.format(i.lower())])]
    df_pp_subst = df_pp[df_pp.Substance.str.lower().str.contains(i.lower())] if not any(df_pp.Substance.str.lower().isin(['{}'.format(i.lower())])) else df_pp[df_pp.Substance.str.lower().isin(['{}'.format(i.lower())])]

    idx_found=df_pp_subst["Substance"].str.lower().str.find(i)
    # bro_found=df_pp_subst["Substance"].str.lower()[25][:idx_found[25]].find("(")
    # brc_found=df_pp_subst["Substance"].str.lower().str.find(")")
    # Abbreviations=

    variants=np.unique(df_pp_subst["Substance"].tolist())
    if len(variants)>1:
        variants
        maxstr=min(df_pp_subst["Substance"].tolist(),key=len)
        df_pp_subst = df_pp_subst[df_pp_subst["Substance"]==maxstr]
    df_pcpp = df_pp_subst[['Substance','Segment Number','Segment Diameter','Energy Parameter',
    'Associating Volume','Associating Energy','H_melt','T_melt','Substance Category','Publication 1st',"Molar Mass","Delta cp","Associating Scheme","Segment Number/Molar Mass","rho0Poly0","ksw","nsw","Tg","rho0"]]
    shouldbenumeric=['Segment Number','Segment Diameter','Energy Parameter','Associating Volume','Associating Energy','H_melt','T_melt',"Molar Mass","rho0Poly0","ksw","nsw","Tg","rho0"]

    for k,valk in enumerate(shouldbenumeric):
     for kk,valkk in enumerate(df_pcpp[valk]):
        if df_pcpp[valk][df_pcpp[valk].index[kk]]!=None and df_pcpp[valk][df_pcpp[valk].index[kk]]!="none":
            if ("kJ/mol" in  df_pcpp[valk][df_pcpp[valk].index[kk]] or "[kJ/mol]" in  df_pcpp[valk][df_pcpp[valk].index[kk]]):
                df_pcpp[valk]=str(float(re.sub("[^0123456789\.]", '', df_pcpp[valk][df_pcpp[valk].index[kk]]))*1000)
            else:
                df_pcpp[valk]=re.sub("[^0123456789\.Ee+-]", '', df_pcpp[valk][df_pcpp[valk].index[kk]])
    #[re.sub('\D', '', 'aas30dsa20')
    dict_pcpp = df_pcpp.to_dict('list')
    dict_pcpp_all_raw[i] = dict_pcpp
  except:
   True
 else:
  print("Missing input!")
  return

# A for loop iterates through each "Substance Category" and if the content is not 'associating' or 'non-polar',
# the key information and the position of the entry is safed as a two entries long list
# (['substance combination', 'position of the entry to be deleted']) in the list called "del_info".
 for k0,v0 in dict_pcpp_all_raw.items():
  for i, value in enumerate(v0['Substance Category']):
   if not (value == 'associating' or value == 'non-polar' or value == 'induced-associating'  #troublesome if none is found
   or value == 'non-associating') and len(del_info)!=len(v0['Substance Category'])-1:
    del_info.append([k0,i])

# A for loop iterates through the list "del_info" to delete all entries which do not belong to the "Substance
# Categroy" 'non-polar' or 'associating' or 'induced-associating' or 'non-associating'. In order to do so, it
# uses the first position "value[0]" as a key, the second postion "value[1]" as an index for the entry which
# should be deleted and this index is substracted by the number of loop in which the for loop currently is. I.e.
# in the first round "i" is zero and in the second it is one. This prevents errors since the indices of the
# entries are changed when one entry is deleted (reduced by 1). The "elif"-condition serves as an aid to check
# if the current substance still belongs to the  substance of the previous loop. If so, the deleting index
# ('value[1]') is reduced by one.
 n = 0
 for i, value in enumerate(del_info):
   if i == 0:
    for j in del_info_key_names_1:
     try:
      del dict_pcpp_all_raw[value[0]][j][value[1]]
     except:
      True

   elif value[0] == del_info[i-1][0]:
    n += 1
    for j in del_info_key_names_1:
     try:
      del dict_pcpp_all_raw[value[0]][j][value[1]-n]
     except:
      True

   else:
    n = 0
    for j in del_info_key_names_1:
     try:
      del dict_pcpp_all_raw[value[0]][j][value[1]-n]
     except:
      True

# A for loop iterates through the dictionary "dict_pciap_all_raw" and extracts the year hidden in the
# "Publication"-String with the help of a "regular expression" ('re'). The extracted year is safed inside
# a list called "year_extracted" and the original "Publication"-String is replaced by the year from this list.
 for k0,v0 in dict_pcpp_all_raw.items():
  for k1,v1 in v0.items():
   if k1 == 'Publication 1st':
    for i, value in enumerate(v1):
     year_extract_info = re.compile(r'\D(?<!-)(\d{4,4})(?!-)\D')
     year_extracted = year_extract_info.findall(" "+str(value)+" ")
     v1[i] = int(year_extracted[0]) if year_extracted!=[] else 1987

# A for loop iterates through all "Publication 1st"-years of each substance. The greatest value of all
# years (in case more than one year exists) is saved as "max_value". Both the substance combination
# (first level key: 'k0') and the indices of all years ('i') except for the year with the greatest value
# (newest publication) is saved in a list which is again saved in a list called "del_info_1".
 for k0,v0 in dict_pcpp_all_raw.items():
  for k1,v1 in v0.items():
   if k1 == 'Publication 1st' and len(v1) > 1:
    max_value = max(v1)
    for i, value in enumerate(v1):
     if value != max_value:
      del_info_1.append([k0,i])

# A For loop iterates through the list "del_info_1" and deletes all information not belonging to the newest
# publication. In order to do so, it uses the first position "value[0]" as a key, the second postion "value[1]"
# as an index for the entry which should be deleted and this index is substracted by the number of loop in which
# the for loop currently is. I. e. in the first round "i" is zero and in the second it is one. This prevents
# errors since the indices of the entries are changed when one entry is deleted (each time reduced by one). The
# "elif"-condition serves as an aid to check if the current substance still belongs to the previous substance.
# If so, the deleting index ('value[1]') is reduced by one.
 n = 0
 for i, value in enumerate(del_info_1):
   if i == 0:
    for j in del_info_key_names_1:
     try:
      del dict_pcpp_all_raw[value[0]][j][value[1]]
     except:
      True
   elif value[0] == del_info_1[i-1][0]:
    n += 1
    for j in del_info_key_names_1:
     try:
      del dict_pcpp_all_raw[value[0]][j][value[1]-n]
     except:
      True
   else:
    n = 0
    for j in del_info_key_names_1:
     try:
      del dict_pcpp_all_raw[value[0]][j][value[1]-n]
     except:
      True

# To rearrange "dict_pcpp_all_raw" in a way that the values are not inside a list but "free" values a new
# dictionary called "dict_pcpp_all" is generated. Inside the new dictionary each key is connceted to the value
# inside the list of the original dictionary "dict_pcpp_all_raw". To prevent overriding, the function
# "setdefault()" is used.
 for k0,v0 in dict_pcpp_all_raw.items():
  for k1,v1 in v0.items():
   for i in v1:
    dict_pcpp_all.setdefault(k0,{}).update({k1:i})

# All values whose type is a string - actually is a float though - is converted to a float.
 for k0,v0 in dict_pcpp_all.items():
  for k1,v1 in v0.items():
   if not isinstance(v1, int):
    try:
     dict_pcpp_all[k0][k1] = float(v1)
    except:
     ValueError


#---------------------------------------------------------------


# The MAF needs to be sorted in a manner that it matches the structure of "subst_comb" (i.e. "a + b", is first
# introduced in the following chapter). The first section of code lines considers positively charged iones
# which are replaced before the earlier mentioned structure is inforced. The second section of code lines
# considers the case that substances are already matching the desired structure and replaces them temporarily
# as well. The now remaining and wrongly structured substances are corrected. Afterwards the replaced positively
# charged ions and already correct substances are reintroduced.
 df_iap['Binary Mixture'] = df_iap['Binary Mixture'].apply(lambda x: x.replace('Na+','Na_Plus'))
 df_iap['Binary Mixture'] = df_iap['Binary Mixture'].apply(lambda x: x.replace('Li+','Li_Plus'))
 df_iap['Binary Mixture'] = df_iap['Binary Mixture'].apply(lambda x: x.replace('Ch2+','CH2_Plus'))
 df_iap['Binary Mixture'] = df_iap['Binary Mixture'].apply(lambda x: x.replace('H3O+','H3O_Plus'))
 df_iap['Binary Mixture'] = df_iap['Binary Mixture'].apply(lambda x: x.replace('NH4+','NH4_Plus'))
 df_iap['Binary Mixture'] = df_iap['Binary Mixture'].apply(lambda x: x.replace('NAD+','NAD_Plus'))
 df_iap['Binary Mixture'] = df_iap['Binary Mixture'].apply(lambda x: x.replace('Mg+','Mg_Plus'))
 df_iap['Binary Mixture'] = df_iap['Binary Mixture'].apply(lambda x: x.replace('Ca2+','Ca2+_Plus'))
 df_iap['Binary Mixture'] = df_iap['Binary Mixture'].apply(lambda x: x.replace('K+','K_Plus'))
 df_iap['Binary Mixture'] = df_iap['Binary Mixture'].apply(lambda x: x.replace('[C_n_mim]+','[C_n_mim]_Plus'))
 df_iap['Binary Mixture'] = df_iap['Binary Mixture'].apply(lambda x: x.replace('[C_n_py]+','[C_n_py]_Plus'))
 df_iap['Binary Mixture'] = df_iap['Binary Mixture'].apply(lambda x: x.replace('[C_n_mpy]+','[C_n_mpy]_Plus'))
 df_iap['Binary Mixture'] = df_iap['Binary Mixture'].apply(lambda x: x.replace('[C_n_mpyr]+','[C_n_mpyr]_Plus'))

 df_iap['Binary Mixture'] = df_iap['Binary Mixture'].apply(lambda x: x.replace(' + ','PLUS'))
 df_iap['Binary Mixture'] = df_iap['Binary Mixture'].apply(lambda x: x.replace('+',' + '))
 df_iap['Binary Mixture'] = df_iap['Binary Mixture'].apply(lambda x: x.replace('PLUS',' + '))

 df_iap['Binary Mixture'] = df_iap['Binary Mixture'].apply(lambda x: x.replace('Na_Plus','Na+',))
 df_iap['Binary Mixture'] = df_iap['Binary Mixture'].apply(lambda x: x.replace('Li_Plus','Li+',))
 df_iap['Binary Mixture'] = df_iap['Binary Mixture'].apply(lambda x: x.replace('CH2_Plus','Ch2+'))
 df_iap['Binary Mixture'] = df_iap['Binary Mixture'].apply(lambda x: x.replace('H3O_Plus','H3O+'))
 df_iap['Binary Mixture'] = df_iap['Binary Mixture'].apply(lambda x: x.replace('NH4_Plus','NH4+'))
 df_iap['Binary Mixture'] = df_iap['Binary Mixture'].apply(lambda x: x.replace('NAD_Plus','NAD+'))
 df_iap['Binary Mixture'] = df_iap['Binary Mixture'].apply(lambda x: x.replace('Mg_Plus','Mg+'))
 df_iap['Binary Mixture'] = df_iap['Binary Mixture'].apply(lambda x: x.replace('Ca2+_Plus','Ca2+'))
 df_iap['Binary Mixture'] = df_iap['Binary Mixture'].apply(lambda x: x.replace('K_Plus','K+'))
 df_iap['Binary Mixture'] = df_iap['Binary Mixture'].apply(lambda x: x.replace('[C_n_mim]_Plus','[C_n_mim]+' ))
 df_iap['Binary Mixture'] = df_iap['Binary Mixture'].apply(lambda x: x.replace('[C_n_py]_Plus','[C_n_py]+'))
 df_iap['Binary Mixture'] = df_iap['Binary Mixture'].apply(lambda x: x.replace('[C_n_mpy]_Plus','[C_n_mpy]+'))
 df_iap['Binary Mixture'] = df_iap['Binary Mixture'].apply(lambda x: x.replace('[C_n_mpyr]_Plus','[C_n_mpyr]+'))

# In a for loop the entries of "subst_input" are put together as substance combinations (combinations with itself
# is avoided by the "if"-condition, first with second, first with third and so on). In this manner all substance
# combination are searched. A new dataframe called "df_iap_subst" is generated. Then, the columns called "kij",
# "Substance Catecory I", "Substance Category II" and "Publication" of "df_iap_subst" are extracted and saved in
# a dataframe called "df_pciap". The dataframe "df_pciap" is converted into a dictionary and is eventually added
# to the final dictionary "dict_pciap_all_raw" which collects all previously extracted columns for each
# substance combination of "subst_input".
 df_iap['Binary Mixture']=df_iap['Binary Mixture'].str.replace(" ","") #troublesome if space in subtsance
 subst_input2=[dict_pcpp_all[val]["Substance"] for i,val in enumerate(dict_pcpp_all)]
 if len(subst_input2) > 1:
  for idxi,i in enumerate(subst_input2):
   for idxn,n in enumerate(subst_input2):
    if not i == n:
     booli=df_iap['Binary Mixture'].str.lower().str.contains(subst_input[idxi])
     boolj=df_iap['Binary Mixture'][booli].str.lower().str.contains(subst_input[idxn])
     df_iap_subst=df_iap[booli][boolj]

     boolki=np.asarray([],dtype=bool)
     boolkj=np.asarray([],dtype=bool)
     for k,val in enumerate (df_iap_subst["Binary Mixture"].str.lower()):
         streidxi=val[val.find(subst_input[idxi])+len(subst_input[idxi])] if val.find(subst_input[idxi])+len(subst_input[idxi])<len(val) else "1"
         streidxj=val[val.find(subst_input[idxn])+len(subst_input[idxn])] if val.find(subst_input[idxn])+len(subst_input[idxn])<len(val) else "1"
         booi=not streidxi.isalpha()
         booj=not streidxj.isalpha()
         boolki=np.hstack((boolki,booi))
         boolkj=np.hstack((boolkj,booj))
     df_iap_subst=df_iap_subst[boolki&boolkj]

     #subst_comb = (n + "+"  + i).replace(" ","")
     #subst_comb2 = (subst_input[idxn]  + "+"  + subst_input[idxi]).replace(" ","")
     subst_comb = (subst_input[idxn]  + "+"  + subst_input[idxi]).replace(" ","")
     #df_iap_subst = df_iap[df_iap['Binary Mixture'].str.lower().isin(['{}'.format(subst_comb.lower())])] if any(df_iap['Binary Mixture'].str.lower().isin(['{}'.format(subst_comb.lower())])) else df_iap[df_iap['Binary Mixture'].str.lower().isin(['{}'.format(subst_comb2.lower())])]
     df_pciap = df_iap_subst[['kij',"ksw",'Substance Category I','Substance Category II','Publication']]
     dict_pciap = df_pciap.to_dict('list')
     dict_pciap_all_raw[subst_comb] = dict_pciap
     #Still makes problems for pvpk25 instead if pvp need to think about sth else
 #print('\n')
 #print(dict_pciap_all_raw)

# The for loop adds the key value information of those entries which are empty (because they do not appear in
# the original MAF) to a list called "del_info_2". In the following for loop this information is used to delete
# the empty entries.
 for k0,v0 in dict_pciap_all_raw.items():
  for k1,v1 in v0.items():
   if not v1:
    del_info_2.append(k0)

 for i in del_info_2:
  try:
   del dict_pciap_all_raw[i]
  except KeyError:
   True

# A for loop iterates through the dictionary "dict_pciap_all_raw" and if either "Substance Category I"
# (substance category of the first appearing substance in the substance combination name, thus '[0]') or
# "Substance Category II" (same logic but the second substance, thus '[2]') occurs, it is checked whether the
# "Substance Category" of the pure substance from "dict_pcpp_all" is contained in the corresponding
# "Substance Category I or II (depending on the case)" from "dict_pciap_all_raw". Since the "Substance Category"
# is always 'associating' or 'non-polar' (or the other cases) (due to the above introduced restriction) but
# the "Substance Category I or II" can be for example 'associating (2B)', the containment is checked rather than
# the exact match. This is done with the aid of a 'regular expression'. To receive either the first substance
# name or the second substance name of the substance combination ('Substance_A' + 'Substance_B') the string is
# split ['.split()'] into three substrings (0: 'Substance_A, 1:'+', 2:'Substance_B'). Depending on the case
# either the first substring (Substance_A -> '[0]') or the last substring (Substance_B -> '[2]') is used. If the
# "Substance Category" is not contained in "Substance Category I or II" the key information and the position of
# the entry is safed as a two entries long list (['Substance combination', Position of the entry to be deleted]) in the list
# called "del_info_3".
 for k0,v0 in dict_pciap_all_raw.items():
  if v0 == 'Substance Category I':
   for i, value in enumerate(v0['Substance Category I']):
    if not re.search(dict_pcpp_all[k0.split()[0]]['Substance Category'], value):
     del_info_3.append([k0,i])
  elif v0 == 'Substance Category II':
   for i, value in enumerate(v0['Substance Category II']):
    if not re.search(dict_pcpp_all[k0.split()[2]]['Substance Category'], value):
     del_info_3.append([k0,i])

# A for loop iterates through the list "del_info_3" to delete all entries which do not belong to a "Substance
# Categroy I or II" 'non-polar' or 'associating' (or the other cases). In order to do so, it uses the first
# position  "value[0]" as a key, the second postion "value[1]" as an index for the entry which should be deleted
# and this index is substracted by the number of loop in which the for loop currently is. I.e. in the first round
# "i" is zero and in the second it is one. This prevents errors since the indices of the entries are changed when
# one  entry is deleted (reduced by one).
 n = 0
 for i, value in enumerate(del_info_3):
  if i == 0:
   for j in del_info_key_names_2:
    try:
     del dict_pciap_all_raw[value[0]][j][value[1]]
    except:
     True
  elif value[0] == del_info_3[i-1][0]:
   n += 1
   for j in del_info_key_names_2:
    try:
     del dict_pciap_all_raw[value[0]][j][value[1]-n]
    except:
     True
  else:
   n = 0
   for j in del_info_key_names_2:
    try:
     del dict_pciap_all_raw[value[0]][j][value[1]-n]
    except:
     True

# A for loop is used to evaluate if there is more than one publication and if all the publications are the same.
# If so, all entries except for the first (index "0") are deleted.
 for k0,v0 in dict_pciap_all_raw.items():
  if len(v0['Publication']) > 1 and len(set(v0['Publication'])) == 1:
   del v0['Publication'][1:]
   del v0['kij'][1:]
   del v0['Substance Category I'][1:]
   del v0['Substance Category II'][1:]

# A for loop iterates through the dictionary "dict_pciap_all_raw" and extracts the year hidden in the
# "Publication"-String with the help of a "regular expression" ('re'). The extracted year is safed inside a
# list called "year_extracted" and the original "Publication"-String is replaced by the year from this list.
 for k0,v0 in dict_pciap_all_raw.items():
  for k1,v1 in v0.items():
   if k1 == 'Publication':
    for i, value in enumerate(v1):
     year_extract_info = re.compile(r'\D(?<!-)(\d{4,4})(?!-)\D')
     year_extracted = year_extract_info.findall(" "+value+" ")
     try:
      v1[i] = int(year_extracted[0])
     except:
      v1[i] = 0

# A for loop iterates through the keys of the first level of "dict_pciap_all_raw" (names of substance
# combinations) and splits the substrings seperated by a space into substrings. Those substrings are then sorted
# to receive a chronological order of the substrings. This is done to be able to indentify substance combination
# names which are just turned vice versa (i.e. "Substance_A + Substance_B" has the same interactions as
# "Substance_B + Substance_A"). In those cases the received "dup_key" (='duplicate key' which is split and then
# sorted) is the same as for 'hidden duplicates' because "Substance_A + Substance_B" as well as "Substance_B +
# Substance_A" are turned to "['+','Substance_A','Substance_B']. The original key (first level) and the
# "dup_key" are both added inside a list to another the list called "check_for_dup_keys". Before adding this
# combination it is checked if "dup_key" is already present in the list "check for_dup_keys". If not, it is
# added and if so, the index of the already excisting "duplicate key" is extracted by using the function
# "next()" which is necessary in order to get the index of an element in a sublist. Eventually, the original key
# (i.e. "Substance_A + Substance_B")  and the orginal "hidden duplicate" (i.e. "Substance_B + Substance_A") are
# added to the list called "dup_keys". Now the information which keys are "hidden duplicates" and at which
# position they are is known/saved.
 for k0,v0 in dict_pciap_all_raw.items():
  dup_key = sorted(k0.split())
  if not dup_key in (item for sublist in check_for_dup_keys for item in sublist):
   check_for_dup_keys.append([k0,dup_key])
  else:
   try:
    dup_key_index = next(index for index, value in enumerate(check_for_dup_keys) if dup_key in value)
    dup_keys.append(k0)
    dup_keys.append(check_for_dup_keys[dup_key_index][0])
   except:
    StopIteration

# The for loop uses the previously generated list "dup_keys" which contains the keys of the first level of the
# dictionary "dict_pciap_all_raw" (names of the substance combinations). However, those keys are connected
# directly in series and to avoid mix-ups a counter "k" is introduced. The counter is reset to 1 after two loops
# because each pair of duplicates is represented exactly by two entries
# (i.e. ['water + ethanol','ethanol + water']). So after two loops a new substance combination begins. For each
# substance combination pair (duplicates) the key of the first level (i = first entry), the publication year
# (value = second entry) and the position of the publication year inside the dictionary-list (n = third entry)
# is added to a list called "dup_pubs". In case two loops are finished (k = 3), the maximum value of the
# second entries is extracted by the function "max()". With a generator expression the index of  this maximum
# value is determined. With the aid of this index, the entry of the maximum value is deleted from "dup_pubs".
# Eventually, another for loop adds the remaining entries of "dup_pubs" to the list "del_info_4" and then clears
# the list "dup_pubs".
 k = 1
 for i in dup_keys:
  if k < 3:
   k += 1
   for n, value in enumerate(dict_pciap_all_raw[i]['Publication']):
    dup_pubs.append([i,value,n])
   if k == 3:
    max_pub = max(i[1] for i in dup_pubs)
    max_pub_index = next(index for index, value in enumerate(dup_pubs) if max_pub in value)
    dup_pubs.pop(max_pub_index)
    for j in dup_pubs:
     del_info_4.append(j)
    dup_pubs.clear()
    k = 1
  else:
   print("An unexpected error occured!")

# The previously generated list "del_info_4" containing the information of all substance combinations (including
# duplicates) is used to delete everything which does not belong to the newest publication.
 n = 0
 for i, value in enumerate(del_info_4):
  if i == 0:
   for j in del_info_key_names_2:
    try:
     del dict_pciap_all_raw[value[0]][j][value[1]]
    except:
     True
  elif value[0] == del_info_4[i-1][0]:
   n += 1
   for j in del_info_key_names_2:
    try:
     del dict_pciap_all_raw[value[0]][j][value[1]-n]
    except:
     True
  else:
   n = 0
   for j in del_info_key_names_2:
    try:
     del dict_pciap_all_raw[value[0]][j][value[1]-n]
    except:
     True

# The for loop adds the key value information of those entries which might be empty (after checking for
# duplicates and deleting all older information than the those from the newest publication) to a list called
# del_info_5. In the following for loop this information is used to delete the empty entries.
 for k0,v0 in dict_pciap_all_raw.items():
  for k1,v1 in v0.items():
   if not v1:
    del_info_5.append(k0)

 for i in del_info_5:
  try:
   del dict_pciap_all_raw[i]
  except KeyError:
   True

# A for loop iterates through all "Publication"-years of each substance combination. The greatest value of all
# years (in case more than one year exists) is saved as "max_value". Both the substance combination ('k0') and
# the indices of all years ('i') except for the year with the greatest value (newest publication) is saved in a
# list which is saved in a list called "del_info_5".
 for k0,v0 in dict_pciap_all_raw.items():
  for k1,v1 in v0.items():
   if k1 == 'Publication' and len(v1) > 1:
    max_value = max(v1)
    for i, value in enumerate(v1):
     if value != max_value:
      del_info_6.append([k0,i])

# All information which belongs to years older than the year of the newest publication is deleted. While
# deleting, the original indices are changing. Errors are prevented by reducing the index from "del_info_4"
# by the number of the current iteration loop. If a new substance combination is reached the counter 'n' is
#  reset to zero.
 n = 0
 for i, value in enumerate(del_info_6):
  if i == 0:
   for j in del_info_key_names_2:
    try:
     del dict_pciap_all_raw[value[0]][j][value[1]]
    except:
     True
  elif value[0] == del_info_6[i-1][0]:
   n += 1
   for j in del_info_key_names_2:
    try:
     del dict_pciap_all_raw[value[0]][j][value[1]-n]
    except:
     True
  else:
   n = 0
   for j in del_info_key_names_2:
    try:
     del dict_pciap_all_raw[value[0]][j][value[1]-n]
    except:
     True

# The for loop fills a new dictionary called "dict_pciap_all_organized" with the names of the substance
# combination (from "dict_pciap_all_raw") and adds a list with sublists for three kinds of kij. The first is a
# constant value (kij_c), the second value is linear-dependent on the temperature (kij_T) and the last value is
# quadratic-dependent on the temperatur (kij_T2). This generates a dictionary where a list of lists (three kinds
# of kijs) belongs to each substance combination.
 for k0,v0 in dict_pciap_all_raw.items():
  dict_pciap_all_organized.update({k0:[['kij_c'],['kij_T'],['kij_T2'],["ksw"]]})

# The for loop adds to the first list (kij_c) of the dictionary "dict_pciap_all_organized" only those entries
# which do not contain the letter "T". To the third list (kij_T2) only those entries are added which contain
# "T²". All other cases need to be linear-dependent on the temperature and thus are added to the the second list
# (kij_T).
 for k0,v0 in dict_pciap_all_raw.items():
  for k1,v1 in v0.items():
    if k1 == 'kij':
     for i in v1:
      if not 'T' in i:
       dict_pciap_all_organized[[k0][0]][0].append(i.replace('²','**2').replace('×','*').replace('−','-')\
       .replace('[K]','').replace('[K²]','').replace('(K)','').replace('K','').replace('^','**').replace(',','.'))
      elif 'T²' in i:
       dict_pciap_all_organized[[k0][0]][2].append(i.replace('²','**2').replace('×','*').replace('−','-')\
       .replace('[K]','').replace('[K²]','').replace('(K)','').replace('K','').replace('^','**').replace(',','.'))
      else:
       dict_pciap_all_organized[[k0][0]][1].append(i.replace('²','**2').replace('×','*').replace('−','-')\
       .replace('[K]','').replace('[K²]','').replace('(K)','').replace('K','').replace('^','**').replace(',','.'))
    elif k1 == 'ksw':
      for i in v1:
       i="0" if i==None else i
       dict_pciap_all_organized[[k0][0]][3].append(i.replace('²','**2').replace('×','*').replace('−','-')\
       .replace('[K]','').replace('[K²]','').replace('(K)','').replace('K','').replace('^','**').replace(',','.'))

# The for loop iterates through the three kinds of kijs and if the lists consists of more than two
# entries (which means more than one value for the same kij) it deletes all entries except of the first
# two entries. If one of the list only consists of one entry (which means there is no value for kij) it
# appends a zero.
 for k0,v0 in dict_pciap_all_organized.items():
  for v1 in v0:
    if len(v1) > 2:
     del v1[2:]
    elif len(v1) == 1:
     v1.append(0)

# The dictionary "dict_pciap_all_organized" has stored the name for kij-types ('kij_c','kij_T','kij_T2') in the
# same list as the values for the kijs. The for loop is used to rearrange the dictionary in way that the
# substance combination is the key on the first level ('k0') and the type of kij (1:kij_c, 2:kij_T, 3:kij_T2)
# the key on the second level ('k1'). The respective value of each kij should be linked to it as a value. To
# achieve this, a new dictionary "dict_pciap_all" is generated. To the key of the first level ('k0') a
# dictionary is linked which is filled with the key of the second level ('k1') connected to its value ('v1').
# To prevent overriding during the for loop "setdefault()" is used. Otherwise, only the last key and value of
# each substance combination of the original dictionary "dict_pciap_all" is added to the new dictionary
# "dict_pciap_all". The other information would be overridded.
 for k0,v0 in dict_pciap_all_organized.items():
  for k1,v1 in v0:
    dict_pciap_all.setdefault(k0,{}).update({k1:v1})

# A for loop evaluates if the value for kij_c is zero and if - at the same time - the value for kij_T is given.
# If so, the value for 298,15 K of the function which is given as a string for kij_T is calculated and set as
# the value for the kij_c. If kij_c and kij_T are zero but kij_T2 has a value, the value for 298,15 K of the
# respective function is calculated and also set as kij_c.
 for k0,v0 in dict_pciap_all.items():
  if v0['kij_c'] == 0 and v0['kij_T'] != 0:
    f = func_kij_temp + v0['kij_T']
    v0['kij_c'] = eval(f)(T)
  elif v0['kij_c'] == 0 and v0['kij_T'] == 0 and v0['kij_T2'] != 0:
    f = func_kij_temp + v0['kij_T2']
    v0['kij_c'] = eval(f)(T)

# All values of kij_c which are not of the type "float" are converted to a float.
# for k0,v0 in dict_pciap_all.items():
  if not isinstance(v0['kij_c'], float):
    v0['kij_c'] = float(v0['kij_c'])
  if not isinstance(v0['ksw'], float):
    v0['ksw'] = float(v0['ksw'])

#Evaluate association scheme
 for i,val in enumerate(dict_pcpp_all):
  tru=[valj for j,valj in enumerate(str(dict_pcpp_all[val]["Associating Scheme"])) if valj.isdigit()]
  tru2=[valj.isalpha() for j,valj in enumerate(str(dict_pcpp_all[val]["Associating Scheme"]))]
  if any(tru2):
   if dict_pcpp_all[val]["Associating Scheme"]=="none" or dict_pcpp_all[val]["Associating Scheme"]==None:
      dict_pcpp_all[val]["Associating Scheme"]=1
   else:
      dict_pcpp_all[val]["Associating Scheme"]=int(dict_pcpp_all[val]["Associating Scheme"][0])//2
  else:
   dict_pcpp_all[val]["Associating Scheme"]=int("".join(tru[0:len(tru)//2]))
  if dict_pcpp_all[val]["Segment Number"]=="none" or dict_pcpp_all[val]["Segment Number"]==None:
      dict_pcpp_all[val]["Segment Number"]=dict_pcpp_all[val]["Segment Number/Molar Mass"]*dict_pcpp_all[val]["Molar Mass"]
  dict_pcpp_all[val]["Molar Mass"]=dict_pcpp_all[val]["Molar Mass"]/1000
  #if dict_pcpp_all[val]["Substance Category"]=="associating" and (dict_pcpp_all[val]["Associating Energy"]=="none" or dict_pcpp_all[val]["Associating Energy"]==None):
  if (dict_pcpp_all[val]["Associating Energy"]=="none" or dict_pcpp_all[val]["Associating Energy"]==None):
   dict_pcpp_all[val]["Associating Energy"]=0
  #if dict_pcpp_all[val]["Substance Category"]=="associating" and (dict_pcpp_all[val]["Associating Volume"]=="none" or dict_pcpp_all[val]["Associating Volume"]==None):
  if (dict_pcpp_all[val]["Associating Volume"]=="none" or dict_pcpp_all[val]["Associating Volume"]==None):
   dict_pcpp_all[val]["Associating Volume"]=0.02
  # if (dict_pcpp_all[val]["H_melt"]!="none" and dict_pcpp_all[val]["H_melt"]!=None):
  #  if ("kJ/mol" in dict_pcpp_all[val]["H_melt"] or "[kJ/mol]" in dict_pcpp_all[val]["H_melt"]):
  #   dict_pcpp_all[val]["H_melt"]=float(re.sub('\D',"",dict_pcpp_all[val]["H_melt"]))
  #dict_pcpp_all[val]["Tg"]=383.9
  #dict_pcpp_all[val]["rho0"]=1180
  #dict_pcpp_all[val]["rho0Poly0"]=1355.5190482692865
  #dict_pcpp_all[val]["ksw"]=3.4376005848126273e-08
  #dict_pcpp_all[val]["name"]=subst_input[i]
  dict_pcpp_all[val]["Substance"]=subst_input[i]
  if dict_pcpp_all[val]["Substance"]=="water":
   sigmawater=2.7927+10.11*np.exp(-0.01775*T)-1.417*np.exp(-0.01146*T)
   dict_pcpp_all[val]["Segment Diameter"]=sigmawater


 return dict_pcpp_all, dict_pciap_all

def get_par(subst_input,T=298.15):
    unique_input=pd.unique(subst_input).tolist()
    subst_output,subst_kij=get_pcpp(unique_input,T=T)

    oldtonew=[("Substance","name"),("Molar Mass","Mi"),("Segment Diameter","sigi"),("Segment Number","mi"),("Energy Parameter","ui"),("Associating Energy","epsAiBi"),("Associating Volume","kapi"),("Associating Scheme","N"),("T_melt","Tm"),("H_melt","deltaHSL"),("Delta cp","cpSL")]
    for key,value in enumerate(subst_output):
     for i,val in enumerate (oldtonew):
      subst_output[value][val[1]]=subst_output[value].pop(val[0])
    purelist=tuple([subst_output[value] for key,value in enumerate(subst_output)])
    notunique=[]
    for vali in subst_input:
        for valj in purelist:
            notunique.append(valj) if vali==valj["name"]  else None
    purelist=tuple(notunique)
    kijnew=dict([(value.replace("+",""),subst_kij[value]["kij_c"]) for key,value in enumerate(subst_kij)])
    kswnew=dict([(value.replace("+",""),subst_kij[value]["ksw"]) for key,value in enumerate(subst_kij)])


    return purelist,kijnew,kswnew


if __name__=="__main__":
    subst_input_raw = input("Enter the substance or the substances. Enter one after another seperated by a semikolon (i.e. ethanol;acetone;1-butanol and so on): ")
    subst_input = subst_input_raw.split(";")
    pure,kij=get_par(subst_input)

    print('\n')
    print(pure)
    print(kij)
