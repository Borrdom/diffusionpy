import xlrd
import matplotlib.pyplot as plt
import numpy as np
filename=r'F:\Dominik Boettcher\DB_KT_25C_ReservoirA_absorptionIDM1original.xls'
book = xlrd.open_workbook(filename, encoding_override = "utf-8")
sheet = book.sheet_by_index(0)
start=39 #Start der der Excel Reihe
nrow=sheet.nrows
nrow2=nrow-start
t_exp = np.asarray([sheet.cell(i+start, 0).value for i in range(nrow2)]) # [s]
t_exp=t_exp*60 #[s]
m_exp = np.asarray([sheet.cell(i+start, 1).value for i in range(nrow2)]) # [mg]
m_exp=m_exp #[g]
a_exp = np.asarray([sheet.cell(i+start, 13).value for i in range(nrow2)]) # [%]
if any(a_exp==87.81):
    tscale=1
else:
    tscale=1
a_exp[a_exp==9.24]=10.0
a_exp[a_exp==29.39]=30.0
a_exp[a_exp==44.49]=45.0
a_exp[a_exp==59.9]=60.0
a_exp[a_exp==73.44]=75.0
a_exp[a_exp==87.81]=90.0



areal_exp = np.asarray([sheet.cell(i+start, 14).value for i in range(nrow2)]) #[%]
Feuchte, indices=np.sort(np.unique(a_exp,return_index=True)) #Extrahiere Sprünge
Feuchte=Feuchte*0.01
#Detect absorbtion and desorption cy#

m1_exp=np.hstack((m_exp[3400:4182],m_exp[14000:]))
t1_exp=np.hstack((t_exp[3400:4182],t_exp[14000:]))
parab=np.poly1d(np.polyfit(t1_exp,m1_exp,2))
#plt.plot(t1_exp,m1_exp)
#plt.plot(t_exp,parab(t_exp))
plt.plot(t_exp,m_exp)
mcorr=m_exp-(parab(t_exp)-m_exp[0])-(m_exp[0]-parab(t_exp)[0])
plt.plot(t_exp,mcorr)