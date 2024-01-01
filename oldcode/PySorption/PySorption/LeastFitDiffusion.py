from matplotlib import pyplot as plt
import numpy as np
#import LeastFitDiffusionFunktionen as LDF
import os
from PySorption.DiffusionFinal import Film
cwd=os.getcwd()
from PySorption.dcc import window
import matplotlib as mpl
from PySorption.DiffusionsmodellPlots import Plotlist
import pandas as pd
import pdb
from PyCSAFT.list_manipulations import unnesting,averagedict,averagelist,averagelisttimeseries
from collections import OrderedDict
from scipy.interpolate import InterpolatedUnivariateSpline
#from sklearn.gaussian_process import GaussianProcessRegressor
#Least Square Regression for data provided by lis
def FitDiffusionData(Film1,Data,Fit=True):
    if Film1.Mode=="MEOS" and Fit==True:
        D12True=input("Vary D12? <<y>> or <<n>>\n")=="y"
        exTrue=input("Vary exponent? <<y>> or <<n>>\n")=="y"
        ETrue=input("Vary E? <<y>> or <<n>>\n")=="y"
        tauTrue=input("Vary tau? <<y>> or <<n>>\n")=="y"

    nJumplist=[]
    # Matrix=Film1.comp[0]
    # water=Film1.comp[1]
    # api=Film1.comp[2]
    keylist=[]
    DlistFick,DlistCranc,wreallist,Errorlist,Feuchtelist,tveclist,wtlist,wtreflist,wrealmittellist,mveclist,sratelist,sraterlist,sraterrlist,ARDlist,ARDrlist,taulist,explist,Elist=([] for i in range(18))
    backend_=mpl.get_backend()
    mpl.use("Agg")
    Data={0:Data,1:Data} if isinstance(Data,str) else Data
    for ppp in range(len(Data)):
        DCranc,tvecsol,wtrefsol,wtsol,Dsol,tausol,expsol,Esol=([] for i in range(8))
        #Film1=Film(Data[ppp],Matrix,water)
        Film1.__init__(Data[ppp],Film1.comp,Film1.wPolyASD)
        #Film1.Mode="NEAT"
        #Film1.Mode="NEAT"
        Film1.ideal=1 #if Fit else 0
        Cranc=False#True if not Fit else False
        #Film1.endpoint=True if Fit else False
        #Film1.exp_boundary()


        # Film1.MY=0
        # Film1.NET=0
        # Film1.RecalcNET=0 #if ppp==0 else 0

        Film1.exp_boundary()
        wreal=Film1.wreal
        Feuchtereal=Film1.Feuchtereal
        nJumplist.append(Film1.nHum)
        mvecsol=[]
        sratesol=[]
        sratersol=[]
        sraterrsol=[]
        ARDsol=[]
        ARDrsol=[]
        for iii in range(nJumplist[ppp]):
            #DWASDFEED=Film1.CalcTerDiffcoeff()
            #Film1.D12=DWASDFEED
            plt.close("all")

            #iii=1 if iii==0 else iii

            Film1.indJ=iii
            Film1.Cranc() if Film1.Mode!="MEOS" and Fit!=True else None
            DCranc=np.append(DCranc,Film1.D12) 

            if Fit:
                Film1.LeastMEOS(D12True,exTrue,ETrue,tauTrue) if Film1.Mode=="MEOS" else Film1.Least()

            else:
                Film1.CalcTerDiffcoeff() if Film1.Mode!="MEOS" else None

            Film1.Diffusion1D() if not Cranc else Film1.Cranc(Fit=Film1.D12)
            Dsol=np.append(Dsol,Film1.D12)
            if Film1.Mode=="MEOS":
                tausol=np.append(tausol,Film1.eta0)
                expsol=np.append(expsol,Film1.exponent)
                Esol=np.append(Esol,Film1.E0)
            tvecsol.append(Film1.tlist[iii]/60)
            mvecsol.append(Film1.mlist[iii]) #if not Film1.des else mvecsol.append(Film1.mlist[iii][0]+Film1.mlist[iii][-1]-Film1.mlist[iii])
            wtrefsol.append(Film1.wtref) #if not Film1.des else wtrefsol.append(Film1.wtref[0]+Film1.wtref[-1]-Film1.wtref)
            wtsol.append(Film1.wt) #if not Film1.des else wtsol.append(Film1.wt[0]+Film1.wt[-1]-Film1.wt)

            ttemp=Film1.tlist[iii]/60
            ARDtemp=np.abs((Film1.wtref-Film1.wt))/(Film1.wtref)
            ARDsol.append(ARDtemp)
            ARDrsol.append(np.average(ARDtemp))
            sratesol.append(Film1.sratel[iii])
            sratersol.append(np.average(Film1.sratel[iii][ttemp>60]))
            sraterrsol.append(np.average(Film1.sratel[iii][0]))
            #sratersol.append(np.average(Film1.sratel[iii][Film1.wtref[-1]*0.99<Film1.wtref]))
            #sraterrsol.append(np.average(Film1.sratel[iii][Film1.wtref[-1]*0.99>Film1.wtref]))

            #wrealmittel=np.append(wrealmittel,Film1.wrealmittelref)Changed on 06.04.20
            #wrealmittel.append(Film1.wrealmittelref)
        DlistFick.append(Dsol)
        wreallist.append(np.asarray(wreal))
        wrealmittellist.append(np.asarray(Film1.wrealmittelref))
        #wrealmittellist.append(wreal)
        Feuchtelist.append(np.asarray(Feuchtereal))
        DlistCranc.append(DCranc)

        tveclist.append(tvecsol)
        mveclist.append(mvecsol)
        sratelist.append(sratesol)
        sraterlist.append(sratersol)
        sraterrlist.append(sraterrsol)
        ARDlist.append(ARDsol)
        ARDrlist.append(ARDrsol)
        wtreflist.append(wtrefsol)
        wtlist.append(wtsol)
        keylist.append(Film1.dickeys)
        if Film1.Mode=="MEOS":
            taulist.append(tausol)
            explist.append(expsol)
            Elist.append(Esol)
    ###### Added 28.04.20

    tvecdict=[]
    wrealdict=[]
    wrealmitteldict=[]
    Feuchtedict=[]
    wtrefdict=[]
    wtdict=[]
    DlistFickdic=[]
    DlistCrancdic=[]
    taudict=[]
    expdict=[]
    Edict=[]
    sratedict=[]
    sraterdict=[]
    sraterrdict=[]
    ARDdict=[]
    ARDrdict=[]
    mpl.use(backend_)
    for ppp in range(len(Data)):
        wrealdict.append(OrderedDict(zip( keylist[ppp],wreallist[ppp])))
        tvecdict.append(OrderedDict(zip( keylist[ppp],tveclist[ppp])))
        wrealmitteldict.append(OrderedDict(zip( keylist[ppp],wrealmittellist[ppp])))
        DlistFickdic.append(OrderedDict(zip( keylist[ppp],DlistFick[ppp])))
        DlistCrancdic.append(OrderedDict(zip( keylist[ppp],DlistCranc[ppp])))
        Feuchtedict.append(OrderedDict(zip( keylist[ppp],Feuchtelist[ppp])))

        wtrefdict.append(OrderedDict(zip( keylist[ppp],wtreflist[ppp])))
        wtdict.append(OrderedDict(zip( keylist[ppp],wtlist[ppp])))
        sratedict.append(OrderedDict(zip( keylist[ppp],sratelist[ppp])))
        sraterdict.append(OrderedDict(zip( keylist[ppp],sraterlist[ppp])))
        sraterrdict.append(OrderedDict(zip( keylist[ppp],sraterrlist[ppp])))
        ARDdict.append(OrderedDict(zip( keylist[ppp],ARDlist[ppp])))
        ARDrdict.append(OrderedDict(zip( keylist[ppp],ARDrlist[ppp])))
        if Film1.Mode=="MEOS":
            taudict.append(OrderedDict(zip( keylist[ppp],taulist[ppp])))
            expdict.append(OrderedDict(zip( keylist[ppp],explist[ppp])))
            Edict.append(OrderedDict(zip( keylist[ppp],Elist[ppp])))


    orderdkeys=keylist[0]
    wreallist=[]
    Feuchtelist=[]
    tveclist=[]
    wtlist=[]
    wtreflist=[]
    wrealmittellist=[]
    taulist=[]
    explist=[]
    Elist=[]
    ARDlist=[]
    ARDrlist=[]
    sratelist=[]
    sraterlist=[]
    sraterrlist=[]
    for ppp in range(len(Data)):
        wreallist.append(np.asarray([wrealdict[ppp][val] for i,val in enumerate(orderdkeys)]))



        tveclist.append(np.asarray([tvecdict[ppp][val] for i,val in enumerate(orderdkeys)]))
        wrealmittellist.append(np.asarray([wrealmitteldict[ppp][val] for i,val in enumerate(orderdkeys)]))

        Feuchtelist.append(np.asarray([Feuchtedict[ppp][val] for i,val in enumerate(orderdkeys)]))

        wtreflist.append(np.asarray([wtrefdict[ppp][val] for i,val in enumerate(orderdkeys)]))
        wtlist.append(np.asarray([wtdict[ppp][val] for i,val in enumerate(orderdkeys)]))
        sratelist.append(np.asarray([sratedict[ppp][val] for i,val in enumerate(orderdkeys)]))
        sraterlist.append(np.asarray([sraterdict[ppp][val] for i,val in enumerate(orderdkeys)]))
        sraterrlist.append(np.asarray([sraterrdict[ppp][val] for i,val in enumerate(orderdkeys)]))
        ARDlist.append(np.asarray([ARDdict[ppp][val] for i,val in enumerate(orderdkeys)]))
        ARDrlist.append(np.asarray([ARDrdict[ppp][val] for i,val in enumerate(orderdkeys)]))

        if Film1.Mode=="MEOS":
            taulist.append(np.asarray([taudict[ppp][val] for i,val in enumerate(orderdkeys)]))
            explist.append(np.asarray([expdict[ppp][val] for i,val in enumerate(orderdkeys)]))
            Elist.append(np.asarray([Edict[ppp][val] for i,val in enumerate(orderdkeys)]))



    ######

    THFaktor_fun,RH_fun,ww_fun=Film1.THFaktor_fun,Film1.RH_fun,Film1.ww_fun#Film1.PCSAFT(Film1.p,Film1.T,Film1.n)
    THFaktorNET_fun,THFaktorSA_fun=Film1.THFaktorNET_fun,Film1.THFaktorSA_fun
    #Film1.sim_boundary()
    #average and standardeviation
    Crancave,Crancstd=averagedict(DlistCrancdic)
    Fickave, Fickstd=averagedict(DlistFickdic)

    wave, wstd=averagedict(wrealdict)
    ARDrave, ARDrstd=averagedict(ARDrdict)
    sraterave, sraterstd=averagedict(sraterdict)
    sraterrave, sraterrstd=averagedict(sraterrdict)

    wmittelave, wmittelstd=averagedict(wrealmitteldict)
    #wmittelave, wmittelstd=wave, wstd
    Feuchteave, Feuchtestd=averagedict(Feuchtedict)

    #tvecave, tvecstd=LDF.averagelisttimeseries(tveclist,tveclist)
    wtave, wtstd=averagelisttimeseries(tveclist,wtlist)
    wtrefave, wtrefstd=averagelisttimeseries(tveclist,wtreflist)
    srateave, sratestd=averagelisttimeseries(tveclist,sratelist)
    ARDave, ARDstd=averagelisttimeseries(tveclist,ARDlist)
    #ARDrave, ARDrstd=averagelisttimeseries(tveclist,ARDrlist)
    
    mvecave, mcvecstd=averagelisttimeseries(tveclist,mveclist)

    if Film1.Mode=="MEOS":
        tauave, taustd=averagedict(taudict)
        expave, expstd=averagedict(expdict)
        Eave, Estd=averagedict(Edict)

    tarr=np.array(tveclist)
    nt=tarr.shape[2]
    for i in range(nt):
        tarr[:,:,i]=tarr[:,:,i]-tarr[:,:,0]
    tendn=tarr[:,:,-1]
    amin=np.argmin(tendn,axis=0)
    tveca=[tarr[val,i,:] for i,val in enumerate(amin)]
    #tvecave=tveclist[0]

    #wtave=wtlist[0]
    #wtrefave=wtreflist[0]
    nHum=int(DlistCranc[0].shape[0])
    Stefanave=np.zeros(DlistCranc[0].shape)
    Stefanstd=np.zeros(DlistCranc[0].shape)
    #wTH=(np.hstack((0,wave))[1:]+np.hstack((0,wave))[:-1])/2

    for jjj in range(nHum):

        Stefanave[jjj]=Fickave[jjj]/THFaktor_fun(wmittelave[jjj]) #changed wTH to wrealmittel
        Stefanstd[jjj]=Fickstd[jjj]/THFaktor_fun(wmittelave[jjj])
        #try 07042022______!!!
        #wpred=np.hstack((0,ww_fun(Feuchteave)))
        #wpredm=wpred[1:]*0.7+wpred[:-1]*0.3
        #Stefanave[jjj]=Fickave[jjj]/THFaktor_fun(wpredm[jjj]) #changed wTH to wrealmittel
        #Stefanstd[jjj]=Fickstd[jjj]/THFaktor_fun(wpredm[jjj])
        #try 07042022______!!!

    DlistStefan=[]
    for i,vali in enumerate(DlistFick):
        temp=np.asarray([valj for j, valj in enumerate(vali)]) if Film1.Mode=="MEOS" else np.asarray([valj/THFaktor_fun(wrealmittellist[i][j]) for j, valj in enumerate(vali)])
        DlistStefan.append(temp)
        DlistFick[i]= np.asarray([valj*THFaktor_fun(wrealmittellist[i][j]) for j, valj in enumerate(vali)]) if Film1.Mode=="MEOS" else DlistFick[i]
    Stefanave, Stefanstd=averagelist(DlistStefan)
    if Film1.Mode=="MEOS":
        Fickave, Fickstd=averagelist(DlistFick)

    #Read Literature Data
    #Tg
    polymer=Film1.molecule
    api=Film1.molecule2
    wPolyASD=Film1.wPolyASD
    wPolyASDstr=str(round(wPolyASD*100)/100) if not int(wPolyASD)==round(wPolyASD*100)/100 else str(int(wPolyASD))
    mol=polymer+"_"+api+"_"+wPolyASDstr if wPolyASD<0.99 else polymer

    THlist=[]
    THNETlist=[]
    THSAlist=[]
    for jjj in range(nHum):
        THlist.append(THFaktor_fun(wmittelave[jjj]))
        THNETlist.append(THFaktorNET_fun(wmittelave[jjj]))
        THSAlist.append(THFaktorSA_fun(wmittelave[jjj]))

    RH_fun=Film1.RH_fun
    wPCSAFT=np.linspace(0,np.max(wreal)*1.1,1000)
    wPCSAFT=np.linspace(0,1,1000)
    FeuchteSAFT=np.reshape(RH_fun(wPCSAFT),1000)
    wPCSAFTbar=1/2*(wPCSAFT[1:]+wPCSAFT[:-1])[1:-1]
    #wPCSAFTbar=(wPCSAFT[1:]*0.7+wPCSAFT[:-1]*0.3)[1:-1]

    #TH=np.ones(wPCSAFTbar.shape[0])
    #TH=THFaktor_fun(wPCSAFTbar)
    #2


    listofdicts={"w2/- ":wreallist,
    "w2ave/- ":[wave],
    "w2std/- ":[wstd],
    "RH/- ":Feuchtelist,
    "RHave/- ":[Feuchteave],
    "RHstd/- ":[Feuchtestd],
    "w2bar/- ":wrealmittellist,
    "w2barave/- ":[wmittelave],
    "w2barstd/- ":[wmittelstd],
    "DCranc/m/s^2 ":DlistCranc,
    "DCrancave/m/s^2 ":[Crancave],
    "DCrancstd/m/s^2 ":[Crancstd],
    "DFick/m/s^2 ":DlistFick,
    "DFickave/m/s^2 ":[Fickave],
    "DFickstd/m/s^2 ":[Fickstd],
    "DStefan/m/s^2 ":DlistStefan,
    "DStefanave/m/s^2 ":[Stefanave],
    "DFStefanstd/m/s^2 ":[Stefanstd],
    "ARD/-":ARDrlist,
    "ARDave/-":[ARDrave],
    "ARDstd/-":[ARDrstd],
    "srate/-":sraterlist,
    "srateave/-":[sraterave],
    "sratestd/-":[sraterstd],
    "sraterr/-":sraterrlist,
    "sraterrave/-":[sraterrave],
    "sraterrstd/-":[sraterrstd]}
    #WriteToExcel(listofdicts,["Tabelle1"],mol=mol)

    df1=pd.json_normalize(listofdicts, sep='_')
    df2=unnesting(df1,df1.columns.tolist(),axis=0)
    df3=df2.explode(df2.columns.tolist())
    excel_path=os.path.join(cwd,"LeastFit","Fit") if Fit else os.path.join(cwd,"LeastFit","Prediction")
    os.mkdir(excel_path) if not os.path.isdir(excel_path) else None
    df3.to_excel(os.path.join(excel_path,mol+".xlsx"))



    #RHaxis=np.asarray([0.10,0.30,0.45,0.60,0.75,0.90])
    #w_fun=InterpolatedUnivariateSpline(Feuchtereal,wreal)
    #waxis=w_fun(RHaxis)
    ### use only usual number of jumps if there are less or more
    # UsualNumberofJumps=6
    # wave=wave[:UsualNumberofJumps]
    # Feuchteave=Feuchteave[:UsualNumberofJumps]
    # Crancave=Crancave[:UsualNumberofJumps]
    # Fickave=Fickave[:UsualNumberofJumps]
    # Stefanave=Stefanave[:UsualNumberofJumps]




    waxis=(np.asarray(wave)*100).round(1)
    RHaxis=np.asarray(Feuchteave)*100
    #RH_fun=InterpolatedUnivariateSpline(wave,Feuchteave,k=2)



    wPCSAFT=np.linspace(0,np.max(wave)*1.1,1000)
    wPCSAFTbar=1/2*(wPCSAFT[1:]+wPCSAFT[:-1])[1:-1]
    TH=THFaktor_fun(wPCSAFTbar)
    THNET=THFaktorNET_fun(wPCSAFTbar)
    THSA=THFaktorSA_fun(wPCSAFTbar)

    #LDF.GammaPlot(TH, mol, wPCSAFTbar,Origin=Film1.Origin)
    Plotlist([wPCSAFTbar,wPCSAFTbar,wPCSAFTbar],[TH,THNET,THSA],[wPCSAFTbar,wPCSAFTbar,wPCSAFTbar],[TH,THNET,THSA],Origin=Film1.Origin,filename="TH_"+mol,xla="w2/-",yla="Gamma/-")
    # LDF.GammaPlot(THNET, mol+"NET", wPCSAFTbar)
    # LDF.GammaPlot(THSA, mol+"SA", wPCSAFTbar)

    wPCSAFT=np.linspace(0,1,1000)
    wPCSAFTbar=1/2*(wPCSAFT[1:]+wPCSAFT[:-1])[1:-1]
    #TH=THFaktor_fun(wPCSAFTbar)
    # THNET=THFaktorNET_fun(wPCSAFTbar)
    # THSA=THFaktorSA_fun(wPCSAFTbar)
    #LDF.GammaPlot(TH,mol+"zoom", wPCSAFTbar,Origin=Film1.Origin)
    #Plotlist([wPCSAFTbar],[TH],[wPCSAFTbar],[TH],Origin=Film1.Origin,filename="THzoom_"+mol,xla="w2/-",yla="Gamma/-")
    # LDF.GammaPlot(THNET, mol+"NET"+"zoom", wPCSAFTbar)
    # LDF.GammaPlot(THSA, mol+"SA"+"zoom", wPCSAFTbar)



    #LDF.IsothermPlot(Feuchteave,Feuchtestd,wave,wstd,mol,FeuchteSAFT,wPCSAFT,Origin=Film1.Origin)
    Plotlist([np.asarray(Feuchteave)],[np.asarray(wave)],[np.asarray(FeuchteSAFT)],[np.asarray(wPCSAFT)],Origin=Film1.Origin,filename="Isotherm_"+mol,xla="RH/-",yla="w2/-")
    for i,val in enumerate(tveclist):
        #LDF.SorpionStepsPlot(tveclist[i],wtlist[i],wtreflist[i],i,mol,wtrefstd=np.zeros_like(wtrefstd),conti=0,Origin=Film1.Origin,des=Film1.des)
        Plotlist(tveclist[i],wtreflist[i],tveclist[i],wtlist[i],Origin=Film1.Origin,filename="Kinetics_"+str(i)+"_"+mol,xla="t/min",yla="w2/-")
        #Plotlist(tveclist[i],sratelist[i],tveclist[i],sratelist[i],Origin=Film1.Origin,filename="srate_"+str(i)+"_"+mol,xla="t/min",yla="srate/min")
        #Plotlist(tveclist[i],ARDlist[i],tveclist[i],ARDlist[i],Origin=Film1.Origin,filename="ARD_"+str(i)+"_"+mol,xla="t/min",yla="ARD/-")
        #LDF.SorpionStepsPlot(tveclist[i],mveclist[i],mveclist[i],i,mol,wtrefstd=[0],conti=1)
    #LDF.SorpionStepsPlot(tveca,wtave,wtrefave,"ave",mol,wtrefstd=wtrefstd,conti=0,Origin=Film1.Origin,des=Film1.des)
    Plotlist(tveca,wtrefave,tveca,wtave,wtrefstd,wtrefstd,Origin=Film1.Origin,filename="Kinetics_ave_"+mol,xla="t/min",yla="w2/-")
    Plotlist(tveca,srateave,tveca,srateave,sratestd,sratestd,Origin=Film1.Origin,filename="srate_ave_"+mol,xla="t/min",yla="srate/min")
    #Plotlist(tveca,ARDave,tveca,ARDrave,ARDstd,ARDstd,Origin=Film1.Origin,filename="ARD_ave_"+mol,xla="t/min",yla="ARD/-")

    #LDF.DiffcoeffPlot(Crancave, Crancstd, wave, wstd, mol, "Cranc",waxis,RHaxis,Origin=Film1.Origin) if Fit else None
    #Plotlist([wave],[Crancave],[wave],[Crancave],[wstd],[Crancstd],Origin=Film1.Origin,filename="Cranc_"+mol,xla="w2/-",yla="D/m/s^2") if Fit else None
    #LDF.DiffcoeffPlot(Fickave, Fickstd, wave, wstd, mol, "Fick", waxis,RHaxis,Origin=Film1.Origin) if Fit else None
    #Plotlist([wave],[Fickave],[wave],[Fickave],[wstd],[Fickstd],Origin=Film1.Origin,filename="Fick_"+mol,xla="w2/-",yla="DV/m/s^2") if Fit else None
    #LDF.DiffcoeffPlot(Stefanave, Stefanstd, wave, wstd, mol, "Stefan",waxis,RHaxis,Origin=Film1.Origin) if Fit else None
    #Plotlist([wave],[Stefanave],[wave],[Stefanave],[wstd],[Stefanstd],Origin=Film1.Origin,filename="Stefan_"+mol,xla="w2/-",yla="DS/m/s^2") if Fit else None
    Plotlist([wave,wave,wave],[Crancave,Fickave,Stefanave],[wave,wave,wave],[Crancave,Fickave,Stefanave],[wstd,wstd,wstd],[Crancstd,Fickstd,Stefanstd],Origin=Film1.Origin,filename="D_"+mol,xla="w2/-",yla="DS/m/s^2") if Fit else None
    if Film1.Mode=="MEOS":
        #LDF.DiffcoeffPlot(tauave, taustd, wave, wstd, mol, "tau",waxis,RHaxis,Origin=Film1.Origin) if Fit else None
        Plotlist([wave],[tauave],[wave],[tauave],[wstd],[taustd],Origin=Film1.Origin,filename="Tau_"+mol,xla="w2/-",yla="tau/s") if Fit else None
        #LDF.DiffcoeffPlot(np.fmax(expave,1E-10), np.fmax(expstd,1E-10), wave, wstd, mol, "exp",waxis,RHaxis,Origin=Film1.Origin) if Fit else None
        Plotlist([wave],[np.fmax(expave,1E-10)],[wave],[np.fmax(expave,1E-10)],[wstd],[np.zeros_like(expave)],Origin=Film1.Origin,filename="EXP_"+mol,xla="w2/-",yla="Lambda/-") if Fit else None
        #LDF.DiffcoeffPlot(Eave, Estd, wave, wstd, mol, "E",waxis,RHaxis,Origin=Film1.Origin) if Fit else None
        Plotlist([wave],[Eave],[wave],[Eave],[wstd],[Estd],Origin=Film1.Origin,filename="E_"+mol,xla="w2/-",yla="E/Pa") if Fit else None
    plt.show()

if __name__=="__main__":
    window1=window()
    Film1=Film(window1.result,window1.components,window1.wPolyASD)
    #FitDiffusionData(Data,water,Matrix)
    FitDiffusionData(Film1,window1.result,Fit=False)
0