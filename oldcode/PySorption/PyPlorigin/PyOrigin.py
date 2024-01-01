import win32com.client
from os.path import join
from os import getcwd
import numpy as np
import itertools
import pandas as pd
import re
nane=np.asarray([[np.nan,np.nan,np.nan],[np.nan,np.nan,np.nan]]).T

def MusterDictOrigin():
    xData={"values":0,"unit":"-","comment":"-"}
    yData={"values":0,"unit":"-","comment":"-"}
    xSim={"values":0,"unit":"-","comment":"-"}
    ySim={"values":0,"unit":"-","comment":"-"}
    xErr={"values":0,"unit":"-","comment":"-"}
    yErr={"values":0,"unit":"-","comment":"-"}
    return xData,yData,xSim,ySim,xErr,yErr
#ERROR!
def TernaryPlot(SLE=nane,LLE=[nane,nane],VLE=[nane,nane,nane],Tg=nane,SLEexp=nane,LLEexp=nane,VLEexp=[nane,nane],Tgexp=nane,molecules=["polymer","solvent","api"],G=nane,Spinodal=[nane,nane]):
    """Plots a ternary Phase Diagramm with two Phase LLE and liquid phase VLE and SLE composotions.
    The composition Data has to be supplied in the form np.array([w1,w2,w3]) and is therefore of shape (3,X)
    where X is the number of Datapoints. If no Data for the keywords is supplied NAN-Values will be filled
    and if all values of a workbook are nan values it is deleted. LLE are a list of the two phasses [(3,X),(3,X)]
    The VLE is a list of isohumidity values [(3,X),...,(3,X)] where the length depends on the humidity
    values supplied. Measurements, Tg and SLE compositions are supplied as (3,X) """

    # My Code usually suplies polymer, solvent,api while the Vorlage is of form api,solvent,polymer
    SLE=SLE[::-1,:]
    LLE=[val[::-1,:] for i,val in enumerate(LLE)]
    Spinodal=[val[::-1,:] for i,val in enumerate(Spinodal)]
    #G=[G[0][::-1,:],G[1]]
    VLE=[val[::-1,:] for i,val in enumerate(VLE)]
    Tg=Tg[::-1,:]
    SLEexp=SLEexp[::-1,:]
    LLEexp=LLEexp[::-1,:]
    VLEexp=[val[::-1,:] for i,val in enumerate(VLEexp)]
    Tgexp=Tgexp[::-1,:]
    ##
    filename="_".join(molecules)
    app = win32com.client.Dispatch("Origin.ApplicationSi")
    app.Execute ("sec -poc 10")
    string2=join(getcwd(),"OriginVorlage")
    string3="Origin_Ternär_Vorlage.opj"
    string4=join(string2,string3)
    app.Load(string4)
    app.Execute ("sec -poc 10")
    cwd = r'%s' % getcwd().replace('\\','/')
    pageName=[app.WorksheetPages[pg].name for pg,value in enumerate(app.WorksheetPages)]
    SLEFlaeche=np.asarray([[1,0,0],SLE[:,0]]).T
    TgFlaeche=np.asarray([[0,0,1],[1,0,0],Tg[:,0]]).T
    SLE=np.hstack((SLE,SLEFlaeche))
    Tg=np.hstack((Tg,TgFlaeche))
    dic={"SLE":SLE,
         "LLE":LLE,
         "VLE":VLE,
         "Tg":Tg,
         "SLEMessung":SLEexp,
         "LLEMessung":LLEexp,
         "VLEMessung":VLEexp,
         "TgMessung":Tgexp,
         "Gibbs": G,
         "Spinodal":Spinodal}
    #[3,:],[2,:],[1,:]
    # consider Rolling
    #VLE with pointing and phi should be as accurate
    w1list=[]
    w2list=[]
    w3list=[]
    SLE1=[]
    SLE2=[]
    SLE3=[]

    for pg,value in enumerate(pageName):

        if pageName[pg]=="LLE":
            for i,vali in enumerate(dic[pageName[pg]][0][0,:]):
                w1list.append(LLE[0][0,i])
                w1list.append(LLE[1][0,i])
                w1list.append("-")
                w2list.append(LLE[0][1,i])
                w2list.append(LLE[1][1,i])
                w2list.append("-")
                w3list.append(LLE[0][2,i])
                w3list.append(LLE[1][2,i])
                w3list.append("-")
            app.PutWorksheet(pageName[pg], w1list, 0, 3)
            app.Execute ("sec -poc 10")
            app.PutWorksheet(pageName[pg], w2list, 0, 4)
            app.Execute ("sec -poc 10")
            app.PutWorksheet(pageName[pg], w3list, 0, 5)
            app.Execute ("sec -poc 10")
            for j,valj in enumerate(dic[pageName[pg]][0][:,0]):
                app.PutWorksheet(pageName[pg], dic[pageName[pg]][0][j,:].tolist(), 0, j)
                app.Execute ("sec -poc 10")
            for j,valj in enumerate(dic[pageName[pg]][1][:,0]):
                app.PutWorksheet(pageName[pg], dic[pageName[pg]][1][j,:].tolist(), 0, 6+j)
                app.Execute ("sec -poc 10")
        elif pageName[pg]=="Spinodal":
            for j,valj in enumerate(dic[pageName[pg]][0][:,0]):
                app.PutWorksheet(pageName[pg], dic[pageName[pg]][0][j,:].tolist(), 0, j)
                app.Execute ("sec -poc 10")
            for j,valj in enumerate(dic[pageName[pg]][1][:,0]):
                app.PutWorksheet(pageName[pg], dic[pageName[pg]][1][j,:].tolist(), 0, 3+j)
                app.Execute ("sec -poc 10")
        elif pageName[pg]=="VLE":
            for i,vali in enumerate(dic[pageName[pg]]):
                for j,valj in enumerate(dic[pageName[pg]][0][:,0]):
                    app.PutWorksheet(pageName[pg], dic[pageName[pg]][i][j,:].tolist(), 0, 4*i+j)
                    app.Execute ("sec -poc 10")
        elif pageName[pg]=="VLEMessung":
            for i,vali in enumerate(dic[pageName[pg]]):
                for j,valj in enumerate(dic[pageName[pg]][0][:,0]):
                    app.PutWorksheet(pageName[pg], dic[pageName[pg]][i][j,:].tolist(), 0, 4*i+j)
                    app.Execute ("sec -poc 10")
        elif pageName[pg]=="SLE":
            for i,val in enumerate(dic[pageName[pg]][:,0]):
                app.PutWorksheet(pageName[pg], dic[pageName[pg]][i,:].tolist(), 0, i)
                app.Execute ("sec -poc 10")

            # for i,vali in enumerate(dic[pageName[pg]][0,:]):
            #     SLE1.append(SLE[0,i])
            #     SLE1.append(1)
            #     SLE1.append("-")
            #     SLE2.append(SLE[1,i])
            #     SLE2.append(0)
            #     SLE2.append("-")
            #     SLE3.append(SLE[2,i])
            #     SLE3.append(0)
            #     SLE3.append("-")
            #SLE1=[dic[pageName[pg]][0,0],1,"-",dic[pageName[pg]][0,-1],1,"-",dic[pageName[pg]][0,0],dic[pageName[pg]][0,-1]]
            #SLE2=[dic[pageName[pg]][1,0],0,"-",dic[pageName[pg]][1,-1],0,"-",dic[pageName[pg]][1,0],dic[pageName[pg]][1,-1]]
            #SLE3=[dic[pageName[pg]][2,0],0,"-",dic[pageName[pg]][2,-1],0,"-",dic[pageName[pg]][2,0],dic[pageName[pg]][2,-1]]
            # app.PutWorksheet(pageName[pg], SLE1, 0, 3)
            # app.Execute ("sec -poc 10")
            # app.PutWorksheet(pageName[pg], SLE2, 0, 4)
            # app.Execute ("sec -poc 10")
            # app.PutWorksheet(pageName[pg], SLE3, 0, 5)
            # app.Execute ("sec -poc 10")
        elif pageName[pg]=="Tg":
            for i,val in enumerate(dic[pageName[pg]][:,0]):
                app.PutWorksheet(pageName[pg], dic[pageName[pg]][i,:].tolist(), 0, i)
                app.Execute ("sec -poc 10")
            #Tg1=[dic[pageName[pg]][0,0],0,"-",dic[pageName[pg]][0,-1],0,"-",dic[pageName[pg]][0,0],dic[pageName[pg]][0,-1]]
            #Tg2=[dic[pageName[pg]][1,0],0,"-",dic[pageName[pg]][1,-1],0,"-",dic[pageName[pg]][1,0],dic[pageName[pg]][1,-1]]
            #Tg3=[dic[pageName[pg]][2,0],1,"-",dic[pageName[pg]][2,-1],1,"-",dic[pageName[pg]][2,0],dic[pageName[pg]][2,-1]]
            #app.PutWorksheet(pageName[pg], Tg1, 0, 3)
            #app.Execute ("sec -poc 10")
            #app.PutWorksheet(pageName[pg], Tg2, 0, 4)
            #app.Execute ("sec -poc 10")
            #app.PutWorksheet(pageName[pg], Tg3, 0, 5)
            #app.Execute ("sec -poc 10")
        elif pageName[pg]=="Gibbs":
            for i,val in enumerate(dic[pageName[pg]][:,0]):
                app.PutWorksheet(pageName[pg], dic[pageName[pg]][i,:].tolist(), 0, i)
                app.Execute ("sec -poc 10")


        else:
            for i,val in enumerate(dic[pageName[pg]][:,0]):
                app.PutWorksheet(pageName[pg], dic[pageName[pg]][i,:].tolist(), 0, i)
                app.Execute ("sec -poc 10")
        try:
            if  pd.isnull(dic[pageName[pg]]).any():
                app.Execute ("sec -poc 10")
                app.Execute("layer -d "+pageName[pg])
                app.DestroyPage(pageName[pg])
                app.Execute ("sec -poc 10")
        except:
            print("Isnan still makes problems")



    cwd = r'%s' % string2.replace('\\','/')
    app.Execute ("sec -poc 10")
    app.Execute ('expGraph type:=jpg path:="'+cwd+'" export:=project overwrite:=Replace filename:="'+filename+'";')
    #strPathName = app.GetSaveAsFilename(filename, "Project files (*.OPJ),*.OPJ", 0)
    app.Execute ("sec -poc 10")
    app.Execute('save -DIX %Y'+filename+'.opju;')
    #fig,ax=plt.subplots()
    #im=plt.imread(cwd+"/"+filename+".jpg")
    #ax.imshow(im)
    #ax.axis("off")
    #plt.show()
    #app.Execute ("sec -poc 10")
    app.Exit()



def Plot(x,y,xs,ys,xe,ye,xlabel,ylabel,filename,Vorlage="Origin.opj"):
    app = win32com.client.Dispatch("Origin.ApplicationSi")
    app.Execute ("sec -poc 10")
    xData,yData,xSim,ySim,xErr,yErr=x["values"],y["values"],xs["values"],ys["values"],xe["values"],ye["values"]
    xu,yu,xsu,ysu,xeu,yeu=x["unit"],y["unit"],xs["unit"],ys["unit"],xe["unit"],ye["unit"]
    xc,yc,xsc,ysc,xec,yec=x["comment"],y["comment"],xs["comment"],ys["comment"],xe["comment"],ye["comment"]

    string2=join(getcwd(),"OriginVorlage")
    string3=Vorlage

    string4=join(string2,string3)
    app.Load(string4)
    print("Origin template was successfully opened.")
    app.Execute ("sec -poc 10")

    pageName=[]
    pageNumber=[]
    before=len(app.WorksheetPages)
    if len(app.WorksheetPages)<len(xData):
        for i in range(len(app.WorksheetPages),len(xData)):
            app.Execute ("sec -poc 10")
            source="Sheet"+str(i+1)
            source2="Reihe"+str(i+1)
            app.Execute ("newbook name:="+source+" sheet:=1")
        

    for pg,value in enumerate(app.WorksheetPages):
        pageName.append(app.WorksheetPages[pg].name)
        pageNumber.append(int(re.findall(r'\d+', pageName[pg])[0]))
    pageName=np.asarray(pageName)[np.argsort(pageNumber)]
    #pageName=np.sort(pageName)
    app.Execute ("sec -poc 10")
    pgl=len(pageName)
    for pg,value in enumerate(pageName):
        x=np.asarray(np.nan) if pg>=len(xData) else xData[pg]
        y=np.asarray(np.nan) if pg>=len(yData) else yData[pg]
        xs=np.asarray(np.nan) if pg>=len(xSim) else xSim[pg]
        ys=np.asarray(np.nan) if pg>=len(ySim) else ySim[pg]
        xe=np.asarray(np.nan) if pg>=len(xErr) else xErr[pg]
        ye=np.asarray(np.nan) if pg>=len(yErr) else yErr[pg]

        if  pd.isnull(x).any() and pd.isnull(y).any()  and pd.isnull(xs).any()  and pd.isnull(ys).any()  and pd.isnull(xe).any()  and pd.isnull(ye).any() :
            app.Execute ("sec -poc 10")
            app.Execute("layer -d "+pageName[pg])
            app.DestroyPage(pageName[pg])
            app.Execute ("sec -poc 10")
        else:
            app.Execute ("sec -poc 10")
            app.Execute ("win -a "+pageName[pg])
            app.Execute ("sec -poc 10")
            app.PutWorksheet(pageName[pg], x.tolist(), 0, 0)
            app.Execute ("sec -poc 10")
            app.Execute ('col(1)[C]$ ='+xc+';')
            app.Execute ("sec -poc 10")
            app.Execute ('col(1)[U]$ ='+xu+';')
            app.Execute ("sec -poc 10")
            app.PutWorksheet(pageName[pg], y.tolist(), 0, 1)
            app.Execute ("sec -poc 10")
            app.Execute ('col(2)[C]$ ='+yc+';')
            app.Execute ("sec -poc 10")
            app.Execute ('col(2)[U]$ = '+yu+';')
            app.Execute ("sec -poc 10")
            app.PutWorksheet(pageName[pg], xs.tolist(), 0, 4)
            app.Execute ("sec -poc 10")
            app.Execute ('col(5)[C]$ ='+xsc+';')
            app.Execute ("sec -poc 10")
            app.Execute ('col(5)[U]$ = '+xsu+';')
            app.Execute ("sec -poc 10")
            app.PutWorksheet(pageName[pg], ys.tolist(), 0, 5)
            app.Execute ("sec -poc 10")
            app.Execute ('col(6)[C]$ ='+ysc+';')
            app.Execute ("sec -poc 10")
            app.Execute ('col(6)[U]$ = '+ysu+';')
            app.Execute ("sec -poc 10")
            app.PutWorksheet(pageName[pg], xe.tolist(), 0, 2)
            app.Execute ("sec -poc 10")
            app.Execute ('col(3)[C]$ ='+xec+';')
            app.Execute ("sec -poc 10")
            app.Execute ('col(3)[U]$ = '+xeu+';')
            app.Execute ("sec -poc 10")
            app.PutWorksheet(pageName[pg], ye.tolist(), 0, 3)
            app.Execute ("sec -poc 10")
            app.Execute ('col(4)[C]$ ='+yec+';')
            app.Execute ("sec -poc 10")
            app.Execute ('col(4)[U]$ = '+yeu+';')
            app.Execute ("sec -poc 10")
            app.Execute ("page.active:=1")
            app.Execute ("sec -poc 10")            
            app.Execute ("plotxy iy:=["+pageName[pg]+"]1!(1,2) o:=[Graph1]1 plot:=201;") if pg>=before else None
            app.Execute ("sec -poc 10")
            app.Execute ("plotxy iy:=["+pageName[pg]+"]1!(5,6) o:=[Graph1]1 plot:=200;") if pg>=before else None
            app.Execute ("sec -poc 10")




    try:
        xmax=np.max(list(itertools.chain(*xData)))
        x0=np.min(list(itertools.chain(*xData)))
        ymax=np.max(list(itertools.chain(*yData)))
        y0=np.min(list(itertools.chain(*yData)))
    except:
        xmax=np.max(list(xData))
        x0=np.min(list(xData))
        ymax=np.max(list(yData))
        y0=np.min(list(yData))


    xmag=10**(np.floor(np.log10(abs(xmax-x0))))
    ymag=10**(np.floor(np.log10(abs(ymax-y0))))
    if xmag==0:
        xmag=1
    if ymag==0:
        ymag=1
    def myround(x, base=5):
        return base*0.5 * np.round(x/(base*0.5))
    xmax=myround(xmax,xmag)#+xemax*2
    ymax=myround(ymax,ymag)#+yemax*2
    x0=myround(x0,xmag)#-xe0*2
    y0=myround(y0,ymag)#-ye0*2
    app.Execute ("sec -poc 10")
    app.Execute ("win -a "+"Graph1")

    app.Execute("laysetscale igp:=Graph1 layer:=1 axis:=x from:="+str(x0)+" to:="+str(xmax))
    app.Execute ("sec -poc 10")
    app.Execute("laysetscale igp:=Graph1 layer:=1 axis:=y from:="+str(y0)+" to:="+str(ymax))
    app.Execute ("sec -poc 10")
    # try:
    #     app.Execute('label -n XB "'+ xlabel+'"')
    #     app.Execute ("sec -poc 10")
    #     app.Execute('label -n YL "' +ylabel+'"')
    #     app.Execute ("sec -poc 10")
    # except:
    #     pass
    app.Execute('save -DIX %Y'+filename+'.opju;')
    app.Execute ("sec -poc 10")
    print("Done! Origin was safely closed.")
    app.Exit()




def Contour(xyz=[nane,nane,nane],filename="filename"):
    app = win32com.client.Dispatch("Origin.ApplicationSi")
    app.Execute ("sec -poc 10")
    string2=join(getcwd(),"OriginVorlage")
    string3="Contour.opj"
    string4=join(string2,string3)
    app.Load(string4)
    print("Origin template was successfully opened.")
    app.Execute ("sec -poc 10")
    cwd = r'%s' % getcwd().replace('\\','/')
    pageName=[app.WorksheetPages[pg].name for pg,value in enumerate(app.WorksheetPages)]
    dic={"xyz":xyz}
    for pg,value in enumerate(pageName):
        if pageName[pg]=="xyz":
            for i,vali in enumerate(dic[pageName[pg]]):
                for j,valj in enumerate(dic[pageName[pg]][0][:,0]):
                    app.PutWorksheet(pageName[pg], dic[pageName[pg]][i][j,:].tolist(), 0, 4*i+j)
                    app.Execute ("sec -poc 10")

    app.Execute ("sec -poc 10")
    app.Execute('save -DIX %Y'+filename+'.opju;')
    app.Execute ("sec -poc 10")
    app.Exit()
    print("Done! Origin was safely closed.")



if __name__=="__main__":
    #example for Plot

    #generate two sets of data with experimental data, xerror, yerror and simulation
    import random
    n=10
    xdata1,xdata2=np.linspace(0,10,n),np.linspace(0,10,n)
    xsim1,xsim2=xdata1,xdata2
    ysim1,ysim2=xdata1*4,xdata2**2
    ydata1,ydata2=ysim1+(np.random.rand(n)-0.5)*10,ysim2+(np.random.rand(n)-0.5)*10
    xerr1,xerr2=np.random.rand(n)*0.5,np.random.rand(n)*0.5
    yerr1,yerr2=np.random.rand(n)*5,np.random.rand(n)*5
    xvec=[xdata1,xdata2] #Provide datasets as lists.
    yvec=[ydata1,ydata2]
    xsvec=[xsim1,xsim2]
    ysvec=[ysim1,ysim2]
    xevec=[xerr1,xerr2]
    yevec=[yerr1,yerr2]
    # Use Muster dictionary to fill in data
    x,y,xs,ys,xe,ye=MusterDictOrigin()
    # Fill in data
    x["values"],y["values"]=xvec,yvec
    xs["values"],ys["values"]=xsvec,ysvec
    xe["values"],ye["values"]=xevec,yevec
    # Fill in units if necessary else default is -
    x["unit"],y["unit"]="min","mg"
    xs["unit"],ys["unit"]="min","mg"
    xe["unit"],ye["unit"]="min","mg"
    #Fill in comments if necessary else default is -
    x["comment"],y["comment"]="time","mass reading"
    xs["comment"],ys["comment"]="time","mass simulation"
    xe["comment"],ye["comment"]="t error","mass error"
    # Choose x and y labels
    xlabel="t/min"
    ylabel="m/%(CRLF)mg"
    #give it a name name
    filename="example"
    Vorlage="Origin.opj" #Name der orgin file in: cwd+"OriginVorlage" Ordner
    Vorlage="Solubilities_Vorlage.opj"
    Plot(x,y,xs,ys,xe,ye,xlabel,ylabel,filename,Vorlage=Vorlage)

    # w2=np.linspace(0,1,100)
    # w1=(0.2)*(1-w2)
    # w3=(0.8)*(1-w2)
    # SLE=np.asarray([w1,w2,w3])

    # Tg=np.asarray([w3,w1,w2])
    # LLE=[SLE,Tg]
    # VLE=[SLE,Tg,np.asarray([w2,w3,w1]),SLE,Tg,np.asarray([w2,w3,w1])]
    # SLEexp=np.asarray([[np.nan,np.nan,np.nan],[np.nan,np.nan,np.nan]]).T
    # LLEexp=SLEexp
    # VLEexp=SLEexp
    # Tgexp=SLEexp
    #TernaryPlot(SLE,LLE,VLE,Tg,SLEexp,LLEexp,VLEexp,Tgexp)
    #RadialContour(rho=[Tg],sigma=[Tg])
