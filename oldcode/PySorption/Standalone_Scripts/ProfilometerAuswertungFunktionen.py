import numpy as np
from scipy.ndimage import median_filter
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import matplotlib.cm as cm
from os.path import join
import os
from scipy.stats import zscore
import csv
from itertools import zip_longest
import pandas as pd
from mpl_toolkits.mplot3d import axes3d

def MusterDictOrigin():
    xData={"values":5,"unit":"-","comment":"-"}
    yData={"values":5,"unit":"-","comment":"-"}
    xSim={"values":5,"unit":"-","comment":"-"}
    ySim={"values":5,"unit":"-","comment":"-"}
    xErr={"values":5,"unit":"-","comment":"-"}
    yErr={"values":5,"unit":"-","comment":"-"}
    return xData,yData,xSim,ySim,xErr,yErr

def ReadData(path1):
    fobj1 = open(path1, errors='ignore').read()
    fobj2=fobj1.replace("\t","")
    fobj3=fobj2.replace("\n","")
    l1=fobj3.find('Sclen')+6
    l2=fobj3.find('Stitch')
    startidx1=fobj3.rfind("SCALED DATA:")+len("SCALED DATA:\n\n")
    endidx1=fobj3.rfind("DEKTAK")
    Data1=fobj3[startidx1:endidx1]
    Data1=Data1+("     ")
    Datared=[]
    j=0
    Length=float(fobj3[l1:l2])
    for i in range(len(Data1)):
        if i>j:
            if Data1[i]!=" ":
                k=0
                j=i
                if j<len(Data1):
                    while k==0:
                        j=j+1
                        if Data1[j]==" ":
                                k=1
                                Datared.append(float(Data1[i:j]))
    Datared=np.asarray(Datared)*1000/930 #The profilometer was wrongly calibrated the calibration standard of thickness 1000µm was measured 930 µm. Unfortunatly, Nobody at the Micro-Nano Chair cares to calibrate this thing...  

    return Datared, Length

def Jumpdetection(data,windowlength,plotrawdata,shift1):
    Dres=data
    data = (data - data.mean())/data.std()

    template = np.concatenate([[-1]*windowlength,[1]*windowlength])
    nData=data.shape[0]
    corr_res =np.correlate(data, template,mode='same')
    corr_res=np.gradient(data)
    #dcorr_res=np.gradient(corr_res)
    corL=corr_res[shift1:int(nData/2)]
    corR=corr_res[int(nData/2):-shift1]
    #corL=ddata[shift1:int(nData/2)]
    #corR=ddata[int(nData/2):-shift1]

    indL=np.argmax(corL)+shift1
    indR=int(nData/2)+np.argmin(corR)
    if plotrawdata==1:
        fig3, ax3=plt.subplots()
        #ax3.plot(corr_res)
        ax3.plot(Dres,label='rawData')
        ax3.plot(corr_res*10,label='Violations')
        ax3.plot([indL,indR],corr_res[[indL,indR]]*10,'rx',label='Detected')
        ax3.set_xlabel('Data Point[-]')
        ax3.set_ylabel('Convolve[-]')
        ax3.legend()
    return indL,indR

def Outlierdetection(xdata,ydata):
    threshold=2.5
    zdata=np.abs(zscore(ydata))
    xres=xdata[zdata<threshold]
    yres=ydata[zdata<threshold]
    return xres,yres

def Jumpdetection2(indL,indR,D1,plotrawdata,windowlength,shift2,name):
    max_indx=indL
    min_indx=indR
    nData=D1.shape[0]
    x=np.linspace(1,nData,nData)
    
    #xnL=x[max_indx-n:max_indx+n]
    #xnR=x[min_indx-n:min_indx+n]
    #LedgeL=D1[max_indx-n:max_indx+n]
    #LedgeR=D1[min_indx-n:min_indx+n]
    #peL=np.polyfit(xnL,LedgeL,1)
    #peR=np.polyfit(xnR,LedgeR,1)
    #veL=np.poly1d(peL)
    #veR=np.poly1d(peR)
    #xL=x[0:max_indx]
    #xR=x[min_indx:-1]
    #indpL=np.argmin(abs(LD1-veL(xL)))
    #indpR=np.argmin(abs(RD1-veR(xR)))+min_indx
    indpL=int(max_indx-windowlength/2)
    indpR=int(min_indx+windowlength/2)
    # indpL=0+100
    # indpR=nData-100
    
    xpL=x[shift2:indpL]
    xpR=x[indpR:-1]
    PL=D1[shift2:indpL]
    PR=D1[indpR:-1]
    xpL,PL=Outlierdetection(xpL,PL)
    xpR,PR=Outlierdetection(xpR,PR)
    PLmed=median_filter(PL,101,mode='nearest')
    PRmed=median_filter(PR,101,mode='nearest')
    Pmed=np.hstack((PLmed,PRmed))
    xp=np.hstack((xpL,xpR))
    P=np.hstack((PL,PR))
    sigma=np.fmax(((Pmed-P)**2/(P.shape[0]-1))**(1/2),1E-1)
    w=1/sigma    
    Poly2=np.polyfit(xp,P,3,w=w)
    vpoly2=np.poly1d(Poly2)    
    correction=vpoly2(x)
    #correction=np.zeros_like(D1)
    Dres=D1-correction
    string=["correction","correctedData"]
    DatatoCSV=[correction,Dres]
    #WriteResultsCSV(string,DatatoCSV,"correction"+name)
    if plotrawdata==1:
        fig2, ax2=plt.subplots()
        ax2.plot(D1,label='rawData')
        ax2.plot(Dres,label='correctedData')
        ax2.plot(correction,label='correction')
        ax2.set_xlabel('Data Points[-]')
        ax2.set_ylabel('Height[nm]')
        ax2.legend()
    return Dres,indpL,indpR

def Dektak(Data,Len,m,turn,windowlength,plotrawdata,shift1,shift2,name,rhoreal):
    Datafilt=median_filter(Data,301,mode='nearest')
    indL, indR=Jumpdetection(Datafilt,windowlength,plotrawdata,shift1) 
    Dres,indpL,indpR=Jumpdetection2(indL,indR,Data,plotrawdata,windowlength,shift2,name)
    
    
    Dres[0:indpL]=0
    Dres[indpR:-1]=0
    correctedData=Dres
    nData=Dres.shape[0]
    nR=int(nData)
    Reduction=5
    l = range(nR-1)
    spacing=l[0::Reduction]
    nR=int(nR/Reduction)
    
    Dres=Dres[spacing]
    thr=2000
    Dresr1=Dres>thr
    Dresr2=Dres[::-1]>thr
    correctedData1=correctedData>thr
    correctedData2=correctedData[::-1]>thr
   
    for i in range(nR-1):
        if Dresr1[i]==False and Dresr1[i+1]==True:
            idx1=i
            break
    for i in range(nR-1):        
        if Dresr2[i]==False and Dresr2[i+1]==True:
            idx2=nR-i
            break
    for i in range(nData-1):
        if correctedData1[i]==False and correctedData1[i+1]==True:
            idxx1=i
            break
    for i in range(nData-1):        
        if correctedData2[i]==False and correctedData2[i+1]==True:
            idxx2=nData-i
            break
    
    
    if (idx1+idx2)%2==1:
        idx2+=1
    
    D1=Dres[idx1:idx2]
    #correctedData=correctedData[idxx1:idxx2]
    nD1=D1.shape[0]
    LD1=np.fmax(Dres[idx1:int((idx2+idx1)/2)],0)
    RD1=np.fmax(Dres[int((idx2+idx1)/2):idx2],0)
    LD1=LD1[::-1]
    D=Len*1000*nD1/nR

    #D=14547048
    R=D/2
    nR1=int(nD1/2)
    rvecL=np.linspace(0,R,nR1)
    rvecR=np.linspace(0,R,nR1)
    rvecmL=np.ones(nR1-1)
    rvecmR=np.ones(nR1-1)
    LD1m=np.ones(nR1-1)
    RD1m=np.ones(nR1-1)
    
    #Data with x values to CSV
    rvec=np.linspace(-Len/2*1E3,Len/2*1E3,np.shape(correctedData)[0])
    string=["Scanlength","correctedData",'Lmittel']
    #rvecwrite=np.linspace(-Len/2*1E3,Len/2*1E3,np.shape(correctedData)[0])
    dpi=96
    fs=(4604/(8*dpi),3602/(8*dpi))
    fig3, ax3=plt.subplots(dpi=dpi,figsize=fs)
    ax3.plot(rvec/1E6,correctedData/1E3)
    Lmittel=(m*1E-6)/rhoreal/(np.pi/4*(D*1E-9)**2)*1E9
    Lmittelvec=np.zeros(np.shape(correctedData)[0])
    Lmittelvec[idxx1:idxx2]=Lmittel*np.ones(int(idxx2-idxx1))
    ax3.plot(rvec/1E6,Lmittelvec/1E3)
    ax3.set_xlabel('ScanLength[mm]')
    ax3.set_ylabel('Height[$\mu m$]')
    #ax3.legend()
    #fig3.savefig(join(os.getcwd(),"Plots",name+"side.jpeg"))
    DatatoCSV=[rvec,correctedData,Lmittelvec]
    #WriteResultsCSV(string,DatatoCSV,"xdata"+name)
    ThicknessPlot(rvec/1E6,correctedData/1E3,Lmittelvec/1E3,name)
    
    
    for i in range(nR1-1):
        rvecmL[i]=rvecL[i+1]**2-rvecL[i]**2
        rvecmR[i]=rvecR[i+1]**2-rvecR[i]**2
        LD1m[i]=(LD1[i]+LD1[i+1])/2
        RD1m[i]=(RD1[i]+RD1[i+1])/2
    
    VolL=np.sum(np.pi/2*LD1m*rvecmL)
    VolR=np.sum(np.pi/2*RD1m*rvecmR)
    Vol=VolL+VolR
    Volumen=Vol*(10**-7)**3
    rho=m/Volumen

    return rho,D,LD1,RD1,nR1,rvecR

def NormalizeOnSameRadius(Dav,D,rvecm,rvec,LD,RD):
    n=len(rvec)
    Li=[]
    Ri=[]
    for kk in range(n):
        rvei=(Dav/D[kk]+1E-4)*rvec[kk]
        LD_fun=interp1d(rvei,LD[kk])
        RD_fun=interp1d(rvei,RD[kk])
        Li.append(LD_fun(rvecm))
        Ri.append(RD_fun(rvecm))
        #string=["Li","Ri"]
        #DatatoCSV=["Li","Ri"]
        #WriteResultsCSV(string,DatatoCSV,name)
    return Li, Ri
   
    
    
def ProfilometerPlot(Li,Ri,turns,nPhi,nAV,rvecm,cmap,dens,name):
    Phi=np.pi
    n=len(Li)
    Qi=[]
    for kk in range(n):        
        turn=turns[kk]
        phivecL=np.linspace(Phi,0,nPhi)
        phivecR=np.linspace(0,Phi,nPhi)
        QL=np.zeros((nPhi,nAV))
        QR=np.zeros((nPhi,nAV))
        LD1=Li[kk]
        RD1=Ri[kk]
        for i in range(nAV):
            for j in range(nPhi):
                QL[j,i]=LD1[i]+phivecL[j]*(RD1[i]-LD1[i])/Phi
                QR[j,i]=RD1[i]+phivecR[j]*(LD1[i]-RD1[i])/Phi
        Q=np.zeros((2*nPhi,nAV))
        for j in range(nPhi):
            Q[j+nPhi,:]=QL[nPhi-j-1,:]
            Q[j,:]=QR[j,:]
        for i in range(nAV):
            Q[:,i]=np.roll(Q[:,i],int(turn/(np.pi)*nPhi))
        Qi.append(Q)
        #fig2, ax2= plt.subplots()

    Qave=np.mean(np.asarray(Qi),axis=0)
    phivec=np.hstack((-phivecL,phivecR))
    RR,PPHI=np.meshgrid(rvecm,phivec)
    X, Y=RR*np.cos(PPHI), RR*np.sin(PPHI)
    scale=1E3
    Qpl=Qave/scale
    Xpl=X/scale**2
    Ypl=Y/scale**2
    #fig1, ax1= plt.subplots()
    
    fig1 = plt.figure(dpi=96, figsize=(4604/(8*96),3602/(8*96)))
    ax1 = fig1.add_subplot(111, projection="3d")
    fig1.subplots_adjust(bottom = 0.6)
    norm = plt.Normalize(np.min(Qpl), np.max(Qpl))
    Xbackground,Ybackground=np.meshgrid(np.hstack((-rvecm[::-1],rvecm)),np.hstack((-rvecm[::-1],rvecm)))
    Qbackground=np.zeros((2*nAV,2*nAV))
    ax1.contourf(Xbackground/scale**2,Ybackground/scale**2,Qbackground/scale,dens,norm=norm,cmap=cmap)
    ax1.contourf(Xpl,Ypl,Qpl,dens,norm=norm,cmap=cmap)
    #ax1.plot_surface
    
    ax1.set_xlabel('X-Axis[mm]')
    ax1.set_ylabel('Y-Axis[mm]')
    
    m = cm.ScalarMappable(cmap=cmap, norm=norm)
    m.set_array(Qpl)
    cbar=plt.colorbar(m)
    cbar.set_label('[$\mu m$]', rotation=0)
    plt.show()
    cwd=os.getcwd()
    #fig1.savefig(join(cwd,"Plots",name+".jpeg"))
    return
 
def WriteResultsCSV(columnname,data,name):
    column_names=[]
    d=[]
    import os
    for i in range(len(columnname)):
        column_names.append(columnname[i])
        d.append(data[i])
    

    cwd=os.getcwd()
    export_data = zip_longest(*d, fillvalue = '')
    filecsv=join(cwd,'Profilometer',name+'.csv')
    with open(filecsv, 'w+', encoding="ISO-8859-1", newline='') as myfile:
        wr = csv.writer(myfile)
        wr.writerow((column_names))
        wr.writerows(export_data)
    myfile.close()
    csv2xlsx(filecsv)
    return

def csv2xlsx(filecsv):
    df_new = pd.read_csv(filecsv)
    writer = pd.ExcelWriter(filecsv[:-4]+'.xlsx')
    df_new.to_excel(writer, index = False)
    writer.save()
    return

def ThicknessPlot(rvec,correctedData,Lmittel,name):
    dpi=96
    fs=(4604/(8*dpi),3602/(8*dpi))
    fig6, ax6 = plt.subplots(figsize=fs,dpi=dpi)
    ax6.plot(rvec,correctedData,'k-',LineWidth=0.5)
    ax6.plot(rvec,Lmittel,'b-',LineWidth=0.5)
    yax6=ax6.set_ylabel(r'L'+'\n'+ r'      [$\mu m$]')
    ax6.set_xlabel(r'x-Axis [mm]')
    ax6.set_xlim([np.min(rvec),np.max(rvec)])
    yax6.set_rotation(0)
    yax6.set_position([150/dpi,100/dpi,222/dpi,138.5/dpi])
    fig6.show()
    cwd=os.getcwd()
    name=os.path.basename(name)
    #fig6.savefig(join(cwd,"Plots",name+'Thickness.jpeg'))
    import PyPlorigin.PyOrigin as OriginVorlagePlot#added 03.03.2020 for origin
    xData,yData,xSim,ySim,xErr,yErr=MusterDictOrigin()
    xlab='L%(CRLF)[\g(mu)m]'
    ylab='x-Axis[mm]'
    x=[]
    y=[]
    xe=[]
    ye=[]
    xs=[]
    ys=[]
    xu="mum"
    yu="mm"
    xc=xlab
    yc=ylab
    x.append(np.asarray(rvec))
    y.append(np.asarray(correctedData))
    xs.append(np.asarray(rvec))
    ys.append(np.asarray(Lmittel))
    xe.append(np.asarray([0]))
    ye.append(np.asarray([0]))
    xData["values"],yData["values"],xSim["values"],ySim["values"],xErr["values"],yErr["values"]=x,y,xs,ys,xe,ye
    xData["unit"],yData["unit"],xSim["unit"],ySim["unit"],xErr["unit"],yErr["unit"]=xu,yu,xu,yu,xu,yu
    xData["comment"],yData["comment"],xSim["comment"],ySim["comment"],xErr["comment"],yErr["comment"]=xc,yc,xc,yc,xc,yc
    OriginVorlagePlot.Plot(xData,yData,xSim,ySim,xErr,yErr,ylab,xlab,name+'Thickness')
    return