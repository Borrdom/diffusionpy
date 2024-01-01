import casadi as cs
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib as mtl
colors=mtl.colors.cnames
from os.path import join
cwd=os.getcwd()
import pandas as pd
from scipy.optimize import least_squares
from scipy.interpolate import InterpolatedUnivariateSpline
import time
from PyCSAFT import get_pcpar
from PyCSAFT.list_manipulations import unnesting
#sns.set()
from PyPlorigin.PyOrigin import TernaryPlot,Plot,MusterDictOrigin
start_time = time.time()
from random import random
from PySorption.DiffusionsmodellPlots import Plotlist
def cs2np(SX):
    dummy=cs.SX.sym('dummy')
    fun=cs.Function('fun',[dummy],[SX])
    arr=np.reshape(np.array(fun(1)),(SX.shape[0]))
    return arr

def InsertExpression(ins,var,expression):
	expression_fun=cs.Function("fun",[var],[expression])
	return expression_fun(ins)


def PyCSAFToutput(sol,expression,des):
	expression_res=[]
	expression_fun=cs.Function("fun",[des],[expression])
	for i,val in enumerate(sol):
		expression_res.append(cs2np(expression_fun(val)))
	return np.asarray(expression_res)

class Mixture:
    def __init__(self,*args,wPolyASD=1,dikij={},dilij={},dikijb={}):
        self.comp=args
        self.ncomp=len(args)
        self.ReInit()
        self.Lit=0
        self.wwASD=0
        self.rho0Poly0=None
        self.rho0API0=1300 #API density is also very sensitive 1320 vs 1300
        self.dictkij=dikij
        self.dictlij=dilij
        self.dictkijb=dikijb
        self.nGplot=20
        self.nVLE=20 #20
        self.nLLE=20
        self.nASD=21
        self.nSLE=20
        self.SpinL=[]
        self.SpinR=[]
        self.rhoASD_fun=[]
        self.nSpin=100
        self.idx=1
        self.idxp=0
        self.idxa=2

    def ReInit(self):
        self.x0=np.zeros(self.ncomp)
        self.x0[0]=1
        self.XAi0=np.zeros(self.ncomp)
        self.XAi0[0]=0.45
        self.eta0=0.34
        self.Z0=1


    def PCSAFT(self,psys=0.03166*1E5*0.5,T=298.15,na=""):
        #from PCSAFTFunctions import comp,InsertExpression,PyCSAFToutput
        k=1.380649*1E-23
        self.k=k
        Nav=6.022E23
        self.Nav=Nav
        Lambda=(10**10)**3/Nav
        npolyI=7
        npolycoef=3
        a0=cs.SX.ones(npolycoef,npolyI)
        b0=cs.SX.ones(npolycoef,npolyI)
        a0[:,0]=[0.9105631445,-0.3084016918,-0.0906148351]
        a0[:,1]=[0.6361281449,0.1860531159,0.4527842806]
        a0[:,2]=[2.6861347891,-2.5030047259,0.5962700728]
        a0[:,3]=[-26.547362491,21.419793629,-1.7241829131]
        a0[:,4]=[97.759208784,-65.255885330,-4.1302112531]
        a0[:,5]=[-159.59154087,83.318680481,13.776631870]
        a0[:,6]=[91.297774084,-33.746922930,-8.6728470368]

        b0[:,0]=[0.7240946941,-0.5755498075,0.0976883116]
        b0[:,1]=[2.2382791861,0.6995095521,-0.2557574982]
        b0[:,2]=[-4.0025849485,3.8925673390,-9.1558561530]
        b0[:,3]=[-21.003576815,-17.215471648,20.642075974]
        b0[:,4]=[26.855641363,192.67226447,-38.804430052]
        b0[:,5]=[206.55133841,-161.82646165,93.626774077]
        b0[:,6]=[-355.60235612,-165.20769346,-29.666905585]
        self.PyCSAFToutput=PyCSAFToutput

        Mi,mi,epsAiBi,kapi,ui,sigi,N=list(map(cs.SX,zip(*[(val["Mi"],val["mi"],val["epsAiBi"],val["kapi"],val["ui"],val["sigi"],val["N"]) for val in self.comp])))
        scomp=[val["name"] for val in self.comp]
        ncomp=len(scomp)
        #ncomp,Mi,mi,ui,epsAiBi,kapi,sigi,N,scomp=comp(*self.comp)
        self.ncomp=ncomp
        self.scomp=scomp
        ntype=2
        dictkij=self.dictkij
        dictlij=self.dictlij
        dictkijb=self.dictkijb
        gij=cs.SX.ones(ncomp,ncomp)
        deltAiBj=cs.SX.ones(ncomp,ncomp)
        epsAiBj=cs.SX.ones(ncomp,ncomp)
        kapij=cs.SX.ones(ncomp,ncomp)
        sigij=cs.SX.ones(ncomp,ncomp)
        dij=cs.SX.ones(ncomp,ncomp)
        uij=cs.SX.ones(ncomp,ncomp)

        #Initializekij
        kij=cs.SX.zeros(ncomp,ncomp)
        lij=cs.SX.zeros(ncomp,ncomp)
        for i in range(ncomp):
            co1=scomp[i]
            for j in range(ncomp):
                if i!=j:
                    try:
                        kij[i,j]=dictkij[co1+scomp[j]] if co1+scomp[j] in dictkij.keys() else dictkij[scomp[j]+co1]
                    except:
                        kij[i,j]=0

                    try:
                        lij[i,j]=dictlij[co1+scomp[j]] if co1+scomp[j] in dictkij.keys() else dictlij[scomp[j]+co1]
                    except:

                        lij[i,j]=0

        # Define decision variables
        self.xi=cs.SX.sym('xi'+na,ncomp)
        xi=self.xi
        self.eta=cs.SX.sym('eta'+na)
        eta=self.eta
        self.XAi=cs.SX.sym('XAi'+na,ncomp)
        XAi=self.XAi
        self.rho=cs.SX.sym('rho'+na)
        rho=self.rho
        self.Zsym=cs.SX.sym("Zsym"+na)
        Zsym=self.Zsym
        self.wi=cs.SX.sym("wi"+na,ncomp)
        wi=self.wi
        self.rhoi=cs.SX.sym('rhoi'+na,ncomp)
        rhoi=self.rhoi

        wieq=xi*Mi/(cs.sum1(xi*Mi))
        #wieq=xi*Mi/mi/(cs.sum1(xi*Mi/mi))
        xieq=wi/Mi/(cs.sum1(wi/Mi))
        #xieq=(wi/Mi*mi)/(cs.sum1(wi/Mi*mi))
        varnlp=cs.vertcat(xi,eta,XAi,Zsym)


        # Hard Chain Contribution
        mibar=cs.sum1(xi*mi)
        di=sigi*(1-0.12*cs.exp(-3*ui/T))
        rhoeq=6/cs.pi*eta*(cs.sum1(mi*xi*di**3))**-1 #mi*xi=phi*mbar und dann mbar*rho=rhophi
        #rhoeq=6/cs.pi*eta*(cs.sum1(xi*di**3))**-1
        veq=1/rhoeq
        vi=cs.gradient(veq,xi)
        etaeq=cs.pi/6*cs.sum1(mi*xi*di**3)
        #etaeq=cs.pi/6*rho*cs.sum1(mi*xi*di**3)
        #etaeq=cs.pi/6*rho*cs.sum1(xi*di**3)
        z3=cs.pi/6*rho*cs.sum1(mi*xi*di**3)
        #z3=cs.pi/6*rho*cs.sum1(xi*di**3)
        z2=cs.pi/6*rho*cs.sum1(mi*xi*di**2)
        #z2=cs.pi/6*rho*cs.sum1(xi*di**2)
        z1=cs.pi/6*rho*cs.sum1(mi*xi*di**1)
        #z1=cs.pi/6*rho*cs.sum1(xi*di**1)
        z0=cs.pi/6*rho*cs.sum1(mi*xi*di**0)
        #z0=cs.pi/6*rho*cs.sum1(xi*di**0)

        #Mixing rules
        def matmultadd(a):
            return (a*cs.SX.ones((a.shape[0],a.shape[0]))).T+a
        dij=matmultadd(di)
        sigij=matmultadd(sigi)/2
        uij=(ui@ui.T)**0.5*(1-kij)
        kapij=(kapi@kapi.T)**0.5*((sigi@sigi.T)**0.5/sigij)**3
        epsAiBj=matmultadd(epsAiBi)/2*(1-lij)
        gij=1/(1-z3)+(di@di.T/dij)*3*z2/(1-z3)**2+((di@di.T)/(dij))**2*2*z2**2/(1-z3)**3
        fhs=1/z0*(3*z1*z2/(1-z3)+z2**3/(z3*(1-z3)**2)+(z2**3/z3**2-z0)*cs.log(1-z3))
        fhc=mibar*fhs-cs.sum1(xi*(mi-1)*cs.log(cs.diag(gij)))
        #fhc=fhs-cs.sum1(xi*(1-1/mi)*cs.log(cs.diag(gij)))

        # Dispersion Contribution
        a=a0[0,:]+(mibar-1)/mibar*a0[1,:]+(mibar-1)/mibar*(mibar-2)/mibar*a0[2,:]
        b=b0[0,:]+(mibar-1)/mibar*b0[1,:]+(mibar-1)/mibar*(mibar-2)/mibar*b0[2,:]


        I1,I2=cs.polyval(a[::-1],z3),cs.polyval(b[::-1],z3)
        #I1,I2=cs.SX.zeros(ncomp),cs.SX.zeros(ncomp)
        #for i in range(ncomp):
        #    a=a0[0,:]+(mi[i]-1)/mi[i]*a0[1,:]+(mi[i]-1)/mi[i]*(mi[i]-2)/mi[i]*a0[2,:] #The integral for Segment based is now for every segment and inside the sum
        #    b=b0[0,:]+(mi[i]-1)/mi[i]*b0[1,:]+(mi[i]-1)/mi[i]*(mi[i]-2)/mi[i]*b0[2,:]
        #    I1[i]=cs.polyval(a[::-1],z3)
        #    I2[i]=cs.polyval(b[::-1],z3)

        Zhc=rho*cs.gradient(fhc,rho)
        C1=(1+Zhc+rho*cs.gradient(Zhc,rho))**(-1) #does this influence the lwidth of the LLE?07.06.20#Hard to say look at the binary LLE first. RIT seem to find the pure water phase
        #sould be ok for segmental
        m2es3mat=(xi@xi.T)*(mi@mi.T)*(uij/T)*sigij**3
        #m2es3mat=(xi@xi.T)*(uij/T)*sigij**3
        m2e2s3mat=m2es3mat*(uij/T)
        m2es3,m2e2s3=cs.sum1(cs.sum2(m2es3mat)),cs.sum1(cs.sum2(m2e2s3mat))
        #m2es3,m2e2s3=cs.sum1(I1*cs.sum2(m2es3mat)),cs.sum1(I2*cs.sum2(m2e2s3mat))
        #m2es3,m2e2s3=cs.sum1(cs.sum2(I1*m2es3mat)),cs.sum1(cs.sum2(I2*m2e2s3mat))

        fdisp=-2*cs.pi*rho*I1*m2es3-cs.pi*rho*mibar*C1*I2*m2e2s3
        #fdisp=-2*cs.pi*rho*m2es3-cs.pi*rho*C1*m2e2s3*mibar
        #fdisp=-2*cs.pi*rho*m2es3-cs.pi*rho*C1*m2e2s3#*mibar

        # Association Contribution
        deltAiBj=gij*kapij*sigij**3*(cs.exp(epsAiBj/T)-1)
        var=cs.vertcat(xi,rho,XAi)
        varrho=cs.vertcat(xi,rhoeq,XAi)
        XAi_eq=(1+cs.sum1(rho*xi*XAi*N*deltAiBj.T).T)**-1
        XAi_eq=(1+cs.sum1(rho*xi/mi*XAi*N*deltAiBj.T).T)**-1
        XAi_eqx=InsertExpression(varrho,var,XAi_eq)
        tm4=cs.sum1((cs.log(XAi)-XAi+1)*N*ntype*xi)
        #tm4=cs.sum1((cs.log(XAi)-XAi+1)*N*ntype*xi/mi)
        tme4=cs.sum1((cs.log(XAi)-1/2*XAi+1/2)*N*ntype*xi)
        #tme4=cs.sum1((cs.log(XAi)-1/2*XAi+1/2)*N*ntype*xi/mi)
        tmp4=cs.sum2(cs.sum1((XAi@XAi.T)*deltAiBj*(N@N.T)*(xi@xi.T)*ntype))
        #tmp4=cs.sum2(cs.sum1((XAi@XAi.T)*deltAiBj*(N@N.T)*(xi@xi.T)*ntype*(mi@mi.T)**-1))
        q=tm4-rho/2*tmp4



        #fassoc=tme4
        #Sum up Contriibutions
        fres=fhc+fdisp+q
        #fres=fhc+fdisp
        #Compute chemical potential another way
        rhosum=cs.sum1(rhoi)
        xinet=rhoi/rho
        rhoieq=xi*rhoeq
        rhoieq2=xi*rho
        fresrhoi=InsertExpression(xinet,xi,fres)
        Fres=InsertExpression(rhosum,rho,fresrhoi*rho)
        Q=InsertExpression(xinet,xi,q*rho)
        Q=InsertExpression(rhosum,rho,Q)
        mureseq=cs.gradient(Fres,rhoi)
        mures=cs.fmax(InsertExpression(rhoieq,rhoi,mureseq),-1E327)
        mures2=InsertExpression(rhoieq2,rhoi,mureseq)
        muid=cs.fmax(cs.log(xi*rhoeq*k*T),-1E327)
        #muid=cs.fmax(cs.log(rhoeq*xi/mi*k*T),-1E327)

        muid2=cs.log(xi*rho*k*T)
        #muid2=cs.log(xi*rho/mi*k*T)
        #dfresdx=InsertExpression(varrho,var,cs.gradient(fres,xi))
        dfresdrho=InsertExpression(varrho,var,cs.gradient(fres,rho))
        #THFaktor=2*rhoi/rho*cs.jacobian(mureseq*rho,rhoi)-cs.sum2(rhoi/rho*cs.jacobian(mureseq,rhoi))
        THFaktor=2*xi*cs.jacobian(muid2+mures2,xi)-cs.sum2(xi*cs.jacobian(muid2+mures2,xi))
        self.Gamma=InsertExpression(varrho,var,THFaktor)
        #Calculate densities
        rhobar=rhoeq/Nav*(10**10)**3
        #rhobar=cs.sum1(rhoeq*xi/mi)/Nav*(10**10)**3
        self.rhobar=rhobar
        rhopol=rhobar*xi*Mi
        #rhopol=(rhoeq*xi/mi)/Nav*(10**10)**3*Mi
        self.rhopol=rhopol
        rhosys=rhobar*cs.sum1(Mi*xi)
        #rhosys=cs.sum1((rhoeq*xi/mi*Mi))/Nav*(10**10)**3
        self.rhosys=rhosys
        #Xi=wi/(1-wi)
        #rhopolX=rhopol/Xi
        rhoblend=rhosys*(1-wieq)#-rhopol
        self.rhoblend=rhoblend
        #Calculate pressure
        Z=1+rhoeq*dfresdrho#+Zassoc

        p=(Zsym*rhoeq*k*T)*(10**10)**3
        #p=(Zsym*rhobar*k*Nav*T)


        #Calculate chemical potential
        self.mures=mures
        self.muid=muid
        self.lnphi=mures-cs.log(Zsym)
        self.lnfi=self.lnphi+cs.log(xi)+cs.log(p)
        self.lnfi=muid+mures
        self.lnai=self.lnphi+cs.log(xi)
        #Later cutoff reduces Breakoffs
        muideq=cs.log(rhoi*k*T)
        #muideq=cs.log(rhoi/mi*k*T)
        self.muges=self.muid+self.mures# Changed it so that chemical potential is never infinite but a big ass number
        self.g=cs.sum1(xi*self.muges)

        CD=cs.SX.eye(self.ncomp)
        #varnlp=cs.vertcat(xi,eta,XAi,Zsym)
        #var=cs.vertcat(xi,rho,XAi)

        #lnlfidni=cornecker/xi+rho*cs.jacobian(mureseq,rhoi)
        lnlfidni=rho*cs.jacobian(muideq+mureseq,rhoi)
        lnlfidniX=2*rho*cs.jacobian(muideq+mureseq,rhoi)-cs.sum2(rho*cs.jacobian(muideq+mureseq,rhoi))
        lnlfidni2=InsertExpression(varrho,var,InsertExpression(rhoieq,rhoi,rho*cs.jacobian(muideq+mureseq-cs.log(rhoi/rho),rhoi)))

        lnlfidniX2=InsertExpression(varrho,var,InsertExpression(rhoieq,rhoi,lnlfidniX))
        #self.Gamma=lnlfidniX2
        #mugesGamma=InsertExpression(cs.vertcat(1-cs.sum1(xi[1:]),xi[1:],rho,XAi),var,muideq+mureseq)
        #for i in range(ncomp):
        #    for j in range(ncomp):
        #        CD[i,j]+=xi[i]*cs.gradient(mugesGamma[i]-cs.log(xi)[i],xi[j])
        #self.Gamma=CD+lnlfidni2*cs.vertcat([xi]*self.ncomp)
        qnn=lnlfidni
        qxx=cs.jacobian(cs.gradient(Q,XAi),XAi)
        qxn=cs.jacobian(cs.gradient(Q,XAi),rhoi)
        qnx=cs.jacobian(cs.gradient(Q,rhoi),XAi)
        ann=qnn-rho*qnx@cs.inv(qxx)@qxn

        #-1E13
        #Stability determinant
        L1=cs.fmax(ann,-1E13)#[1:,1:] #Isotherm points seem to be unstable# Without fid they are all unstable
        drhoidxi=cs.jacobian(rhoieq,xi)

        self.L1=InsertExpression(varrho,var,InsertExpression(rhoieq,rhoi,L1))
        #self.Gamma=self.L1

        detL1=cs.det(L1)#L1[0,0]*L1[1,1]-L1[0,1]*L1[1,0]#cs.det(L1) #Basically right
        self.detL1=InsertExpression(varrho,var,InsertExpression(rhoieq,rhoi,detL1))#[1:,1:]
        #Crirical determinant
        ddetL1drhoi=rhoeq*cs.gradient(detL1,rhoi)
        M1=cs.horzcat(ddetL1drhoi,L1[:,1:])
        detM1=cs.det(M1)
        self.detM1=InsertExpression(varrho,var,InsertExpression(rhoieq,rhoi,detM1))


        #
        self.p=p
        self.Z=Z
        self.xieq=xieq
        self.wieq=wieq
        self.Jp=(Zsym-psys/(rhoeq*(10**10)**3*k*T))*100 #Is it Zsym or Z which is better 03.06.20
        #self.Jp=(Zsym-psys/(rhobar*k*Nav*T))*100
        self.Jx=(XAi-XAi_eqx)*100 #Also Multiplied by 100
        self.JZ=(Z-Zsym)*100 #Removed the factor 100
        self.Jbal=(cs.sum1(xi)-1)*100#Factor of one hundred seems to reduce"other solutions"
        self.varnlp=varnlp
        self.varxieq=cs.vertcat(xieq,eta,XAi,Zsym)
        self.varwnlp=cs.vertcat(wi,eta,XAi,Zsym)
        self.fres=fres
        self.XAi_eqx=XAi_eqx
        self.varscan=cs.vertcat(eta,XAi,Zsym)
        self.InsertExpression=InsertExpression
        self.vi=vi*Nav/(10**10)**3/Mi
        self.vmol=1/rhobar*1E6
        # Für NET
        self.rho0Poly=cs.SX.sym("rho0Poly")

        self.rho0Blend=cs.SX.sym("rho0Blend")
        #rho0Poly=self.rho0Poly
         #wozu brauche ich Znet !
        #cpoly=rho0Poly/Mi[0]/Lambda
        #ci=rho0Poly/Mi/Lambda
        #idxvec=np.asarray(list(range(self.ncomp)))
        #muresnetrholist=[InsertExpression(rhoieq[idxvec[idxvec!=i]],rhoi[idxvec[idxvec!=i]],mureseq) for i in range(ncomp)]
        #muresnetrholist2=[InsertExpression(ci[i],rhoi[i],muresnetrholist[i]) for i in range(ncomp)]
        #muresnetrho=InsertExpression(cs.vertcat(cpoly,rhoieq[1:]),rhoi,mureseq)
        #muresnetrho=muresnetrholist2[0]
        #self.mugesnetrho=[self.muid+muresnetrholist2[i]  for i in range(ncomp)]
        self.Jnetrho=(self.rho0Poly-self.rhopol)

        self.Jnetrho=(self.rho0Blend-self.rhoblend)

       #Net-SAFT extention
    def NETPCSAFT(self,psys=0.03166*1E5*0.5,T=298.15):
        idx=self.idx
        idxp=self.idxp
        rhoPolysys=self.rhoPolysys(psys=psys,T=T)
        #rho0Poly=self.rho0Poly
        varnlp=self.varnlp
        #print(rhoPolysys)
        InsertExpression=self.InsertExpression
        #mugesnetrho=self.mugesnetrho[idxp]
        #mugesnet=InsertExpression(cs.vertcat(varnlp,rhoPolysys),cs.vertcat(varnlp,self.rho0Poly),mugesnetrho)
        mugesnet=self.muges
        #Jnetrho=self.Jnetrho[idxp]
        Jnetrho=self.Jnetrho[idx]
        #Jnet=InsertExpression(cs.vertcat(varnlp,rhoPolysys),cs.vertcat(varnlp,self.rho0Poly),Jnetrho)
        Jnet=InsertExpression(cs.vertcat(varnlp,rhoPolysys),cs.vertcat(varnlp,self.rho0Blend),Jnetrho)
        self.mugesnet=mugesnet
        self.Jnet=Jnet
        Znet=1
        self.JZnet=(Znet-self.Zsym)
        return

    def PhaseInfo(self,res,success,psys=0.03166E5*0.1,T=298.15,NET=0,gammacalc=0):
        xi=res[0:self.ncomp].T.full()[0]
        if gammacalc==1:
            if not any(xi==1):
                sol=np.asarray([self.Pure(psys=psys,liq=True,T=T,com=i)["lnphii"][i] for i in range(self.ncomp)])
            else:
                sol=np.zeros(self.ncomp)
        else:
            sol=np.zeros(self.ncomp)
        muges=self.muges if NET==0 else self.mugesnet
        wi=self.PyCSAFToutput([res],self.wieq,self.varnlp)[0]
        mui=self.PyCSAFToutput([res],muges,self.varnlp)[0]
        lnphii=self.PyCSAFToutput([res],self.lnphi,self.varnlp)[0]

        #lnphii2=(mui-self.mupure)+self.lnphipure if not any(xi==1) else lnphii
        #lnphii=(mui-self.mupure)+self.lnphipure if NET==1 else lnphii
        #p=self.PyCSAFToutput([res],self.p,self.varnlp)[0]
        detL1=self.PyCSAFToutput([res],self.detL1,self.varnlp)[0]
        detM1=self.PyCSAFToutput([res],self.detM1,self.varnlp)[0]
        g=self.PyCSAFToutput([res],self.g,self.varnlp)[0]
        lnfi=self.PyCSAFToutput([res],self.lnfi,self.varnlp)[0]
        rho=self.PyCSAFToutput([res],self.rhosys,self.varnlp)[0]
        eta=self.PyCSAFToutput([res],self.eta,self.varnlp)[0]
        XAi=self.PyCSAFToutput([res],self.XAi,self.varnlp)[0]
        vmol=self.PyCSAFToutput([res],self.vmol,self.varnlp)[0]
        muxi=mui-cs.log(xi).full().T[0]
        vi=self.PyCSAFToutput([res],self.vi,self.varnlp)[0]
        TH=self.PyCSAFToutput([res],cs.reshape(self.Gamma,(self.ncomp*self.ncomp,1)),self.varnlp).T
        TH=np.reshape(TH,self.ncomp*self.ncomp)
        psim=self.PyCSAFToutput([res],self.p,self.varnlp)[0]
        #print(TH)
        #print(TH.shape)
        #print(self.Gamma.shape)
        #print(xi.shape)
        pH2OLV=1 if not hasattr(self,"pH2OLV") else self.pH2OLV
        RH=psys/pH2OLV #if not any(xi==1) else [1]
        #if NET==1:
        #    print(self.PyCSAFToutput([res],self.Z,self.varnlp)[0])
        #print(self.PyCSAFToutput([res],self.THFaktordet,self.varnlp)[0])
        lngammai=lnphii-sol if not any(xi==1) else np.log(xi)

        #print(mui-(sol["mui"]+cs.log(xi)+cs.log(gammai))) if not any(xi==1) else 1
        result={"xi":xi,
                "wi":wi,
                "mui":mui,
                "lnphii":lnphii,
                "lngammai":lngammai,
                "muxi":muxi,
                "success":success,
                "p":psys,
                "RH":RH,
                "Stability":detL1,
                "g": g,
                "Crit":detM1,
                "lnfi": lnfi,
                "rho":rho,
                "XAi":XAi,
                "eta":eta,
                "THf":TH,
                "vi":vi,
                "psim":psim,
                "vmol":vmol}
        return result


    def rhoPolysys(self,psys=0.03166*1E5*0.5,T=298.15):
        idx=self.idx
        idxp=self.idxp
        wPolyASD=np.atleast_1d(self.wPolyASD)
        idxp=np.atleast_1d(idxp)
        idx=np.atleast_1d(idx)
        kswvec=np.asarray([val["ksw"] if val["ksw"]!=None else 0 for i,val in enumerate(self.comp)])
        rho0Poly0vec=np.asarray([val["rho0Poly0"] if val["rho0Poly0"]!=None else self.rhopures[i] for i,val in enumerate(self.comp)])
        #rhoF=((1-self.wPolyASD)/rhoAPIpure+self.wPolyASD/rho0Poly0)**-1
        idxvec=np.asarray(list(range(self.ncomp)))
        idxa=np.hstack((idxp,idx))
        vlast=(1-sum(wPolyASD))/rho0Poly0vec[np.delete(idxvec,idxa)] if np.delete(idxvec,idxa).tolist()!=[] else 0
        rhoF=(vlast+sum([val1/rho0Poly0vec[val2] for i,(val1,val2) in enumerate(zip(wPolyASD,idxp))]))**-1
        kswmix=sum([kswvec[val2]*val1*rhoF/rho0Poly0vec[val2] for i,(val1,val2) in enumerate(zip(wPolyASD,idxp))])
        nsw=2#self.comp[int(idxp)]["nsw"]
        rhoFPol=self.wPolyASD*rhoF#Mixing rule API

        p_fun=(1-kswmix*(psys/self.pH2OLV)**nsw)
        #Dilatation approach confirmed by Sarti DOI:10.1021/acs.iecr.9b05453
        print(rhoF)
        return rhoF*p_fun

    def VLE(self, psys=0.03166*1E5*0.5,T=298.15,NET=0):

        """Calculates the VLE of component 1 and 2 where 1 is assumed to have no vapor pressure"""
        #from PCSAFTFunctions import PyCSAFToutput
        idx=self.idx
        idxp=self.idxp

        self.PCSAFT(psys=psys,T=T)
        self.pH2OLV=self.PureVaporPressure(T=T) if not hasattr(self, 'pH2OLV') else self.pH2OLV
        self.NETPCSAFT(psys=psys,T=T) if NET==1 else None
        Purevap=self.Pure(psys=psys,liq=False,T=T,com=idx)
        #Pureliq=self.Pure(psys=psys,liq=True,T=T,com=idx)
        muv1=Purevap["mui"][idx]
        #phiv1=np.exp(Purevap["lnphii"][idx])

        #Jmu=(muv1-self.muges[1])/100
        mul1=self.muges[idx] if NET==0 else self.mugesnet[idx]
        Jmu=(muv1-mul1) #Changed from /100 to /10 at 02.06.20
        #Jmu=(cs.log(psys/self.pH2OLV)-(cs.log(self.xi[idx])+self.lnphi[idx]-philpure+cs.log(psys/self.p)))#)
        #Jmu=(cs.log(psys*phiv1)-self.lnfi[idx]) #ideale Gasphase
        #Jmu=(cs.log(self.xi[idx])+self.lnphi[idx])
        JZ=self.JZ if NET==0 else self.JZnet
        Jp=self.Jp if NET==0 else self.Jnet
        #Jp=self.Jp if NET==0 else

        wPolyASD=np.atleast_1d(self.wPolyASD)
        wieq=self.wieq
        #Changed from wieq[0] to wieq[2] must be accurate fopr binary system 04.06.20
        w1eq=(wPolyASD)*(1-wieq[idx])
        w2eq=wieq[idx]
        w3eq=1-cs.sum1(w1eq)-w2eq
        weq=cs.vertcat(w1eq,w2eq,w3eq)
        idxvec=np.asarray(list(range(self.ncomp)))
        idxa=np.asarray([idxp,idx])
        weq=cs.SX.zeros(self.ncomp)
        weq[idxp]=w1eq

        weq[idx]=w2eq
        if np.delete(idxvec,idxa).tolist()!=[]:
            weq[np.delete(idxvec,idxa)]=w3eq
        x3eq=self.InsertExpression(weq,self.wi,self.xieq)
        #significant better for pvpritonavirwater
        self.JASD=(self.xi-x3eq)*100  #Factor of 100 heare increase decreases "other solutions"
        JASD=self.JASD[idxp] if self.ncomp>2 else []
        #JASD=self.wPolyASD*(1-self.wieq[idx])
        J=cs.vertcat(Jp,self.Jx,Jmu,self.Jbal,JZ,JASD)
        res,success=self.SolveNLP(self.varnlp,True,J)
        result=self.PhaseInfo(res,success,psys,T=T,gammacalc=1)
        print(self.x0)
        return result

    def WriteToExcel(self,listofdicts,nameofsheets,Datatype="LLE"):
        def DictToSheet(info):
            dfinfo=pd.DataFrame.from_dict(info)
            dfinfoun=unnesting(dfinfo,dfinfo.columns.tolist(),axis=0)
            return dfinfoun
        listofsheets=[DictToSheet(valdict) for i,valdict in enumerate(listofdicts)]
        components="_".join([val["name"] for i,val in enumerate(self.comp)])
        filename=Datatype+"_"+components
        pathname=join(cwd,"PyCSAFTSheets",filename)
        writer = pd.ExcelWriter(pathname+'.xlsx', engine='xlsxwriter')
        [listofsheets[i].to_excel(writer,sheet_name=valname) for i,valname in enumerate(nameofsheets)]
        [listofsheets[i].to_csv(pathname+"_"+valname+'.csv') for i,valname in enumerate(nameofsheets)]
        #listofsheets[0].to_excel(writer,sheet_name="example")
        listofcomponent=DictToSheet(self.comp)
        listofcomponent.to_excel(writer,sheet_name="PureParameters")
        dictkij=self.dictkij
        dictkij["name"]="kij"
        dictlij=self.dictlij
        dictlij["name"]="lij"
        listofkij=DictToSheet([dictkij,dictlij])
        listofkij.to_excel(writer,sheet_name="BinaryParamers")
        writer.save()



    def SLE(self, psys=1.013E5,wwASD=0,T=298.15):
        """Calculates the SLE of component 1 2 and 3 where 1 and 2 are assumed to have no melting Point"""
        #from PCSAFTFunctions import PyCSAFToutput
        idx=self.idx
        idxa=self.idxa
        # idxp=self.idxp
        self.PCSAFT(psys=psys,T=T)
        self.pH2OLV=self.PureVaporPressure(T=T) if not hasattr(self, 'pH2OLV') else self.pH2OLV
        # xi=self.xi
        # wi=self.wi
        # xieq=self.xieq
        # wieq=self.wieq
        # IE=self.InsertExpression
        # w2=wwASD
        # w1=(1-wwASD)/2
        # w3=1-w1-w2
        # w0i=cs.vertcat(w1,w2,w3)
        # wi0=np.zeros(self.ncomp)
        # #wseed=np.vstack((w1seedvec,w2seedvec,w3seedvec))
        # wi0=self.idx_fun(wi0,w0i,idx)
        # x0i=IE(w0i,wi,xieq)

        # self.PCSAFT(psys=psys,T=T)
        # self.JSLE=self.xi[idx]-x0i[idx]
        self.JSLE=self.wieq[idx]-wwASD if self.ncomp>2 else []
        #muv3=self.Pure(psys=psys,liq=True,com=2,T=T)["mui"][2] #Chemical potential of pure API
        lnf30=self.Pure(psys=psys,liq=True,com=idxa,T=T)["lnfi"][idxa]
        #lnphi3pure=self.lnphipure[2]

        deltaHSL=self.comp[idxa]["deltaHSL"]
        Tm=self.comp[idxa]["Tm"]
        cpSL=self.comp[idxa]["cpSL"]

        T=T
        R=self.k*self.Nav
        lnfi3=self.lnfi[idxa]
        SLE=-deltaHSL/(R*T)*(1-T/Tm)+cpSL/R*(Tm/T-1-cs.log(Tm/T))
        Jmu=(SLE-(lnfi3-lnf30))/10
        #Jmu=(SLE-(self.lnfi[2]))/10
        J=cs.vertcat(self.Jp,self.Jx,Jmu,self.Jbal,self.JZ,self.JSLE)
        res,success=self.SolveNLP(self.varnlp,True,J)
        result=self.PhaseInfo(res,success,psys,T=T,gammacalc=1)

        return result
    def idx_fun(self,wi0,wseed,idx):
        idxvec=np.asarray(list(range(self.ncomp)))
        if not (idx>self.ncomp-1):
            wi0[idx]=wseed[0]
        for i,vali in enumerate(idxvec[idxvec!=idx]):
            wi0[vali]=wseed[i+1]
        return wi0
    def Spinodal(self,psys=1.013E5,T=298.18,ws=0):
        idxp=self.idxp
        w2seedvec=np.linspace(1E-4,(1-ws)-1E-4,self.nSpin)
        left=np.geomspace(1E-4,(1-ws)/2,self.nSpin//2+1)
        right=(1-ws)-left[::-1][2::]
        w2seedvec=np.hstack((left,right))
        w1seedvec=np.ones_like(w2seedvec)*ws
        w3seedvec=np.fmin(np.fmax(1-w1seedvec-w2seedvec,0),1)
        wi0=np.zeros((self.ncomp,self.nSpin))
        wseed=np.vstack((w1seedvec,w2seedvec,w3seedvec))
        wi0=self.idx_fun(wi0,wseed,idxp) if self.ncomp>2 else np.vstack((w2seedvec,w3seedvec))
        Stability=np.asarray([self.ConcAsInput(psys=psys,liq=True,wres=wi0[:,i],T=T)["Stability"] for i,val in enumerate(w1seedvec)]).T[0]
        roots=InterpolatedUnivariateSpline(w2seedvec, Stability).roots()
        wnull=w2seedvec[np.argmin(Stability)]
        leftside=roots[wnull-roots>0]
        rightside=roots[wnull-roots<0]
        UB2=leftside[np.argmin(abs(wnull-leftside))] if np.any(leftside) else 0
        LB2=rightside[np.argmin(abs(wnull-rightside))] if np.any(rightside) else 1
        UB3=(1-ws)-UB2
        LB3=(1-ws)-LB2
        wnull3=(1-ws)-wnull
        WLEFT=np.zeros(self.ncomp)
        WRIGHT=np.zeros(self.ncomp)
        WNULL=np.zeros(self.ncomp)
        WLEFT=self.idx_fun(WLEFT,np.hstack((ws,UB2,UB3)),idxp)
        WRIGHT=self.idx_fun(WRIGHT,np.hstack((ws,LB2,LB3)),idxp)
        WNULL=self.idx_fun(WNULL,np.hstack((ws,wnull,wnull3)),idxp)

        return WLEFT,WRIGHT,WNULL,Stability
    def LLE(self, psys=1.013E5,ws=0,T=298.15):
        """Calculates the SLE of component 1 2 and 3 where 1 and 2 are assumed to have no melting Point"""
        idxp=self.idxp
        #self.PCSAFT(psys=psys,T=T,na="L1")
        self.pH2OLV=self.PureVaporPressure(T=T) if not hasattr(self, 'pH2OLV') else self.pH2OLV

        WLEFT,WRIGHT,WNULL,Stability=self.Spinodal(psys=psys,T=T,ws=ws)# if False else np.asarray([0,1,0]),np.asarray([0,0,1]),wi0 # Too many function calls

        self.SpinL.append(WLEFT)
        self.SpinR.append(WRIGHT)
        self.PCSAFT(psys=psys,T=T,na="L1")
        varnlpI=cs.vertcat(self.xi,self.eta,self.XAi,self.Zsym)
        JpI=self.Jp
        JxI=self.Jx
        JbalI=self.Jbal
        JZI=self.JZ
        JI=cs.vertcat(JpI,JxI,JbalI,JZI)
        wiI=self.wieq
        xiI=self.xi
        mugesI=self.muges
        lnphiI=self.lnphi

        self.PCSAFT(psys=psys,T=T,na="L2")
        #mugesII=self.muges/1000 #mugesnet aber mit RH=1
        mugesII=self.muges
        JpII=self.Jp
        JxII=self.Jx
        JbalII=self.Jbal
        JZII=self.JZ
        varnlpII=cs.vertcat(self.xi,self.eta,self.XAi,self.Zsym)
        xiII=self.xi
        wiII=self.wieq
        lnphiII=self.lnphi

        Jmu=(mugesI-mugesII)
        for i in range(self.ncomp):
            Jmu[i]=(cs.conditional(xiI[i]<1E-12,[mugesI[i]],mugesII[i])-cs.conditional(xiII[i]<1E-12,[mugesII[i]],mugesI[i]))

        c2gibbs= lambda a,b,c=0: (1/2*(a+2*b)/(a+b+c),3**0.5/2*(a)/(a+b+c))
        dist=lambda x1,x2,y1,y2: ((x1-x2)**2+(y1-y2)**2)**0.5
        w2seed=(1-ws)/2 #is assumed to be in the middle# Make it critical point or wnull
        w3seed=1-ws-w2seed #Try to restrict phase fraction to 0 1 seems not to reduce nan
        wi0=np.zeros(self.ncomp)

        wi0=self.idx_fun(wi0,np.hstack((ws,w2seed,w3seed)),idxp)# if self.ncomp>2 else np.vstack((w2seed,w3seed))
        wI=cs.SX.ones(self.ncomp)
        wII=cs.SX.ones(self.ncomp)
        wI=self.idx_fun(wI,wiI,idxp) if self.ncomp>2 else wI
        wII=self.idx_fun(wII,wiII,idxp) if self.ncomp>2 else wII

        #xF,yF=c2gibbs(ws,w2seed,w3seed) #give component with index i
        xF,yF=c2gibbs(*tuple(wi0.tolist()))
        xL1,yL1=c2gibbs(*tuple([wI[i] for i in range(self.ncomp)]))
        #xL1,yL1=c2gibbs(wiI[0],wiI[1],wiI[2])
        xL2,yL2=c2gibbs(*tuple([wII[i] for i in range(self.ncomp)]))
        #xL2,yL2=c2gibbs(wiII[0],wiII[1],wiII[2])
        d1=dist(xF,xL1,yF,yL1)
        dF=dist(xL1,xL2,yL1,yL2)
        xhi1=d1/dF #Maybe it was the limit on feed fraction. No API is more stable but less accurate

        JLLE=xhi1*wI[0]+(1-xhi1)*wII[0]-ws if self.ncomp>2 else []#if ws==1E-4 else (wiII[idx]-ws)
        J=cs.vertcat(JI,JpII,JxII,Jmu,JbalII,JZII,JLLE) #Assoziation von API sorgt für Probleme


        wi=self.wi
        IE=self.InsertExpression
        xieq=self.xieq
        self.x0iL1=np.fmax(IE(WLEFT,wi,xieq).full().T[0],1E-310)
        self.x0iL2=np.fmax(IE(WRIGHT,wi,xieq).full().T[0],1E-310)
        self.x0ibo=np.fmax(IE(WNULL,wi,xieq).full().T[0],1E-310)

        g=cs.vertcat(-self.detL1/1E8+1)# if ws>1E-1 else gL1#this helps quite  alot#With this binary is also safe
        f=1
        noLLE=np.ones(self.ncomp)
        base={"xi":noLLE,
                "wi":noLLE,
                "mui":noLLE,
                "lnphii":noLLE,
                "lngammai":noLLE,
                "muxi":noLLE,
                "success":0,
                "p":np.asarray([0]),
                "RH":np.asarray([0]),
                "Stability":np.asarray([0]),
                "g": np.asarray([0]),
                "Crit":np.asarray([0]),
                "lnfi": noLLE,
                "rho":np.asarray([1000]),
                "XAi":noLLE}
        if any(Stability<0):
            res,success=self.SolveNLP(cs.vertcat(varnlpI,varnlpII),True,J,f=f,g=g) # #Problems when ws gets too small
            result1=self.PhaseInfo(cs.vertcat(res[0:varnlpI.shape[0]]),success,psys=psys,T=T)
            result2=self.PhaseInfo(res[varnlpI.shape[0]:varnlpI.shape[0]+varnlpII.shape[0]],success,psys=psys,T=T)
        else:

            result1=base
            result2=base

        #res,success=self.SolveNLP(self.varnlp,True,J)
        #result=self.PhaseInfo(res,success,psys,T=T)
        return result1,result2

    def PureVaporPressure(self,T=298.15):
        """Calculates vapor pressure pf component 1"""
        idx=self.idx
        def ObjFun(p0):
            self.PCSAFT(psys=p0,T=T)
            muv1=self.Pure(psys=p0,liq=False,T=T,com=idx)["mui"][idx]
            mul1=self.Pure(psys=p0,liq=True,T=T,com=idx)["mui"][idx]
            return (muv1-mul1)*10
        res=least_squares(ObjFun,0.03*1E5)
        #experimental 11.01.21
        self.pH2OLV=res.x
        return res.x



    def boundsandinits(self,x,liq=True):
        """Sets up standard boundaries and initials for SolveNLP"""

        lbx=[]
        ubx=[]
        x0=[]

        xnames=[x[i]for i in range(x.shape[0])]
        for i in range(x.shape[0]):
            if "xi_" in str(x[i]):# or "wi_0" in str(x[i]):
                idxc=int(str(x[i])[-1])
                x0c=self.x0[idxc]
                #lbx.append(self.x1lb)
                lbx.append(0)
                #ubx.append(self.x1ub)
                ubx.append(1)
                x0.append(x0c)

            # elif "xi_1" in str(x[i]) or "wi_1" in str(x[i]):
            #     x02=self.x02
            #     lbx.append(self.x2lb)
            #     ubx.append(self.x2ub)
            #     x0.append(x02)

            # elif "xi_2" in str(x[i]) or "wi_2" in str(x[i]):
            #     x03=self.x03
            #     lbx.append(self.x3lb)
            #     ubx.append(self.x3ub)
            #     x0.append(x03)

            elif "XAi_"in str(x[i]):
                idxc=int(str(x[i])[-1])
                XAi0c=self.XAi0[idxc]
                lbx.append(0)
                ubx.append(1)
                x0.append(XAi0c)
            # elif "XAi_1"in str(x[i]):
            #     XAi02=self.XAi02
            #     lbx.append(0)
            #     ubx.append(1-0)
            #     x0.append(XAi02)
            # elif "XAi_2"in str(x[i]):
            #     XAi03=self.XAi03
            #     lbx.append(0)
            #     ubx.append(1-0)
            #     x0.append(XAi03)
            elif "eta" in str(x[i]):
                if liq:
                    eta0=self.eta0 #for densities higher than 600 K etalb muist be smaller
                    lbx.append(0.1) #Changed from 0.33 to 0.1 03.06.20
                    ubx.append(1) #Changed from 0.55 to 0.43 19.04.20 0.47 is not problem for LLE
                    x0.append(eta0)
                else:
                    lbx.append(0)
                    ubx.append(0.1)
                    x0.append(10**(-10))
            elif "xiL2_0" in str(x[i]) or "wi_0" in str(x[i]):
                x01=self.x0iL2[0]
                lbx.append(0)
                ubx.append(1)
                x0.append(x01)

            elif "xiL2_1" in str(x[i]) or "wi_1" in str(x[i]):
                x02=self.x0iL2[1]
                lbx.append(self.x0ibo[1])
                ubx.append(1)
                x0.append(x02)

            elif "xiL2_2" in str(x[i]) or "wi_2" in str(x[i]):
                x03=self.x0iL2[2]
                lbx.append(0)
                ubx.append(self.x0ibo[2])
                x0.append(x03)
            elif "xiL1_0" in str(x[i]) or "wi_0" in str(x[i]):
                x01=self.x0iL1[0]
                lbx.append(0)
                ubx.append(1)
                x0.append(x01)

            elif "xiL1_1" in str(x[i]) or "wi_1" in str(x[i]):
                x02=self.x0iL1[1]
                lbx.append(0)
                ubx.append(self.x0ibo[1])
                x0.append(x02)

            elif "xiL1_2" in str(x[i]) or "wi_2" in str(x[i]):
                x03=self.x0iL1[2]
                lbx.append(self.x0ibo[2])
                ubx.append(1-0)
                x0.append(x03)

            elif "XAiL2_0" in str(x[i]) or "XAiL1_0"in str(x[i]):
                XAi01=random()
                lbx.append(0)
                ubx.append(1-0)
                x0.append(XAi01)
            elif "XAiL2_1" in str(x[i]) or "XAiL1_1" in str(x[i]):
                XAi02=random()
                lbx.append(0)
                ubx.append(1-0)
                x0.append(XAi02)
            elif "XAiL2_2" in str(x[i]) or "XAiL1_2"in str(x[i]):
                XAi03=random()
                lbx.append(0)
                ubx.append(1-0)
                x0.append(XAi03)
            elif "K" in str(x[i]):
                lbx.append(0)
                ubx.append(1E13)
                x0.append(1)
            elif "etaL2" in str(x[i]) or "etaL1" in str(x[i]):
                if liq:
                    eta0=self.eta0
                    lbx.append(0.33) #Changed from 0.33 to 0.1 03.06.20 # changed from 0.33 to 0.3 03.07.20
                    ubx.append(0.55) #Changed from 0.55 to 0.43 19.04.20 0.47 is not problem for LLE
                    x0.append(eta0)
                else:
                    lbx.append(0)
                    ubx.append(0.1)
                    x0.append(10**(-10))

            else:
                Z0=self.Z0
                lbx.append(0)
                ubx.append(10000) #reduction from inf to 100 lead to less spikes #Changed from 1000 to 10000 03.07.20 better binary LLE
                x0.append(Z0)
        return lbx, ubx,x0


    def SolveNLP(self,x,liq,J,f=1,g=cs.DM([])):
        """Solves NLP. x=decision variables,liq=bool weather liquid or gas, J=objective"""
        opts = {}
        opts["nlpsol"] = "ipopt"
        #opts['iteration_callback'] = mycallback
        #opts["verbose"]=True
        opts["nlpsol_options"] = {"ipopt.hessian_approximation":"limited-memory", #exact hessians seem to work now?
            "ipopt.acceptable_tol" : 1E-8, #1E-6 for LLE
            "ipopt.tol" : 1E-8,
            "ipopt.print_level" : 0,
            "ipopt.max_iter": 1000, #1000 in LLE problematic
            "ipopt.bound_relax_factor": 0,
            "ipopt.warm_start_init_point":       "yes",
            "ipopt.warm_start_bound_push":       1e-8,
            "ipopt.warm_start_bound_frac":       1e-8,
            "ipopt.warm_start_slack_bound_frac": 1e-8,
            "ipopt.warm_start_slack_bound_push": 1e-8,
            "ipopt.warm_start_mult_bound_push":  1e-8,
            "ipopt.linear_solver":  "mumps",
            "ipopt.fast_step_computation": "yes",
            "print_time":0}
            #"ipopt.recalc_y": "no",
            #"ipopt.recalc_y_feas_tol": 1E-6,
            #"ipopt.acceptable_dual_inf_tol": 1E-6,
            #"ipopt.hessian_constant":"yes",
            #"ipopt.jac_d_constant":"yes"}#hessian constant good for LLE
            #"ipopt.nlp_scaling_method":"no"}#,
            #"iteration_callback": mycallback}
            #"ipopt.constr_viol_tol":1E-10,
            #"ipopt.acceptable_constr_viol_tol":1E-10,
            #"ipopt.sb":"yes",
            #"print_time":0} #supresses ipopt output completly
            #"verbose": True,
            #"gather_stats": True,
            #"calc_g":True}
            #"ipopt.aacceptable_constr_viol_tol": 1E-6}
            #"ipopt.derivative_test": "first-order"}
        optss=opts["nlpsol_options"]
        self.optss=optss
        lbx,ubx,x0=self.boundsandinits(x,liq)
        xnames=[str(x[i]) for i in range(x.shape[0])]
        breakc=2
        if g.shape[0]>0:
            print("LLE!")
            breakc=20 #repeating is helping
        lbg=[0]*J.shape[0]
        ubg=[0]*J.shape[0]
        J=cs.vertcat(J,g) if "xi" in xnames[0]  else J #21.12.20 i do not remember what this was useful for# i think it was for the pure case where i insert all xi
        lbg=lbg+g.shape[0]*[-cs.inf] if "xi" in xnames[0] else lbg
        ubg=ubg+g.shape[0]*[0] if "xi" in xnames[0] else ubg
        J_fun=cs.Function("J_fun",[x],[J])
        g_fun=cs.Function("g_fun",[x],[g])
        J_fun(x0)
        g_fun(x0)
        nlp = {'x':x, 'f':f, 'g':J} # 03.06.20 Instable Points are rejected
        S = cs.nlpsol('S', 'ipopt', nlp,optss)

        res = S(x0=x0,lbx=lbx, ubx=ubx,lbg=lbg,ubg=ubg)
        retstat=S.stats()["return_status"]
        x_opt = res['x']
        count=0

        print(S.stats()["return_status"])
        while retstat!="Solve_Succeeded" : # 29.05.20 found out about retstat

            count+=1 if retstat!="Maximum_Iterations_Exceeded" else 1

            ran=10**((-200)*random()) #Chooses a random concentration in log fashion seems to be better for polymers
            ran2=10**((-200)*random()) #-130 was ok
            self.x0[0]=ran
            self.XAi0[0]=ran2
            for i in range(1,self.ncomp):
                self.x0[i]=(1-ran)/(self.ncomp-1)
                self.XAi0[i]=(1-ran2)/(self.ncomp-1)
            #self.x02=(1-ran)/2
            #self.x03=(1-ran)/2 #Starting at high water concentrations seems to slow down conversio
            self.eta0=(0.55-0.33)*random()+0.33 #Only free molar ratios 280520
            self.Z0=random()*1+0.01
            #self.XAi01=random() #23.06.20
            #self.XAi02=random()
            #self.XAi03=random()
            #self.XAi01=(1-ran2)/2 #23.06.20
            #self.XAi02=ran2
            #self.XAi03=(1-ran2)/2
            #print(self.x1lb)
            #print(self.x1ub)
            #self.ReInit() #Resets boundaries in LLE Routine
            lbx,ubx,x0=self.boundsandinits(x,liq)
            nlp = {'x':x, 'f':f, 'g':J}
            S = cs.nlpsol('S', 'ipopt', nlp,optss)
            res= S(x0=x0,lbx=lbx, ubx=ubx,lbg=lbg,ubg=ubg)
            retstat=S.stats()["return_status"]
            print(S.stats()["return_status"])
            # 0 leads to nan in JX. 1E-323 is the boundary# Still water Jmu2 leads  -inf possibility
            # Leads to infinity of chemical potential X of zero
            #Only one condition seems not to be fullfilled
            #J_fun(x_opt)
#DM([0, 0, 5.55112e-015, 1.11022e-014, 0, -4.38859e-012, -3.38813e-019, 0, -2.08167e-015, -5.55112e-015, 84.9586, 7.10543e-015, -2.84217e-013, 0, -1.499e-012, -2.77556e-017, -9.69782e+006])
            #x_opt
#DM([2.62281e-162, 0.60048, 0.39952, 0.467509, 0.935666, 0.140399, 0.577768, 0.000332801, 0.000110031, 0.99989, 9.5881e-008, 0.450501, 0.804794, 0.0477105, 0.293464, 2.90661e-005, 343.16, 343.16, 343.16])
            x_opt = res['x']
            if count>breakc:
                print("BreakOUT!")
                print("BreakOUT!")
                print("BreakOUT!")
                print("BreakOUT!")
                print("BreakOUT!")
                print("BreakOUT!")
                break

        success=1 if S.stats()["return_status"]=="Solve_Succeeded" else 0
        #Does the initials effect LLE (23.06.20) not really
        if x.shape[0]>2:
            for i in range(x.shape[0]):
                if "xi_" in str(x[i]):
                    idxc=int(str(x[i])[-1])
                    self.x0[idxc]=x_opt[i].full()[0]
                elif "XAi_" in str(x[i]):
                    idxc=int(str(x[i])[-1])
                    self.XAi0[idxc]=x_opt[i].full()[0]
                elif "eta" in str(x[i]):
                    self.eta0=x_opt[i].full()[0]
                else:
                    self.Z0=x_opt[i].full()[0]
        if g.shape[0]>0:
            self.ReInit()
        return x_opt,success


    def Pure(self,psys=0.03166*1E5*0.1,liq=True,com=1,T=298.15):
        #self.PCSAFT(psys=psys,T=T) #changed after binarLLE change i dont know why
        #It seems to make problem with the LLE. Needed for lnphiipure
        #from PCSAFTFunctions import PyCSAFToutput
        xres=cs.SX.zeros(self.ncomp)
        xres[com]=1
        xpure=cs.DM.zeros(self.ncomp)
        xpure[com]=1
        self.varpure=cs.vertcat(xres,self.eta,self.XAi,self.Zsym)
        J=cs.vertcat(self.Jp,self.Jx,self.JZ)
        Jscan=self.InsertExpression(self.varpure,self.varnlp,J)
        #res,success=self.SolveNLP(self.varscan,liq,Jscan)
        res,success=self.SolveNLP(self.varscan,liq,Jscan)

        sol=cs.vertcat(xpure,res)
        result=self.PhaseInfo(sol,success,psys=psys,T=T)
        if liq==True:
            self.lnphipure=result["lnphii"]
            self.mupure=result["mui"]
            self.rhow=result["rho"]
        else:
            self.lnphivap=result["lnphii"]

        return result#self.PyCSAFToutput([sol],self.muges,self.des)[0]

    def ConcAsInput(self,psys=0.03166*1E5*0.1,liq=True,wres=np.asarray([0,1,0]),T=298.15):
        self.PCSAFT(psys=psys,T=T)
        wi=self.wi
        xieq=self.xieq
        IE=self.InsertExpression
        xres=IE(wres,wi,xieq)
        self.varpure=cs.vertcat(xres,self.eta,self.XAi,self.Zsym)
        J=cs.vertcat(self.Jp,self.Jx,self.JZ)
        Jscan=self.InsertExpression(self.varpure,self.varnlp,J)
        res,success=self.SolveNLP(self.varscan,liq,Jscan)
        x_opt=res
        res=cs.vertcat(xres,x_opt)
        result=self.PhaseInfo(res,success,psys,T=T,gammacalc=1)
        return result


    def GPlot(self,psys=0.03166*1E5*1,T=298.15):
        wscan=np.linspace(0,1,self.nGplot)
        w2vec=np.asarray([])
        w1vec=np.asarray([])
        for i,val in enumerate(wscan):
            temp=wscan[wscan-1+val<0]
            w1vec=np.hstack((w1vec,np.ones_like(temp)*val))
            w2vec=np.hstack((w2vec,temp))
        w3vec=np.fmin(np.fmax(1-w1vec-w2vec,0),1)
        #w3vec[w3vec<1E-8]=0
        wvec=np.asarray([w1vec,w2vec,w3vec]) if self.ncomp>2 else np.asarray([w1vec,w2vec]) #wvec[X,Y] X is numper of components Y is number of points
        sol=[]
        for i,val in enumerate(wvec[1,:]):
            sol.append(self.ConcAsInput(psys=psys,liq=True,wres=wvec[:,i],T=T))
        return sol




    def rhoASD(self,psys=0.03166*1E5*0.1,liq=True,T=298.15,n=1000):
        idxp=self.idxp
        self.PCSAFT(psys=psys,T=T) #PCSAFT call within PCSAFT call leads to problems
        J=cs.vertcat(self.Jp,self.Jx,self.JZ)
        sol=[]
        idxvec=np.asarray(list(range(self.ncomp)))
        wres=np.zeros((self.ncomp,n))
        w1res=np.linspace(0,1,n)
        w3res=(1-w1res)
        w2res=np.zeros(n)
        # wres[idx,:]=w2res
        # wres[idxvec[idxvec!=idx][0],:]=w1res.T
        # try:
        #     wres[idxvec[idxvec!=idx][1],:]=w3res.T
        # except:
        #     pass
        wres=self.idx_fun(wres,np.vstack((w2res,w1res,w3res)),idx)
        xres=[]
        for i in range(n):
            xres.append(self.InsertExpression(wres[:,i],self.wi,self.xieq))
            Jscan=self.InsertExpression(cs.vertcat(xres[i],self.eta,self.XAi,self.Zsym),self.varnlp,J)
            res,success=self.SolveNLP(self.varscan,liq,Jscan)
            x_opt=res
            sol.append(cs.vertcat(xres[i],x_opt))


        rhoASD=self.PyCSAFToutput(sol,self.rhosys,self.varnlp).T[0]
        self.rho0API0=rhoASD[0]
        self.rho0Poly0=rhoASD[-1] if self.rho0Poly0==None else self.rho0Poly0

        rhoid=((1-w1res)/self.rho0API0+w1res/self.rho0Poly0)**-1
        rhoASDid=((1-w1res)/rhoASD[0]+w1res/rhoASD[-1])**-1
        VE=1/rhoASD-1/rhoASDid
        rhoASDNET=(1/rhoid+VE)**-1
        #fig,ax=plt.subplots()
        #ax.plot(w1res,rhoASDNET,'k-')
        #ax.plot(w1res,rhoid,'r-')
        self.rhoASD_fun=InterpolatedUnivariateSpline(w1res,rhoASDNET)

        return self.rhoASD_fun
    def rho0Pure(self,psys=0.03166*1E5*0.1,liq=True,T=298.15):

        self.PCSAFT(psys=psys,T=T) #PCSAFT call within PCSAFT call leads to problems
        J=cs.vertcat(self.Jp,self.Jx,self.JZ)
        sol=[]
        xres=[]
        for i in range(self.ncomp):
            wres=np.zeros((self.ncomp,self.ncomp))
            wres[i,i]=1
            xres.append(self.InsertExpression(wres[:,i],self.wi,self.xieq))
            Jscan=self.InsertExpression(cs.vertcat(xres[i],self.eta,self.XAi,self.Zsym),self.varnlp,J)
            res,success=self.SolveNLP(self.varscan,liq,Jscan)
            x_opt=res
            sol.append(cs.vertcat(xres[i],x_opt))
        rhopures=self.PyCSAFToutput(sol,self.rhosys,self.varnlp).T[0]
        self.rhopures=rhopures


    def GetExpData(self):
        try:
            polymer=self.scomp[0]
            #polymer="pvpva64" if polymer=="pvpva" else polymer
            api=self.scomp[2]
            #api="indometacin" if api=="indomethacin" else api
            wPolyASDstr=str(round(self.wPolyASD*100)/100) if not int(self.wPolyASD)==round(self.wPolyASD*100)/100 else str(int(self.wPolyASD))
            molec=polymer+"_"+api+"_"+wPolyASDstr if self.wPolyASD<0.99 else polymer
            molec=api if self.wPolyASD<0.01 else molec

            path1=join(cwd,"LeastFit","Fit",molec+".xlsx")
            self.molec=molec
            if self.Lit>0:
                raise Exception('Experimental Data might be available but literature Data was required')
            Data=pd.read_excel(path1).dropna(how="all")
            self.wexp=np.asarray(Data["w2ave/- 0"].dropna().astype(np.float64))
            self.RHexp=np.asarray(Data["RHave/- 0"].dropna().astype(np.float64))
            self.wstd=np.asarray(Data["w2/- 0"].dropna().astype(np.float64))
            self.RHstd=np.asarray(Data["RHstd/- 0]"].dropna().astype(np.float64))
        except:
            print("No experimental Data found")
            #print("Try searching for literature Data found")

            #try:
            #    path1=join(cwd,"Isotherm",molec+".csv")
            #    Data=pd.read_csv(path1).dropna(how="all")
            #    self.wexp=np.asarray(Data[Data.keys().tolist()[1]].dropna().astype(np.float64))
        #        self.RHexp=np.asarray(Data[Data.keys().tolist()[0]].dropna().astype(np.float64))
        #    except:
        #        self.wexp=[]
        #        self.RHexp=[]



    def Isotherm(self,wPolyASD,T=298.15,niso=30,NET=False):
        #PVPFilm=Mixture(*args,T=298.15,rho0Poly0=rho0Poly00,ksw=ksw0,wPolyASD=wPolyASD)
        idxp=self.idxp
        idx=self.idx
        self.wPolyASD=wPolyASD
        self.pH2OLV=self.PureVaporPressure(T=T) if not hasattr(self, 'pH2OLV') else self.pH2OLV
        #self.rhoASD(psys=self.pH2OLV,T=T)# if not hasattr(self, 'rhoASD_fun') else None
        self.rho0Pure(psys=self.pH2OLV,T=T) if not hasattr(self, 'rhopures') else None
        psys=np.linspace(self.pH2OLV*0.01,self.pH2OLV,niso) # Changed to 03.06.20 checkingwhat happens with VLE RH>1
        psysnet=np.linspace(self.pH2OLV*0.01,self.pH2OLV,niso)

        sol=[self.VLE(psys=val,T=T) for i,val in enumerate(psys)]
        self.ReInit()
        solnet=[self.VLE(psys=val,T=T,NET=1) for i,val in enumerate(psysnet)] if NET else sol
        self.ReInit()
        sucnet=np.asarray([val["success"] for i,val in enumerate(solnet)])
        resnet=np.asarray([val["wi"][idx] for i,val in enumerate(solnet)])
        res=np.asarray([val["wi"][idx] for i,val in enumerate(sol)])
        suc=np.asarray([val["success"] for i,val in enumerate(sol)])
        stab=np.asarray([val["Stability"] for i,val in enumerate(sol)])
        rho=np.asarray([val["rho"] for i,val in enumerate(sol)])

        lnxes=np.asarray([cs.log(val["xi"][idx]) for i,val in enumerate(sol)])
        lnxesnet=np.asarray([cs.log(val["xi"][idx]) for i,val in enumerate(solnet)])
        RHGAMMAPHI=np.asarray([np.exp((val["lngammai"][idx]))*val["xi"][idx] for i,val in enumerate(sol)])
        #RHGAMMAPHINET=np.asarray([np.exp((val["lngammai"][idx]))*(val["xi"][idx] for i,val in enumerate(solnet)])

        #lnpphies=np.asarray([val["lnphii"][1]+cs.log(val["RH"]) for i,val in enumerate(sol)])
        lnmues=np.asarray([val["mui"][idx] for i,val in enumerate(sol)])
        lnmuesnet=np.asarray([val["mui"][idx] for i,val in enumerate(solnet)])
        TH=1+np.gradient(lnmues-lnxes,lnxes)
        THnet=1+np.gradient(lnmuesnet-lnxesnet,lnxesnet)
        #sol["TH"],solnet["TH"]=TH,TH
        for i,val in enumerate(sol):
            sol[i]["TH"]=TH[i]
        for i,val in enumerate(solnet):
            solnet[i]["TH"]=THnet[i]
        wPolyASDstr=str(round(self.wPolyASD*100)/100) if not int(self.wPolyASD)==round(self.wPolyASD*100)/100 else str(int(self.wPolyASD))
        typ="Isotherm" if niso>=100 else "VLE"
        #self.WriteToExcel([sol,solnet],["SAFT","NET"],Datatype=typ+"_wpF_"+wPolyASDstr)
        self.GetExpData()
        #import Thermoplot
        yla="$w_w$[-]"
        xla="$p_{w}/p_{0w}^{LV}$ [-]"
        fig,ax=plt.subplots()

        ax.set_xlabel(xla)
        ax.set_ylabel(yla)


        def SplineCorrection(x,y,suc):
            return InterpolatedUnivariateSpline(x[suc==1],y[suc==1])(x)
        wwSAFT=SplineCorrection(psys/self.pH2OLV,res,suc)
        wwNET=SplineCorrection(psys/self.pH2OLV,resnet,sucnet)

        ax.plot(psys/self.pH2OLV,res,'--',label="PC-SAFT")
        ax.plot(psysnet/self.pH2OLV,resnet,'-.',label="NET-GP")
        ax.plot(psys/self.pH2OLV,wwSAFT,'--',label="PC-SAFT")
        ax.plot(psysnet/self.pH2OLV,wwNET,'-.',label="NET-GP")
        ax.plot(RHGAMMAPHI,wwSAFT,'k-',label="Gamma-Phi-Konzept")
        ax.plot(self.RHexp,self.wexp,'kx') if hasattr(self,"RHexp") else None
        ax.legend()

        self.res=wwSAFT
        self.resnet=wwNET
        self.resges=np.hstack((wwNET[wwSAFT<wwNET],wwSAFT[wwSAFT>wwNET]))
        #self.wg=self.GordonTaylorW(self.wPolyASD)
        self.RHsim=(psys/self.pH2OLV)
        #self.wSLE=np.asarray(wSLE)
        #self.RHSLE=np.asarray(RHSLE)

        ##################### For Flory Huggins
        # from FloryHuggins import FloryFilm
        # w2=np.linspace(0,wwSAFT[-1][-1],100)
        # F1=FloryFilm(self.comp[2],self.comp[1],wwSAFT[-1][-1])
        # RHFlory=F1.FloryRH(w2, 2.272)
        # ax.plot(RHFlory,w2,'r-')
        # THF=F1.THFaktor_fun(w2)
        # gammaf=F1.gamma2_fun(w2)
        # x2=RHFlory/gammaf
        # muf=np.log(x2*gammaf)

        # wiF=np.asarray([1-w2,w2,np.zeros_like(w2)])
        # base=[{"success":np.asarray([1]),"wi":np.asarray([1-val,val,np.zeros_like(val)]), "RH":RHFlory[i],"mui":np.asarray([muf[i],muf[i],muf[i]]),"TH":THF[i]} for i,val in enumerate(w2)]
        # self.WriteToExcel([sol,solnet,base],["SAFT","NET","Flory"],Datatype=typ+"_wpF_"+wPolyASDstr)
        # ################################
        self.WriteToExcel([sol,solnet],["SAFT","NET"],Datatype=typ+"_wpF_"+wPolyASDstr)
        return sol,solnet
    def MultIsothermOrigin(self,T=298.15):
        idx=self.idx
        idxp=self.idxp
        wPolyASDvec=np.asarray([1,0.8,0.5,0])
        wPolyASDvec=np.asarray([1,0.9,0.8,0.7,0.6,0.5]) if self.Lit else wPolyASDvec
        wPolyASDvec=np.asarray([0.9,0.8,0.7,0.6]) if self.Lit else wPolyASDvec
        wPolyASDvec=np.asarray([1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]) if self.Lit else wPolyASDvec
        #### Origin Plot of ISothermal Data###
        resvec=[]
        resnetvec=[]
        wexpvec=[]
        RHexpvec=[]
        RHsimvec=[]
        wstdvec=[]
        RHstdvec=[]
        SLE=[]
        wgvec=[]
        for i,val in enumerate(wPolyASDvec):
            self.Isotherm(val,T=T)
            self.wstd=np.zeros_like(self.wexp) if self.Lit else self.wstd
            self.RHstd=np.zeros_like(self.RHexp) if self.Lit else self.RHstd
            resnetvec.append(self.res*100) # To Include SAFT and NET-SAFT
            resnetvec.append(self.resnet*100)
            RHexpvec.append(self.RHexp*100)
            wexpvec.append(self.wexp*100)

            RHexpvec.append(self.RHsim*100)  # To Include instable points
            wexpvec.append(self.resges*100)

            RHsimvec.append(self.RHsim*100)
            RHsimvec.append(self.RHsim*100)
            wstdvec.append(self.wstd*100)
            RHstdvec.append(self.RHstd*100)
            wstdvec.append(np.zeros_like(self.wstd))
            RHstdvec.append(np.zeros_like(self.RHstd))
            #wgvec.append(self.wg*100)

        #import OriginVorlagePlot#added 03.03.2020 for origin

        xlab='RH[%]'
        ylab='w\-(w)[%]'
        x=RHexpvec
        y=wexpvec
        xe=wstdvec
        ye=RHstdvec
        xs=RHsimvec
        ys=resnetvec
        # wgx=np.asarray([])
        # wgy=np.asarray([])
        # for i,val in enumerate(wPolyASDvec):
        #     wgy=np.hstack((wgy,np.asarray([wgvec[i],wgvec[i],'-'])))
        #     wgx=np.hstack((wgx,np.asarray([0,100,'-'])))
        # ys.append(wgy)
        # xs.append(wgx)
        xData,yData,xSim,ySim,xErr,yErr=MusterDictOrigin()
        xu="%"
        yu="%"
        xc=xlab
        yc=ylab
        xData["values"],yData["values"],xSim["values"],ySim["values"],xErr["values"],yErr["values"]=x,y,xs,ys,xe,ye
        xData["unit"],yData["unit"],xSim["unit"],ySim["unit"],xErr["unit"],yErr["unit"]=xu,yu,xu,yu,xu,yu
        xData["comment"],yData["comment"],xSim["comment"],ySim["comment"],xErr["comment"],yErr["comment"]=xc,yc,xc,yc,xc,yc
        Plot(xData,yData,xSim,ySim,xErr,yErr,xlab,ylab,'GGW'+self.molec,Vorlage="Solubilities_Vorlage.opj")
        return

    def GordonTaylorW(self,wPolyASD,T=298.15):
        Tg0=np.asarray([val["Tg"] for i,val in enumerate(self.comp)])
        rho0=np.asarray([val["rho0"] for i,val in enumerate(self.comp)])
        wF3=1-wPolyASD
        k1=Tg0[1]*rho0[1]/(Tg0[2]*rho0[2])
        k2=Tg0[1]*rho0[1]/(Tg0[0]*rho0[0])
        wg=(k1*wF3*(Tg0[2]-T)+k2*(1-wF3)*(Tg0[0]-T))/((T-Tg0[1])+wF3*k1*(Tg0[2]-T)+(1-wF3)*k2*(Tg0[0]-T))
        self.wg=wg
        return wg
    def WasserFeed(self,T=298.15,RH=0.9,NET=1):
        idx=self.idx
        idxp=self.idxp
        wPolyASDvec=np.linspace(0,1,30)
        self.pH2OLV=self.PureVaporPressure(T=T) if not hasattr(self, 'pH2OLV') else self.pH2OLV
        psys=self.pH2OLV*RH
        res=[]
        for i,val in enumerate(wPolyASDvec):
            self.ReInit()
            self.wPolyASD=val
            sol=self.VLE(psys=psys,T=T,NET=NET)
            res.append(sol["wi"][idx])

        fig,ax=plt.subplots()
        ax.plot(1-wPolyASDvec,np.asarray(res),'ko')
        ax.set_xlabel("DL/-")
        ax.set_ylabel("ww/-")
        return 1-wPolyASDvec,np.asarray(res)

    def Isohumide(self,T=298.15,main=True,NET=False,RHsup=np.asarray([0]),ASDsup=np.asarray([0])):
        idx=self.idx
        idxp=self.idxp
        wPolyASDvec=np.linspace(0,1,self.nASD) if (ASDsup==0).all() else ASDsup#[::sup
        resvec=[] #wPolyASDxRhvecxncomp

        self.pH2OLV=self.PureVaporPressure(T=T) if not hasattr(self, 'pH2OLV') else self.pH2OLV
        self.rho0Pure(psys=self.pH2OLV,T=T) if not hasattr(self, 'rhopures') else None
        VLEexpvec=[]
        for i,vali in enumerate(wPolyASDvec):
            self.wPolyASD=vali

            RHvec=np.linspace(0.05,1,self.nVLE) if (RHsup==0).all() else RHsup
            #RHvec=
            psys=self.pH2OLV*RHvec

            sol=[self.VLE(psys=valj,T=T) for j,valj in enumerate(psys)]
            res=[val["wi"] for i,val in enumerate(sol)]
            self.ReInit()
            if NET==True:
                solnet=[self.VLE(psys=valj,T=T,NET=1) if vali!=0 else self.VLE(psys=valj,T=T) for j,valj in enumerate(psys)]
                resnet=[val["wi"] for i,val in enumerate(solnet)]
                res=np.asarray([res[j] if res[j][1]>resnet[j][1] else resnet[j] for j,valj in enumerate(psys)])
            #res=np.hstack((resnet[resnet>res],res[resnet<res]))

            resvec.append(res)
            self.GetExpData()
            if len(self.wexp)>0:
                w1exp=(self.wPolyASD)*(1-self.wexp)
                w2exp=self.wexp
                w3exp=(1-self.wPolyASD)*(1-self.wexp)
                RHexp=self.RHexp
                wexp=np.asarray([w1exp,w2exp,w3exp])
                VLEexpvec.append(wexp)
                self.wexp=[]
        resvec=np.asarray(resvec)
        VLE=[resvec[:,i,:].T for i, val in enumerate(resvec[0,:,0])]
        VLEexp=VLEexpvec#[VLEexpvec[:,i,:].T for i, val in enumerate(VLEexpvec[0,:,0])]
        TernaryPlot(VLE=VLE,VLEexp=VLEexp,molecules="Isohumide") if main else None
        return VLE,VLEexp
    def GordonTaylor(self,T=298.15):
        wPolyASDvec=np.linspace(0,1,self.nASD)
        wTgvec=[]
        for i,vali in enumerate(wPolyASDvec):
            wTg=self.GordonTaylorW(vali)
            wTgvec.append(np.asarray([vali*(1-wTg),wTg,(1-vali)*(1-wTg)]))
        TG=np.asarray(wTgvec)
        return TG
                #Low Vapor Pressures Stability determinant is smaller zero hence to few DOF
    def MiscibilityGap(self,T=298.15,main=True,name="name"):
        idxp=self.idxp
        idx=self.idx
        self.pH2OLV=self.PureVaporPressure(T=T) if not hasattr(self, 'pH2OLV') else self.pH2OLV
        LLE1=[]
        LLE2=[]
        LLE1vw=np.asarray([0,0,0])
        LLE2vw=np.asarray([0,1,0])
        LLE1v={"wi":np.asarray([0,0,0])}
        LLE2v={"wi":np.asarray([0,1,0])}
        w1ubvec=np.linspace(0,0.6,self.nLLE)[::-1] #Reversing helps
        #w1ubvec=np.hstack((np.linspace(0,0.01,15),np.linspace(0.01,1,5)))
        xlimit=1E-2
        w1ubvec[-1]=xlimit # Makes minor diferences in the polymer rich phase instead if LLEbinary
        #smallws=10**np.linspace(np.log10(xlimit),np.log10(1E-3),8)
        #w1ubvec=np.hstack((w1ubvec,smallws))
        LLE1pinfo,LLE2pinfo=[],[]

        for i,vali in enumerate(w1ubvec):

            LLE1v,LLE2v=self.LLE(psys=self.pH2OLV,ws=vali)# if vali!=0 else self.LLEbinary(psys=self.pH2OLV)  #detL1 wont help here too fe DOF
            LLE1vw=LLE1v["wi"]
            LLE2vw=LLE2v["wi"]
            LLE1.append(LLE1v["wi"]) #if LLE1vw[2]>=LLE2vw[2]+0.0001 else None
            LLE2.append(LLE2v["wi"]) #if LLE1vw[2]>=LLE2vw[2]+0.0001 else None
            LLE1pinfo.append(LLE1v)
            LLE2pinfo.append(LLE2v)
        LLE1=np.asarray(LLE1).T
        LLE2=np.asarray(LLE2).T
        LLE1[np.isnan(LLE1)]=0
        LLE2[np.isnan(LLE2)]=0
        self.WriteToExcel([LLE1pinfo,LLE2pinfo],["Phase1","Phase2"],Datatype="LLE")
        Spinodal=[np.asarray(self.SpinL).T,np.asarray(self.SpinR).T]
        TernaryPlot(LLE=[LLE1,LLE2],Spinodal=Spinodal,molecules=["MiscibilityGap",name]) if main else None
        return LLE1,LLE2,Spinodal

    def SolubilityLine(self,psys=1E5,T=298.15,main=True):
        idxa=self.idxa
        #self.pH2OLV=self.PureVaporPressure(T=T) if not hasattr(self, 'pH2OLV') else self.pH2OLV
        wscan=np.linspace(0,1,self.nSLE)
        SLEnet=[self.SLE(psys=psys,wwASD=valj) for j,valj in enumerate(wscan)] #if SLEbool else VLE
        resSLE=np.asarray([val["wi"] for i,val in enumerate(SLEnet)])
        self.ReInit()
        resvecSLE=np.asarray(resSLE)
        SLE=np.asarray([resvecSLE[:,i] for i, val in enumerate(resvecSLE[0,:])])
        self.WriteToExcel([SLEnet,SLEnet],["SLE1","SLE2"],Datatype="SLE")
        TernaryPlot(SLE=SLE,molecules="SolubilityLine") if main else None
        return SLE
    def StabilityMap(self,psys=1E5,T=298.15,main=True):
        #self.pH2OLV=self.PureVaporPressure(T=T) if not hasattr(self, 'pH2OLV') else self.pH2OLV
        sol=self.GPlot(psys=psys)
        wG=np.asarray([val["wi"] for i,val in enumerate(sol)])
        wG=wG[:,::-1]
        GG=np.asarray([np.asarray([val["Stability"][0]]) for i,val in enumerate(sol)])
        GG[np.isnan(GG)]=1E17
        G=np.hstack((wG,GG)).T
        TernaryPlot(G=G,molecules="StabilityMap") if main else None
        return G

    def GibbsDiagramm(self,T=298.15):
        LLE1,LLE2,Spinodal=self.MiscibilityGap(T=T,main=False)
        VLE,VLEexp=self.Isohumide(T=T,main=False)
        SLE=self.SolubilityLine(T=T,main=False)
        TG=self.GordonTaylor(T=T)
        molecules=[val["name"] for i,val in enumerate(self.comp)]
        G=self.StabilityMap(T=T,main=False)
        TernaryPlot(VLE=VLE,SLE=SLE,LLE=[LLE1,LLE2],Tg=TG.T,molecules=molecules,G=G,VLEexp=VLEexp,Spinodal=Spinodal)
        return
    def CalcSpace(self,liq=True,psys=1E5,T=298.15,n=100):
        import PyCSAFT.DasDennisSpacing
        wSpace=DasDennisSpacing.DasDennis(n,self.ncomp)
        results=[self.ConcAsInput(psys=psys,T=T,wres=wSpace[:,i],liq=liq) for i,val in enumerate(wSpace[0,:])]
        self.WriteToExcel([results],["GammaScan"],Datatype="GammaScan")
    def CalcSpaceDrugload(self,liq=True,psys=1E5,T=298.15,n=10):
        if self.ncomp>2:
            DLvec=np.linspace(0,1,n+1)
        else:
            DLvec=np.asarray([1])

        ww=np.linspace(0,1,n+1)
        wamat=np.asarray([(1-ww)*val for i,val in enumerate(DLvec)])
        wavec=wamat.flatten()
        wwvec=np.asarray([ww for i,val in enumerate(DLvec)]).flatten()
        wpvec=1-wwvec-wavec
        if self.ncomp>2:
            wSpacetemp=np.asarray([wpvec,wwvec,wavec])
            wSpace=np.zeros_like(wSpacetemp)
            wSpace[self.idxp,:]=wSpacetemp[0,:]
            wSpace[self.idx,:]=wSpacetemp[1,:]
            wSpace[self.idxa,:]=wSpacetemp[2,:]
        else:
            wSpace=np.asarray([wavec,wwvec])

        results=[self.ConcAsInput(psys=psys,T=T,wres=wSpace[:,i],liq=liq) for i,val in enumerate(wSpace[0,:])]
        self.WriteToExcel([results],["GammaScan"],Datatype="GammaScan")
        return results
    def TgGibbs(self,psys=1E5,T=298.15,main=True):
        TG=self.GordonTaylor(T=T)
        molecules=[val["name"] for i,val in enumerate(self.comp)]
        TernaryPlot(Tg=TG.T,molecules=molecules)

def TW(Film,p,Tvec):
    Film.idxp=0
    Film.idxa=2
    w2vec=[]
    w1vec=[]
    for i,val in enumerate(Tvec):
        Phase1=Film.SLE(psys=p,T=val,wwASD=0)
        w1vec.append(Phase1["wi"][2])

    w1vec=np.asarray(w1vec)
    fig,ax=plt.subplots()
    ax.plot(w1vec,Tvec-273.15,'r-')
    ax.set_xlabel("wa/-")
    ax.set_ylabel("T/°C")
    return w1vec,Tvec-273.15
# def IsothermSwipe(pol,sol,api,kij):
#     PVPFilm1=Mixture(pol,sol,api,wPolyASD=1,dikij=kij)
#     PVPFilm1.Lit=1
#     wPolyASDvec=np.linspace(0.5,1,6)
#     figs,axs=zip(*[PVPFilm1.Isotherm(val) for i,val in enumerate(wPolyASDvec)])
#     from combineFigs import merge
#     figI,axI=merge(*tuple(figs))
#     return figI,axI

if __name__=="__main__":
    T=298.15
    pure, kij=get_pcpar.get_par(["pvpva64","water","ritonavir"],T=T)
    #
    # pure, kij=get_pcpar.get_par(["pvpva64","water","nifedipin"],T=T)
    pol,water,api=pure
    #kij["waterindomethacin"]=-0.025 old one from Fortran.
    #kij["indomethacinwater"]=-0.025
    # water["sigi"]=2.79707 #Camaretti Wasser
    # water["mi"]=1.204659
    # water["ui"]=353.95
    # water["kapi"]=0.04509
    # water["epsAiBi"]=2425.67
    kij["hpmcasnaproxen"]=-0.0289
    kij["naproxenhpmcas"]=-0.0289

    kij["waterindomethacin"]=-0.025
    kij["indomethacinwater"]=-0.025
    kij["waterritonavir"]=-0.0485 # An llevbare LLE trend für 25°C
    kij["ritonavirwater"]=-0.0485
    #kij["waterritonavir"]=-0.0605 # An Wasersorption 25°C
    #kij["ritonavirwater"]=-0.0605
    kij["waterpvpva64"]=-0.134 #vs 0.134
    kij["pvpva64water"]=-0.134
    kij["waterpvp"]=-0.128 #vs 0.134
    kij["pvpwater"]=-0.128
    #pol["rho0Poly0"]=1140
    #pol["ksw"]=0.06
    Film=Mixture(pol,water,api,dikij=kij)
    #Tvec=np.linspace(0+273.15,170+273.15,100)

    #res=[Film.SLE(T=val) for i,val in enumerate(Tvec)]
    #wnap=np.asarray([val["wi"][2] for i,val in enumerate(res)])
    #pd.DataFrame(Tvec-273.15,wnap*100).to_excel("fil.xlsx")
    #Film.CalcSpaceDrugload(psys=1E5,T=298.15,n=100)
    Film.Isotherm(1,niso=30,NET=0)

    plt.show()
    #Film.LLE(ws=0)
    #plt.show()
    # DL,ww=Film.WasserFeed(T=T,RH=0.9,NET=1)
    #Film.Isotherm(1,niso=100,NET=1)


    #wNET,RHsim=Film.resnet.T[0]*100,Film.RHsim.T[0]*100
    #wSA,RHsim,wnet=Film.res.T[0]*100,Film.RHsim.T[0]*100,Film.resnet.T[0]*100
    #cop=np.vstack((RHsim,wSA)).T
    #copNET=np.vstack((RHsim,wNET)).T
    #pd.DataFrame(wSA,RHsim).to_excel("fil.xlsx")
    #Film.Isotherm(0.8,niso=30,NET=True)


    # pol["rho0Poly0"]=1337.3391800196282
    # pol["ksw"]=0.244317136081968#2.437429751297531e-08


# #250621 Plot von muAPI für Andreas
    # LactoseFilm=Mixture(water,api,dikij=kij)
    # LactoseFilm.idxa=1
    # LactoseFilm.idx=0
    # SL=LactoseFilm.SLE(psys=1E5,wwASD=0,T=298.15)
    # wSLE=SL["wi"][0]
    # muSLE=SL["mui"][1]
    #
    #
    # pvec=np.linspace(1E2,1E5,1)
    # fig,ax=plt.subplots()
    # ax.plot(wSLE,muSLE,'ko',label="SLE")
    # for i,valj in enumerate(pvec):
    #     results=LactoseFilm.CalcSpaceDrugload(liq=True,psys=valj,T=298.15,n=100)
    #     wl,wr,wn,spin=LactoseFilm.Spinodal(valj)
    #
    #     mu=np.asarray([val["mui"][1] for i,val in enumerate(results)])
    #     w=np.asarray([val["wi"][0] for i,val in enumerate(results)])
    #     mu=np.nan_to_num(mu,neginf=-80)
    #     ax.plot(w,mu,label="p="+"{:.1e}".format(valj)+"/Pa")
    #     ax.set_xlabel("ww/-")
    #     ax.set_ylabel("muAPI/-")
    #     ax.plot([wl[1],wl[1]],[np.min(mu),np.max(mu)])
    #     ax.plot([wr[1],wr[1]],[np.min(mu),np.max(mu)])
    # x=[w,np.asarray([wl[1],wl[1]]),np.asarray([wr[1],wr[1]])]
    # y=[mu,np.asarray([np.min(mu),np.max(mu)]),np.asarray([np.min(mu),np.max(mu)])]
    # xs=x
    # ys=y
    # Plotlist(x,y,xs,ys,xe=0,ye=0,Origin=True,filename="default",xla="?",yla="?")
    # plt.show()
    # LactoseFilm=Mixture(pol,water,api,dikij=kij)
    # LactoseFilm.wPolyASD=0
    # VL,VLNET=LactoseFilm.Isotherm(wPolyASD=0,T=298.15,NET=0)
    #
    # wVLE=VL[-1]["wi"][1]
    # muVLE=VL[-1]["mui"][2]
    # ax.plot(wVLE,muVLE,'kx',label="VLE RH=1")
    # ax.plot(wl[1],muVLE,'ko',label="VLE RH=1")
    # ax.plot(wr[1],muVLE,'ko',label="VLE RH=1")
    # ax.legend()
  #250621 Plot von muAPI für Andreas

    # Tvec=np.linspace(0,900,50)
    # res=[]
    # for i,val in enumerate(Tvec):
    #     pure, kij=get_pcpar.get_par(["hpmcas","water","nifedipin"],T=val)
    #     pol,water,api=pure
    #     LactoseFilm=Mixture(pol,water,api,dikij=kij)
    #     res.append(LactoseFilm.Pure(psys=1E5,liq=True,com=2,T=val))
    # rhores=np.asarray([val["rho"] for i,val in enumerate(res)])
    # etares=np.asarray([val["eta"] for i,val in enumerate(res)])
    # fig2,ax2=plt.subplots()
    # ax2.plot(Tvec,rhores)
    # fig3,ax3=plt.subplots()
    # ax3.plot(Tvec,etares)
    # LactoseFilm.idx=1
    # LactoseFilm.idxp=0
    # LactoseFilm.idxa=2
    # #RHsup=np.asarray([0.4316,0.7425,0.8434,0.9358])
    # #LactoseFilm.Isohumide(T=298.15,main=True,NET=True,RHsup=RHsup)
    # #LactoseFilm.MiscibilityGap()
    # LactoseFilm.SolubilityLine()
    # Tvec=np.linspace(0,180,20)+273.15
    # w1vec,Tvec2=TW(LactoseFilm,1E5,Tvec)
    # LactoseFilm.Isotherm(1,NET=True)
    # LactoseFilm.Isotherm(0.9,NET=True)
    # LactoseFilm.Isotherm(0.8,NET=True)
    # LactoseFilm.Isotherm(0.75,NET=True)
    # LactoseFilm.Isotherm(0.7,NET=True)
    # LactoseFilm.Isotherm(0.6,NET=True)


    #LactoseFilm.pH2OLV=LactoseFilm.PureVaporPressure(T=T)
    #xbin=np.linspace(0,1,100)
    #ergebnis=[LactoseFilm.Gammabinär(psys=LactoseFilm.pH2OLV*0.56,T=T,xbin=val) for i,val in enumerate(xbin)]
    #result=LactoseFilm.SLE(psys=1E3,T=T)
    #LactoseFilm.WriteToExcel([ergebnis,ergebnis],["GammaScan","GammaScan"],Datatype="GammaScan")
    #LactoseFilm.Isotherm(T=T,wPolyASD=1)

    #ergenis=LactoseFilm.SolubilityLine(T=T)
    #LLE Berechung für pvpva+wasser+naproxen bricht für das NET-GP kij zusammen
    #ethanol and acetone as solvents
    #pol,water,api=pure
    #
    #pure=tuple([api,water,api])
    #api["rho0Poly0"]=api["rho0"]
    #rhopoly und ksw in parameter datenbank oder zusätzliche datei
    #figI,axI=IsothermSwipe(pvpva64,water,indometacin)
    #kij["waterpvp"],kij["pvpwater"]=-0.13,-0.13
    #pvp["rho0poly0"]=10000*0.152139843258244
    #pvp["ksw"]=10**-10.250345956
    #pvp["ksw"]=0.4/(3166**4)
    #pol["rho0Poly0"]=1318.4
    #pol["ksw"]=2E-8
    #pol["rho0Poly0"]=0.140664719*1E4
    #pol["ksw"]=10**-7.439689065
    #pol["rho0Poly0"]=0.131915931176997*1E4
    #pol["ksw"]=10**-7.842464988
    #pol["rho0Poly0"]=1337.3391800196282
    #pol["ksw"]=2.437429751297531e-08
    #pol["ksw"]=0#PVA good Fit
    #pol["rho0Poly0"]=1090
    #pvp["ksw"]=0.3/(3166)
    #pvp["rho0Poly0"]=1415
    #pvp["kapi"]=0.02
#    pure, kij=get_pcpar.get_par(["pvpva64","water","naproxen"],T=T)
#    kijvec=np.linspace(-.15647,-.130,2)
#    kijvec=[-.123]
#    for i,val in enumerate(kijvec):
#        kij["pvpva64water"]=val
#        kij["waterpvpva64"]=val
    # kij["waterpvpva"]=-0.138327172646232
    # kij["pvpvawater"]=-0.138327172646232
    #kij["pvpwater"]=-0.13069043
    #kij["waterpvp"]=-0.13069043
    # kij["pvpwater"]=-0.128 #-0.148 excelent for 30°C
    # kij["waterpvp"]=-0.128
    # kij["waterindomethacin"]=-0.025
    # kij["indomethacinwater"]=-0.025
    # kij["waterritonavir"]=-0.0605
    # kij["ritonavirwater"]=-0.0605
    #pol["rho0Poly0"]=1450
    #pure=(pol,water,api)

#        PVPFilm1=Mixture(*pure,dikij=kij)#.VLE(psys=3000,T=298.15)
#        PVPFilm1.MiscibilityGap(T=T,name=str(kij["pvpva64water"]))
    #LactoseFilm=Mixture(*pure,wPolyASD=1,dikij=kij)


    #PVPFilm1.VLE()
    #liquid=PVPFilm1.VLE(psys=3000,T=T,idxp=0,idx=1)
    #solid=Film.SLE(psys=3000,wwASD=.4,T=T)
    #PVPFilm1.SolubilityLine
    #PVPFilm1.Isotherm(0.8,T=T)
    #We first determine wether a component is a trace componenr
    #by asking xi<1E-17. Then we calculate the component by
    #lnxi=muiK-lnfi0

    # pure2, kij2=get_pcpar.get_par(["water","acetone","indomethacin"])
    # kij2["waterindomethacin"]=-0.025
    # kij2["indomethacinwater"]=-0.025
    # BinarySLE=Mixture(*pure2,wPolyASD=0,dikij=kij2)
    # BinarySLE.PCSAFT(psys=1.013E5,T=298.15)

    # dictL=BinarySLE.SLE(psys=1.013E5,wwASD=0,T=298.15,i=2)

    #PVPFilm1.Isotherm(1,idxp=0)
    #PVPFilm1.Isotherm(.5)
    #PVPFilm1.Isotherm(1,niso=30)
    # PVPFilm1.Isotherm(0.5,niso=30)
    # #PVPFilm1.Isotherm(0,niso=30)
    # wNET,RHsim=PVPFilm1.resnet.T[0]*100,PVPFilm1.RHsim.T[0]*100
    # wSA,RHsim,wnet=Film.res.T[0]*100,Film.RHsim.T[0]*100,Film.resnet.T[0]*100
    # cop=np.vstack((RHsim,wSA)).T
    # copNET=np.vstack((RHsim,wNET)).T
    #PVPFilm1.Isotherm(0.8,niso=100)
    #PVPFilm1.Isotherm(0.5,niso=100)

    #PVPFilm1.Lit=1
    #IsothermSwipe(*pure,kij)
    #PVPFilm1.Isotherm(0.05)
    #PVPFilm1.Isotherm(0)
    #PVPFilm1.Isotherm(0.2)
    #print(PVPFilm1.LLEbinary())
    #rhoASD=PVPFilm1.rhoASD()
    #IsothermSwipe(pvp,water,indometacin,kij)

    #Film.GibbsDiagramm(T=T)

    #PVPFilm1.Lit=1
    #Film.GibbsDiagramm(T=T)
   # PVPFilm1.Isotherm(1,T=T)
    #PVPFilm1.StabilityMap(T=T)
    #wl,wr,wn,spin=PVPFilm1.Spinodal(T=T,idx=2)
    #PVPFilm1.SolubilityLine(T=T,idx=2)
    #PVPFilm1.Isohumide(T=T)
    # import timeit
    # tic=timeit.default_timer()

    # toc=timeit.default_timer()
    # print(toc-tic)
    #247.29s
    # pure, kij=get_pcpar.get_par(["pvp","water"],T=T)
    # kij["waterpvp"]=-0.128
    # kij["pvpwater"]=-0.128
    # PVPFilm1=Mixture(*pure,dikij=kij)#.VLE(psys=3000,T=298.15)
    # Tvec=np.linspace(150,200,20)
    # def LLET(p,Film,Tvec):
    #     Film.idxp=2
    #     # Film.idx=1
    #     w2vec=[]
    #     w1vec=[]
    #     for i,val in enumerate(Tvec):
    #         Phase1,Phase2=Film.LLE(psys=p,T=val,ws=0)
    #         w1vec.append(Phase1["wi"][0])
    #         w2vec.append(Phase2["wi"][0])
    #     w1vec=np.asarray(w1vec)
    #     w2vec=np.asarray(w2vec)
    #     fig,ax=plt.subplots()
    #     ax.plot(w1vec,Tvec-273.15,'r-')
    #     ax.plot(w2vec,Tvec-273.15,'b-')
    #     ax.set_xlabel("wp/-")
    #     ax.set_ylabel("T/°C")
    # LLET(5000,PVPFilm1,Tvec)
    # Phase1,Phase2=PVPFilm1.LLE(psys=7169,T=T,ws=0)
    #print(Phase1)
    #print(Phase2)

    #PVPFilm1.Isotherm(1,T,idxp=0)
    # PVPFilm1.nGplot=50
    # from OriginVorlagePlot import TernaryPlot
    # sol=PVPFilm1.GPlot(psys=0.03166*1E5)
    # wG=np.asarray([val["wi"] for i,val in enumerate(sol)])
    # wG=wG[:,::-1]
    # GG=np.asarray([np.asarray([val["Stability"][0]]) for i,val in enumerate(sol)])
    # #GG=np.asarray([val["Stability"] for i,val in enumerate(sol)]) #Funktioniert scheint aber keine Druckabhängigkeit zu haben. EVtl. das Problem?
    # GG[np.isnan(GG)]=1E17
    # G=np.hstack((wG,GG)).T

    # TernaryPlot(G=G)

    #PVPFilm1.Isotherm(0.9,)
    #PVPFilm1.MultIsothermOrigin(T=T)
    #PVPFilm1.LLE(i=1)
