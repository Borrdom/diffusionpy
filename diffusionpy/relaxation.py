import numpy as np

def Tggt(wi,rhoi,Tgi,q):
    nc=rhoi.shape[0]
    qmat=np.zeros((nc,nc))
    qmat[np.triu_indices(nc, k=1)]=q
    Excess=np.asarray([np.sum(np.prod(np.meshgrid(wi[i,:],wi[i,:]),axis=0)*qmat) for i,val in enumerate(wi[:,0])])
    Ideal=np.sum(wi*1/rhoi,1)/np.sum(wi*1/rhoi/Tgi,1)
    return Ideal+Excess
    #Triangular numbers for binary interaction parameters


def MEOS(self,rho2II_his,rho2II,activity,rho2,jMS_fun,E0,tau0,exponent,L0,rho02,rho,wi,Tg0,rho0):
    nJ=int(self.nJ)
    E0=self.E0 #E modulus of first Spring
    tau0=self.tau0
    eta0=E0*tau0
    exponent=self.exponent #if hasattr(self,'exponent') else 0
    EJ=E0#cs.MX.ones(nJ)*E0 #Switched? E0 is maxwell E1 is eq 23.07
    etaJ=eta0
    tauJ=etaJ/EJ
    #PDE Parameters
    nz=20  #Anzahl an räumlichen Diskretisierungspunkten
    nt=300  #Anzahl an zeitlichen Diskretisierungspunkten
    # Define Casadi Variables
    sigmaJ=np.zeros((nz+1,nJ))
    sigma=np.sum((sigmaJ.T*EJ).T,axis=1)
    dsigma=np.diff(sigma)
    deltaz=L0/nz
    v2=1/rho02

    epsilon=v2*rho2II
    X2II=rho2II/rho
    w2II=X2II/(X2II+1)
    #epsilon=v2_fun(w2II)*rho2II
    #epsilon=(rho2II+self.rho0ASD)/rho_fun(w2II)-1

    WL=np.exp(-rho2II/(1+epsilon)*exponent)
    #WL=cs.exp(-w2II*exponent)
#____________________________

    Tg=Tggt(wi,rho0,Tg0,0)

    C1=exponent
    #C3=C1
    #freevolume=C2+C1*(Tg-Tg[0])+C3*(self.T-Tg)
    #WL=10**(-C1*(1/freevolume-1/C2))
    if False:
        Tswitch=self.Tswitch
        WL=10**(-(self.T-Tg)/(C2+self.T-Tg)*C1)
        WL=cs.conditional((self.T+C2-Tg)<=0,[WL],1E17)
        C3=self.C3/2.308/self.R
        #C3=self.C3
        C4=-(self.T-Tswitch)/(C2+self.T-Tswitch)*C1
        AR=10**(-C3*(1/(Tg)-1/Tswitch)+C4)
        WL=cs.conditional((Tg)>Tswitch,[WL],AR)
        self.WLF=np.fmin(10**(-(self.T-self.Tgspan)/(C2+self.T-self.Tgspan)*C1)*self.eta0,1E17)
        self.WLF[(self.T+C2-self.Tgspan)<=0]=1E17
        self.AR=10**(-C3*(1/(self.Tgspan)-1/Tswitch)+C4)*self.eta0
        self.WLFreal=np.fmin(10**(-(self.T-self.Tgreal)/(C2+self.T-self.Tgreal)*C1)*self.eta0,1E17)
        self.WLFreal[(self.T+C2-self.Tgreal)<=0]=1E17
        self.ARreal=10**(-C3*(1/(self.Tgreal)-1/Tswitch)+C4)*self.eta0
    wspan=5
    WLF=np.exp(-wspan*self.rho0ASD*exponent)
    AR=np.exp(-wspan*self.rho0ASD*exponent)

    
    wGGW1=5
    initialcond=-rho2+wGGW1+1/v2*np.log(activity2/activity)*RV/np.sum(EJ)

    #initialcond=-rho2+wGGW1+1/v2_fun(w2)*cs.log(activity2/activity)*RV/cs.sum1(EJ)
    def solve(): return 1

    initialrho2=solve(initialcond)
    activity2=1
    activity_fun=lambda a: a
    initialsigma=np.log(activity2/activity_fun(initialrho2))*RV
    sigmaJ_his=np.zeros((nz+1,nt+1,nJ))
    rho2II_his[nz,:]=initialrho2
    for J in range(nJ):
        sigmaJ_his[nz,:,J]=initialsigma*np.atleast_1d(E0)[J]/np.sum(E0)
    R=8.145
    T+298.15
    M2=18.015/1000.
    v2=1/1000.

    RV=R*T*1/M2*1/v2
    MDF=1/RV*dsigma/deltaz #Mechanical driving force
    a2,a1,r1,r2=1,2,3,4
    ma=(a2-a1)/(r2-r1)
    ma_fun= lambda a: ma(a)

    #etaWL=cs.fmin((etaJ@WL.T).T,self.etamax)
    etaWL=(etaJ@WL.T).T
    #etaWL=cs.fmin((etaJ@WL.T).T,(etaJ@AR.T).T)
    #tauWL=etaWL/EJ
    a2II=a2
    drho2dt_hist=np.sum(1/etaWL[-1,:].T*sigmaJ[-1,:].T*EJ**2/RV,axis=0)/(1/a2II[-1]*ma_fun(rho2II[-1])+np.sum(EJ)*v2/RV)
    #drho2dt_hist=cs.sum1(1/WL[-1]*1/tauJ*sigmaJ[-1,:].T*EJ/RV)/(1/a2II[-1]*ma_fun(rho2II[-1])+cs.sum1(EJ)*v2_fun(w2II[-1])/RV)
    dj2II=0.
    drhodtNF=np.stack(dj2II/deltaz,drho2dt_hist)
    dsigmaJdt=np.zeros(nz+1,nJ)
    for i in range(nJ):
        dsigmaJdt[:,i]=-1/etaWL[:,i]*sigmaJ[:,i]*EJ[i]+drhodtNF*v2
        #dsigmaJdt[:,i]=-1/tauJ[i]*sigmaJ[:,i]*1/WL+drhodtNF*v2_fun(w2II)

    sigmaJvec=np.reshape(sigmaJ,((nz+1)*nJ,1))
    dsigmaJdtvec=np.reshape(dsigmaJdt,((nz+1)*nJ,1))
    xvec=np.stack(rho2II,sigmaJvec)
    fvec=np.stack(drhodtNF,dsigmaJdtvec)

    x0sigma=np.stack(sigmaJ_his[:,0,:]/E0,((nz+1)*nJ,1))
    x0=np.stack(rho2II_his[:,0],x0sigma)
    return