import sympy as sym
from sympy import Function,Symbol
from sympy import init_session,init_printing,latex,simplify,expand,factor,collect,cancel,trigsimp,pi,sqrt


p=sym.Symbol("p")
m=sym.Symbol("m")
o=sym.Symbol("Omega_e")
c=sym.Symbol("c")
wp=sym.Symbol("w_p")

theta=sym.Symbol("theta")
sint=sym.sin(theta)
cost=sym.cos(theta)

f=Function("f")(p,theta)

alphapp=Function("alpha_pp")(p,theta)
alphatp=Function("alpha_tp")(p,theta)
alphatt=Function("alpha_tt")(p,theta)
B=Function("B_w")(p,theta)


tau=Function("tau^-1")(p,theta)
vph=Function("v_ph")(p,theta)
v=p/m
tau_exp=pi/4*o**2*(vph-v*cost)/vph*B
wkplus_exp=1/2/o*(wp*p*cost/(c*m))**2+o+((wp**2*p*cost/(2*c**2*m*o))**2+(wp/c)**2)**(1/2)
wkminus_exp=1/2/o*(wp*p*cost/(c*m))**2+o-((wp**2*p*cost/(2*c**2*m*o))**2+(wp/c)**2)**(1/2)


alphapp_exp=tau*p**2*sint**3*vph**2
alphatp_exp=tau*p*sint**2*vph*(vph*cost-v)
alphatt_exp=tau*sint*(vph*cost-v)**2


a=alphapp*sym.diff(f,p)+alphatp*sym.diff(f,theta)
b=alphatp*sym.diff(f,p)+alphatt*sym.diff(f,theta)
dfdt= 1/p**2/sint*(sym.diff(a,p)+sym.diff(b,theta))
dfdt2=dfdt.subs({alphapp:alphapp_exp,alphatp:alphatp_exp,alphatt:alphatt_exp})

dalphapp_dp=sym.diff(alphapp_exp,p)
dalphatp_dp=sym.diff(alphatp_exp,p)
dalphatp_dtheta=sym.diff(alphatp_exp,theta)
dalphatt_dtheta=sym.diff(alphatt_exp,theta)

dtau_dp=sym.diff(tau_exp,p)
dwkplus_dp=sym.diff(wkplus_exp,p)
dwkminus_dp=sym.diff(wkminus_exp,p)
dtau_dtheta=sym.diff(tau_exp,theta)
dwkplus_dtheta=sym.diff(wkplus_exp,theta)
dwkminus_dtheta=sym.diff(wkminus_exp,theta)


def WriteToLaText(name,expression):
    out_file = open(name,"w")
    #out_file.write(latex(trigsimp(factor(simplify(expand(expression))))).replace("{\\left(p,\\theta \\right)}",""))
    out_file.write(latex(expression).replace("{\\left(p,\\theta \\right)}",""))
    out_file.close()

WriteToLaText("dfdt.txt",dfdt)
WriteToLaText("dalphapp_dp.txt",dalphapp_dp)
WriteToLaText("dalphatp_dp.txt",dalphatp_dp)
WriteToLaText("dalphatp_dtheta.txt",dalphatp_dtheta)
WriteToLaText("dalphatt_dtheta.txt",dalphatt_dtheta)

WriteToLaText("dtau_dp.txt",dtau_dp)
WriteToLaText("dwkplus_dp.txt",dwkplus_dp)
WriteToLaText("dwkminus_dp.txt",dwkminus_dp)
WriteToLaText("dtau_dtheta.txt",dtau_dtheta)
WriteToLaText("dwkplus_dtheta.txt",dwkplus_dtheta)
WriteToLaText("dwkminus_dtheta.txt",dwkminus_dtheta)


#out_file2 = open("dfdt2.txt","w")
#out_file2.write(latex(trigsimp(factor(simplify(expand(dfdt2))))).replace("{\left(p,\theta \right)}",""))
#out_file2.close()