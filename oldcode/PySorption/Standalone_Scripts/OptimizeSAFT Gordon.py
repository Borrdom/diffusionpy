RH=0.1
val=self.pH2OLV*RH
sol=self.VLE(psys=val,T=T)
wi=sol["wi"]
w2=wi[1]
J=-self.wPolyASD
w2g=self.GordonTaylor(self.wPolyASD)

g=w2-w2g

max(j) s.t g