import casadi as cs
class Complex(object):
    def __init__(self, re,im):
        self.re = re
        self.im = im

    # adding two objects
    def __add__(self, o):
        return Complex(self.re + o.re,+self.im + o.im)
    def __sub__(self, o):
        return Complex(self.re - o.re,+self.im - o.im)
    def __mul__(self, o):
        return Complex((self.re * o.re)-(self.im*o.im),(self.re*o.im)+(self.im * o.re))
    def __truediv__(self, o):
        den=o.re * o.re + o.im * o.im
        num1=self.re * o.re+self.im * o.im
        num2=self.im * o.re-self.re * o.im
        return Complex(num1/den,+num2/den)
    def __pow__(self,o):
        abs=(self.re*self.re+self.im*self.im)**0.5
        abso=abs**o
        theta=cs.arctan(self.im/self.re)

        return Complex(abso*cs.cos(theta*o),abso*cs.sin(theta*o))

 

A=Complex(5,1)
B=Complex(1,2)
a=5+1j
b=1+2j


C=A*B
D=A/B
c=a*b
d=a/b
E=A**3
e=a**3
print(C.re-c.real)
print(C.im-c.imag)

print(D.re-d.real)
print(D.im-d.imag)

print(E.re-e.real)
print(E.im-e.imag)