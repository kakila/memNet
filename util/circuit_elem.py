'''
 Copyright (C) 2013 - Juan Pablo Carbajal

 This progrm is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.
'''

# Author: Juan Pablo Carbajal <ajuanpi+dev@gmail.com>

from __future__ import division

import sympy as sp
from sympy import sin,exp,log

class circ_elem ():

  def __init__(self):
    self.args = dict()
    self.vars = dict(alg={},dif={})
    self.params = {}
    self.equation = dict()

  def set_variables (self,kw={}):
    pass

  def states(self,kind=[]):
    if not kind:
      return self.vars
    elif kind[0].lower() == "a":
      return self.vars["alg"]
    elif kind[0].lower() == "d":
      return self.vars[dif]

  def dynamics(self):
    data = dict()
    for k,v in self.equation.iteritems():
      data[str(k)] = str(v)
    return data

  def parameters(self):
    return self.params

  def arguments(self):
    return self.args.keys()

class Resistor (circ_elem):

  def __init__(self,V="V",I="I",R="R"):
    circ_elem.__init__(self)
    self.args = dict(V=V,I=I,R=R)
    self.set_variables(self.args)

  def set_variables (self,args):
    self.params = {"R": args["R"]}

    V,I,R = [args[x] for x in ["V","I","R"]]
    V,I,R = sp.symbols(" ".join([V,I,R]))

    self.equation = {I: I*R-V}

class Capacitor (circ_elem):

  def __init__(self,V="V",I="I",C="C",q="q",Cr="Cr"):
    circ_elem.__init__(self)
    self.args = dict(V=V,I=I,C=C,q=q,Cr=Cr)
    self.set_variables(self.args)

  def set_variables (self,args):
    self.vars["dif"] = {"q": args["q"]}
    self.params = dict(C=args["C"], Cr=args["Cr"])

    V,I,C,q,Cr = [args[x] for x in ["V","I","C","q","Cr"]]
    V,I,C,q,Cr = sp.symbols(" ".join([V,I,C,q,Cr]))

    self.equation = {q: I, I: q - C*(V-Cr*I)}

class SchottkyTunnel_volatile (circ_elem):

  def __init__(self,V="V",I="I",w="w",a="a",b="b",d="d",g="g",f="f",l="l",n="n"):
    circ_elem.__init__(self)
    self.args = dict(V=V,I=I,w=w,a=a,b=b,d=d,g=g,f=f,l=l,n=n)
    self.set_variables(self.args)

  def set_variables (self,args):
    self.vars["dif"] = {"w": args["w"]}
    self.params = dict([(x,args[x]) for x in ["a","b","d","g","f","l","n"]])

    V,I,w,a,b,d,g,f,l,n = [args[x] for x in ["V","I","w","a","b","d","g","f","l","n"]]
    V,I,w,a,b,d,g,f,l,n = sp.symbols(" ".join([V,I,w,a,b,d,g,f,l,n]))

    sinh = sp.Function("sinh")
    exp  = sp.Function("exp")

    schottky = a * (1 - exp (-b*V))
    tunnel  = g * sinh (d*V)
    self.equation = {\
    w: l * sinh (n*V) - f*w, \
    I: (1-w) * schottky + w * tunnel}

class HPsin_volatile (circ_elem):

  def __init__(self,V="V",I="I",r="r",m="m",R="R", w="w",C="C",f="f",e="e",s="s",l="l"):
    circ_elem.__init__(self)
    self.args = dict(V=V,I=I,r=r,m=m,R=R, w=w,C=C,f=f,e=e,s=s,l=l)
    self.set_variables(self.args)

  def set_variables (self,args):
    self.vars["dif"] = {"w": args["w"]}
    self.params = dict([(x,args[x]) for x in ["m","R","r","C","f","e","s","l"]])

    V,I,r,m,R,w,C,f,e,s,l = [args[x] for x in ["V","I","r","m","R","w","C","f","e","s","l"]]
    V,I,r,m,R,w,C,f,e,s,l = sp.symbols(" ".join([V,I,r,m,R,w,C,f,e,s,l]))

    PI = sp.Symbol("M_PI")
    power = sp.Function("pow")

    wt   = sp.tanh(w)
    ws   = power (wt,s)
    sin3 = power (sin(PI*ws), 3)
    w5   = power (wt,5)

    M = R - (R-r)/2 * (sp.tanh(w/2 + C)+1)
    self.equation = {w: m*r*I - f*power(sin3 + l*w5, e), I: I*M-V}

class HPtanh_volatile (circ_elem):

  def __init__(self,V="V",I="I",r="r",m="m",R="R",f="f", w="w",s="s",e="e",l="l",wz="wz",C="C"):
    circ_elem.__init__(self)
    self.args = dict(V=V,I=I,r=r,m=m,R=R,l=l, w=w,s=s,e=e,f=f,wz=wz,C=C)
    self.set_variables(self.args)

  def set_variables (self,args):
    self.vars["dif"] = {"w": args["w"]}
    self.params = dict([(x,args[x]) for x in ["m","R","r","l","s","e","f","wz","C"]])

    V,I,r,m,R,w,l,e,s,f,wz,C = [args[x] for x in ["V","I","r","m","R","w","l","e","s","f","wz","C"]]
    V,I,r,m,R,w,l,e,s,f,wz,C = sp.symbols(" ".join([V,I,r,m,R,w,l,e,s,f,wz,C]))

    M = R - (R-r)/2 * (sp.tanh(w/2 + C)+1)

    fw   = (1-l)*(sp.tanh(e*(w-wz-s))/2+1) - (sp.tanh(e*(w-wz+s))/2 +1) + l/2
    self.equation = {\
    w: m*r*I - f*w, \
    I: I*M-V}

class HP_TiO2 (circ_elem):

  def __init__(self,V="V",I="I",r="r",m="m",R="R", w="w",l="l"):
    circ_elem.__init__(self)
    self.args = dict(V=V,I=I,r=r,m=m,R=R, w=w,l=l)
    self.set_variables(self.args)

  def set_variables (self,args):
    self.vars["dif"] = {"w": args["w"]}
    self.params = dict([(x,args[x]) for x in ["m","R","r","l"]])

    V,I,r,m,R,w,l = [args[x] for x in ["V","I","r","m","R","w","l"]]
    V,I,r,m,R,w,l = sp.symbols(" ".join([V,I,r,m,R,w,l]))

    M = R - (R-r)/2 * (sp.tanh(m*w/2)+1)
    self.equation = {w: I - l*(w + 10.0/m), I: I*M-V}

class HP_TiO (circ_elem):

  def __init__(self,V="V",I="I",r="r",m="m",R="R", w="w",l="l"):
    circ_elem.__init__(self)
    self.args = dict(V=V,I=I,r=r,m=m,R=R, w=w,l=l)
    self.set_variables(self.args)

  def set_variables (self,args):
    self.vars["dif"] = {"w": args["w"]}
    self.params = dict([(x,args[x]) for x in ["m","R","r","l"]])

    V,I,r,m,R,w,l = [args[x] for x in ["V","I","r","m","R","w","l"]]
    V,I,r,m,R,w,l = sp.symbols(" ".join([V,I,r,m,R,w,l]))

    positive = sp.Function("positive")
    ws = positive(w)
    M = r*ws + R*(1-ws)
    self.equation = {w: m*I*ws*(1-ws)-l*ws, I: I*M-V}

class HP_TiOabs (circ_elem):

  def __init__(self,V="V",I="I",r="r",m="m",R="R", w="w",l="l"):
    circ_elem.__init__(self)
    self.args = dict(V=V,I=I,r=r,m=m,R=R, w=w,l=l)
    self.set_variables(self.args)

  def set_variables (self,args):
    self.vars["dif"] = {"w": args["w"]}
    self.params = dict([(x,args[x]) for x in ["m","R","r","l"]])

    V,I,r,m,R,w,l = [args[x] for x in ["V","I","r","m","R","w","l"]]
    V,I,r,m,R,w,l = sp.symbols(" ".join([V,I,r,m,R,w,l]))

    positive = sp.Function("positive")
    fabs = sp.Function("fabs")
    ws = positive(w)
    M = r*ws + R*(1-ws)
    self.equation = {w: m*fabs(I)*ws*(1-ws)-l*ws, I: I*M-V}

class MemCapac (circ_elem):
  """HP TiO memristor model.
   Voltage-driven memristor with parallel capacitance.
  """

  def __init__(self,V="V",I="I",r="r",m="m",R="R", w="w",l="l", C="C",Cr="Cr",q="q",mI="mI"):
    circ_elem.__init__(self)
    self.args = dict(V=V,I=I,r=r,m=m,R=R, w=w, l=l,C=C,Cr=Cr,q=q,mI=mI)
    self.set_variables(self.args)

  def set_variables (self,args):
    self.vars["alg"] = {"mI": args["mI"]}
    self.vars["dif"] = dict(w=args["w"],q=args["q"])
    self.params = dict([(x,args[x]) for x in ["m","R","r","l","C","Cr"]])

    V,I,r,m,R,w,l,C,q,mI,Cr = [args[x] for x in ["V","I","r","m","R","w","l","C","q","mI","Cr"]]
    V,I,r,m,R,w,l,C,q,mI,Cr = sp.symbols(" ".join([V,I,r,m,R,w,l,C,q,mI,Cr]))

    M = r*w + R*(1-w)
    self.equation = { \
    w: m*mI*w*(1-w) - l*w, \
    q: I - mI, \
    I: q - C*(V-Cr*(I-mI)), \
    mI: mI*M-V }

class MemCapacS (circ_elem):
  """HP TiO memristor model.
   Voltage-driven memristor with series capacitance.
  """

  def __init__(self,V="V",I="I",r="r",m="m",R="R", w="w",l="l", C="C",q="q",mV="mV"):
    circ_elem.__init__(self)
    self.args = dict(V=V,I=I,r=r,m=m,R=R, w=w, l=l,C=C,q=q,mV=mV)
    self.set_variables(self.args)

  def set_variables (self,args):
    self.vars["alg"] = {"mV": args["mV"]}
    self.vars["dif"] = dict(w=args["w"],q=args["q"])
    self.params = dict([(x,args[x]) for x in ["m","R","r","l","C"]])

    V,I,r,m,R,w,l,C,q,mV = [args[x] for x in ["V","I","r","m","R","w","l","C","q","mV"]]
    V,I,r,m,R,w,l,C,q,mV = sp.symbols(" ".join([V,I,r,m,R,w,l,C,q,mV]))

    M = r*w + R*(1-w)
    self.equation = { \
    w: m*I*w*(1-w) - l*w, \
    q: I, \
    I: q - C*(V-mV), \
    mV: I*M-mV}
