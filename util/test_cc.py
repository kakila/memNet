import codeCreator as cc
import numpy as np

ipath="../include/"
cpath="../src/"

inp = cc.INPUT("input_test",hpath=ipath,cpath=cpath)
rhs = cc.RHS("rhs_test",hpath=ipath,cpath=cpath)
ics = cc.ICS("ics_test",hpath=ipath,cpath=cpath)
mn  = cc.MAIN("main_test",hpath=ipath,cpath=cpath)

var = ["$x","$y"]

print "**** Writing input file ****"
inputs = {"$Vs0":"$A*sin($w0*t)","$Vs1":"$A*sin($w1*t)"}
i_params = {"$A":1,"$w0":2*np.pi,"$w1":2*np.pi*3}
body = dict (vars=var, inputs=inputs, params=i_params)
inp.writeout(body=body,Ninputs=len(inputs),params=i_params)

print "**** Writing rhs file ****"
rhs_f = {"$x":["$y",True],"$y":["-$K*$x - $D*$y",True]}
r_params = {"$K":2*(2*np.pi)**2,"$D":0.1}
body = dict (vars=var, rhs=rhs_f, params=r_params)
rhs.writeout(body=body,params=r_params)

print "**** Writing ics file ****"
ics_v = {"$x":[1.0,0.0],"$y":[0.0,0.0]}
body = dict (vars=var, ics=ics_v, params={})
ics.writeout(body=body,params={})

print "**** Writing main file ****"
mn.writeout(Neq=len(rhs_f),
rhsfile = rhs.path["h"],
inputsfile = inp.path["h"],
outfile = "test_out.dat",
wmode = "w",
tspan= [0,1,5e-2],
tolerances={"vars":var}
)

print "**** Creating Makefile ****"
mk = cc.Makefile(inputs=inp.name,rhs=rhs.name,ics=ics.name,main=mn.name)
mk.writeout()
