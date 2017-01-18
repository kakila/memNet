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

from subprocess import check_call
from subprocess import CalledProcessError
from os import devnull

import codeCreator as cc

class codeMantainer ():

  def __init__(self,name, hpath="./include/", cpath="./src/"):
    self.hpath=hpath
    self.cpath=cpath
    p = dict(hpath=self.hpath,cpath=self.cpath)
    self.inp = cc.INPUT("in_"+name,**p)
    self.rhs = cc.RHS("rhs_"+name,**p)
    self.ics = cc.ICS("ics_"+name,**p)
    self.mn  = cc.MAIN("main_"+name,**p)

    self.name = dict(inputs=self.inp.name,
                     rhs=self.rhs.name,
                     ics=self.ics.name,
                     main=self.mn.name)

    self.mk = [] # placeholder for makefile
    self.updateMakefile()

    self.make = "make"
    #self.mk_options = [r'--file={}'.format(self.mk.name), r'-j7']
    self.mk_options = [r'--file={}'.format(self.mk.name)]

  def run(self):
    with open(devnull, "w") as f:
      check_call ([self.make]+self.mk_options+["all"], stdout=f)

    try:
      check_call ("./"+self.name["main"])
    except CalledProcessError as e:
      print "Oops, simulation failed.(code {})".format(e.returncode)
      print e.output
      return 1
    return 0

  def updateMakefile(self):
    self.mk = cc.Makefile(hpath=self.hpath,cpath=self.cpath, **self.name)
    self.mk.writeout()

  def updateInput(self,**kwargs):
    print "Writintg inputs to {}\n".format(self.name["inputs"])
    v = False
    if kwargs.has_key("verbose"):
      v = kwargs.pop("verbose")

    needed = dict.fromkeys(["inputs","params","vars"])
    parse_needed (needed, kwargs)

    body = dict (vars=needed["vars"],
                 inputs=needed["inputs"],
                 params=needed["params"])
    self.inp.writeout(body=body,
                      Ninputs=len(needed["inputs"]),
                      verbose=v,
                      params=needed["params"])

  def updateRHS(self,**kwargs):
    print "Writintg RHS to {}\n".format(self.name["rhs"])
    v = False
    if kwargs.has_key("verbose"):
      v = kwargs.pop("verbose")

    needed = dict.fromkeys(["rhs","params","vars","input_name"])
    parse_needed (needed, kwargs)

    self.rhs.writeout(body=needed,
                      inputsfile = self.inp.path["h"],
                      verbose=v,
                      params=needed["params"]
                     )

  def updateICS(self,**kwargs):
    print "Writintg initial conditions to {}\n".format(self.name["ics"])
    v = False
    if kwargs.has_key("verbose"):
      v = kwargs.pop("verbose")

    needed = dict.fromkeys(["ics","vars"])
    parse_needed (needed, kwargs)

    needed["params"]={}
    if kwargs.has_key("params"):
      needed["params"] = kwargs["params"]

    self.ics.writeout(body=needed,
                      verbose=v,
                      params=needed["params"]
                     )

  def updateMain(self,**kwargs):
    print "Writintg main file to {}\n".format(self.name["rhs"])
    v = False
    if kwargs.has_key("verbose"):
      v = kwargs.pop("verbose")

    tol = dict()
    wmode="w"
    if kwargs.has_key("tol"):
      tol = kwargs.pop("tol")
    if kwargs.has_key("wmode"):
      wmode = kwargs.pop("wmode")

    needed = dict.fromkeys(["Neq","outfile","tspan","vars","period"])
    parse_needed (needed, kwargs)

    tol["vars"] = needed["vars"]

    needed.pop("vars")
    self.mn.writeout(
                     rhsfile = self.rhs.path["h"],
                     wmode = wmode,
                     tolerances=tol,
                     verbose=v,
                     **needed
                    )

  def updateLib(self,opt=[],name=[]):
    if not opt:
      opt = self.mk_options
    if not name:
      #name = [v for k,v in self.name.iteritems() if k !="main"]
      with open(devnull, "w") as f:
        check_call ([self.make]+self.mk_options+["all"], stdout=f)
    else:
      libs = libify (name)
      for l in libs:
        with open(devnull, "w") as f:
          check_call ([self.make]+self.mk_options+[l], stdout=f)

  def cleanBuild(self,opt="clean"):
    with open(devnull, "w") as f:
      check_call ([self.make]+self.mk_options+[opt], stdout=f)

def libify(name):
  if type(name) is not list:
    name=[name]

  return ["lib"+n+".so" for n in name]

def parse_needed(needed,kw):
  ngiven = len(needed)
  for k,v in kw.iteritems():
    if k in needed.keys():
      needed[k] = v
      ngiven -=1
    else:
      raise ValueError ("Unknown key {} for inputs file".format(k))
  if ngiven:
    missing = [k for k in needed.keys() if not needed[k]]
    raise ValueError ("Needed key(s) {} not given.".format(", ".join(missing)))
