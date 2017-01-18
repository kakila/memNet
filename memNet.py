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

import numpy as np
import networkx as nx

from util.codeMantainer import codeMantainer as cm
from util.circuit_elem import circ_elem

class memristorNetwork ():


    def __init__(self, name="memNet", **kwargs):

      self.name  = name;

      self.graph = nx.MultiDiGraph()

      self.vars  = dict(alg=[],dif=[],net_i=[],net_v=[],inp_v=[],inp_i=[]);

      self.equations = dict (KVL=[], KCL=[], Dyn={})

      self.parameters = dict (inputs={}, elems={})

      self.ics = dict()

      self.options = dict()

      self.tspan = []

      tmp = dict()
      if kwargs.has_key("hpath"):
        tmp["hpath"] = kwargs["hpath"]
      if kwargs.has_key("cpath"):
        tmp["cpath"] = kwargs["cpath"]

      self.__mantainer__ = cm(name,**tmp)

      self.isDAE = True

      self.var_order = ("dif","alg","net_v","net_i","inp_i")

      #### Initialization

      for k in self.__dict__.keys():
        if kwargs.has_key(k):
          method = getattr(self, "__proc_"+k, lambda x: x)
          method(kwargs[k])

    def run(self):
      fail = self.__mantainer__.run()
      if not fail:
        x = np.loadtxt(self.options["outfile"])
        t = x[:-1,0]
        data = dict(zip(self.getVars(),x[:-1,1:].T))
        res  = dict(zip(self.getVars(), \
                        [[y,yp] for y,yp in zip(x[-2,1:].T,x[-1,1:].T)] ))
      else:
        t=[]
        data=[]
        res=[]

      return t,data,res,fail

    def writeout(self,outfile="out.dat",verbose=False, kind=[],**kwargs):
      self.__mantainer__.updateMakefile()
      v = self.getVars()

      if kind=="inputs" or not kind:
        ## Inputs
        self.__mantainer__.updateInput(
                    vars=v,
                    inputs=self.getInputs(),
                    params=self.parameters["inputs"],
                    verbose=verbose
                    )
      if kind=="rhs" or not kind:
        ## Equations
        self.__mantainer__.updateRHS(
                  vars=v,
                  rhs=self.getEquations(),
                  params=self.parameters["elems"],
                  input_name=self.getInputs().keys(),
                  verbose=verbose
                 )
      if kind=="ics" or not kind:
        ## Initial conditions
        self.__mantainer__.updateICS(
                  vars=v,
                  ics=self.ics,
                  verbose=verbose
                 )

      if kind=="main" or not kind:
        self.options["outfile"] = outfile
        ## Main file
        self.__mantainer__.updateMain(
                   vars=v,
                   Neq=len(v),
                   outfile=outfile,
                   tspan=self.tspan,
                   period=self.options["period"],
                   verbose=verbose,
                   **kwargs
                  )

    def compile(self,name=[]):
      self.__mantainer__.updateLib(name)

    def cleanup(self,opt="clean"):
      self.__mantainer__.cleanBuild(opt)

    def draw (self):
      """ Draw the network assuming grid distribution """

      N = len([v for v in self.getVars("net_v") if v[1:3]=="Vo"])
      G = max(self.graph.nodes())
      nlayer = G/(G-N)
      nside  = N/(nlayer-1)

      pos_ = dict([(x[0],x[1]["pos"]) for x in self.graph.nodes(data=True) if x[0]>=0])
      var   = self.getVars(kind="net_v")
      ntype = {"in":[],"out":[],"gnd":[]}
      node_ = [x for x in self.graph.nodes() if x>=0]
      edge_ = self.graph.edges()

      # Classify nodes and edges
      for v in var:
        e = self.getVarEdge([v])[0]
        ii = []
        if e[0]>=0 and e[0] in node_:
          ii = e[0]
        elif e[1]>=0 and e[1] in node_:
          ii = e[1]

        if ii==[]:
          continue

        if "s" in v or "i" in v:
          # input
          ntype["in"].append(ii)
          node_.remove(ii)
          # Remove readout edges since resitance is really big
          if e[0]<0 or e[1]<0:
            edge_.remove(e)

        elif "o" in v:
          # output
          ntype["out"].append(ii)
          node_.remove(ii)
          # Remove readout edges since resitance is really big
          if e[0]<0 or e[1]<0:
            edge_.remove(e)

        elif "g" in v:
          # output
          ntype["gnd"].append(ii)
          node_.remove(ii)

      # Add new grounds for visualization
      G = nx.DiGraph()
      G.add_edges_from(edge_)
      resistors=[]

      for e in edge_:
        if e[0]==-1 and (e[1] in ntype["in"]):
          # Add connections to L_GND: inputs
          G.add_edge("L_GND",e[1])
          G.remove_edge(e[0],e[1])
          resistors.append(("L_GND",e[1]))
        if e[1]==-1 and (e[0] in ntype["gnd"]):
          # Add connections to R_GND: sinks
          G.add_edge(e[0],"R_GND")
          G.remove_edge(e[0],e[1])
          resistors.append((e[0],"R_GND"))
      edge_ = G.edges()

      node_ = ntype["in"] + ntype["out"] + ntype["gnd"] + ["L_GND","R_GND"]
      nn = len (node_)
      color_= [0]*nn
      for i,v in enumerate(node_):
          if v in ntype["in"]:
            color_[i] = "r"
          elif v in ntype["out"]:
            color_[i] = "y"
          elif v in ntype["gnd"]:
            color_[i] = "b"
      color_ += ["k"]*2

      pos_["L_GND"] = (-2,(nside-1)/2.0)
      pos_["R_GND"] = (nlayer+1,(nside-1)/2.0)

      G.remove_node(-1)
      nx.draw (G,pos=pos_,\
               edgelist=resistors, \
               arrows=False,\
               nodelist=[],width=2, \
               edge_color="b", node_color="k",with_labels=True)

      tmp = [e for e in edge_ if e not in resistors]

      nx.draw (G,pos=pos_,edgelist=tmp,\
               nodelist=node_,node_color=color_,width=1)


    ##### Setters ######
    def addElems (self, e, model=[]):

      # Add edges
      self.graph.add_edges_from (e)

      ed = self.graph.edges(data=True)

      if model:
        for x in model:
          if not isinstance (x,circ_elem):
            raise ValueError ("{}, model must be of type circ_elem.".format(x.__class__))

      if len(model) == 1:
        # Get root names for variables and parameters
        var = model[0].states()
        par = model[0].parameters()

        for u,v,d in self.graph.edges_iter(data=True):
          if (u,v) in e and not d:
            i = ed.index((u,v,d))

            # Rename model base variables
            alg = decoVar(var["alg"],i)
            dif = decoVar(var["dif"],i)
            V   = decoVar(["V"],i)
            I   = decoVar(["I"],i)



            # Parameters
            par = dict([(k,kd) for k,kd in zip(par, decoVar(par,i))])
            for k in par.values():
              self.parameters["elems"][k] = []

            vnew = dict(
                     [(k,x) for k,x in zip(var["alg"].keys()+var["dif"].keys(),alg+dif)] +\
                     [(k,x) for k,x in zip(["V","I"],V+I)] +\
                    par.items())

            # Update the model names
            model[0].set_variables(vnew)

            # Get renamed equations
            eqs = model[0].dynamics()

            # Update edge data
            d["alg"]   = alg
            d["dif"]   = dif
            d["eqs"]   = eqs
            d["net_v"] = V
            d["net_i"] = I
            d["pars"]  = par

            # Reset model base variables
            model[0].set_variables(model[0].args)

      elif len(model) == len(e):
        for u,v,d in self.graph.edges_iter(data=True):
          if (u,v) in e and not d:
            j = e.index((u,v))
            e.pop(j)
            m = model.pop(j)
            var = m.states()
            par = m.parameters()

            i = ed.index((u,v,d))

            # Rename model base variables
            alg = decoVar(var["alg"].values(),i)
            dif = decoVar(var["dif"].values(),i)
            V,I = m.args["V"],m.args["I"]
            V   = decoVar([V],i)
            I   = decoVar([I],i)

            # Parameters
            par = dict([(k,kd) for k,kd in zip(par.keys(), decoVar(par.values(),i))])
            for k in par.values():
              self.parameters["elems"][k] = []

            vnew = dict( \
                     [(k,x) for k,x in zip(var["alg"].keys()+var["dif"].keys(),alg+dif)] +\
                     [(k,x) for k,x in zip(["V","I"],V+I)] +\
                    par.items())

            # Update the model names
            m.set_variables(vnew)

            # Get renamed equations
            eqs = m.dynamics()

            # Update edge data
            d["alg"]   = alg
            d["dif"]   = dif
            d["eqs"]   = eqs
            d["net_v"] = V
            d["net_i"] = I
            d["pars"]  = par

            # Reset model base variables
            m.set_variables(m.args)

      else:
        raise ValueError ("You must give one model or as many models as elements.")


      # Update object variables
      self.__update_vars()

    ###
    def addInput (self, e, func=[],param={}):

      # Add elements
      self.graph.add_edges_from (e)
      ed = self.graph.edges(data=True)
      n = len(self.vars["inp_v"])
      i = n
      if len(func) == len(e):
        for u,v,d in self.graph.edges_iter(data=True):
          if (u,v) in e and not d:
            d["alg"] = decoVar(["Vs"],i)
            d["net_i"] = decoVar(["Is"],i)
            d["source"] = func[e.index((u,v))]
            i+=1
      else:
        raise ValueError ("You must give as many functions as inputs.")

      self.__update_vars()

      # TODO add parameters to the edge and in the book keeping
      # add them to the vector. Also the user would need to know the name.
      #write set parameter function
      if param:
        for k in param.keys():
          self.parameters["inputs"][k] = param[k]

    ##### Book keeping ######
    def __update_vars(self):
      self.vars["alg"]   = []
      self.vars["dif"]   = []
      self.vars["inp_v"] = []
      self.vars["inp_i"] = []
      self.vars["net_v"] = []
      self.vars["net_i"] = []
      for e in self.graph.edges_iter(data=True):
        ## Elem variables
        if len(e) >=2:
          # not sources
          if not e[2].has_key("source"):
            self.vars["alg"] += [x for x in e[2]["alg"] if x]
            self.vars["dif"] += [x for x in e[2]["dif"] if x]
            self.vars["net_v"]  += [x for x in e[2]["net_v"] if x]
            self.vars["net_i"]  += [x for x in e[2]["net_i"] if x]
          # sources
          else:
            self.vars["inp_v"] += [x for x in e[2]["alg"] if x]
            self.vars["inp_i"] += [x for x in e[2]["net_i"] if x]

    def save_pickle(self):
      from cPickle import dump
      from io import open
      f = open(self.name+".pkl",'wb')
      dump (self.__dict__,f,-1)
      print "Wrote file {}".format(f.name)
      f.close()

    def load_pickle(self,name=[]):
      from cPickle import load
      from io import open

      if not name:
        name=self.name+".pkl"

      f = open(name,'rb')
      d = load (f)
      print "Loaded file {}".format(f.name)
      f.close()
      for k,v in d.iteritems():
        setattr(self,k,v)

    ##### Equations ####
    def buildKirchhof(self):
      """Assemble Kirchhof equations for the network.
      """

      from util.cycle_space import chords, cycle_space_matrix, rref

      ### Kirchhof voltage law for the network
      V = cycle_space_matrix(self.graph,*chords(self.graph)).T
      if np.linalg.matrix_rank(V) < V.shape[0]:
        print "Reducing KVL"
        tmp = rref (V.copy())
        # Delete zero rows
        print "Before reduction ", V.shape
        M = np.delete (V, np.where((tmp==0).all(axis=1)),axis=0)
        print "Before reduction ", V.shape

      eqV = self.__buildKL__(V,"Voltage")

      ### Kirchhof current law for the network
      C = nx.incidence_matrix(self.graph, oriented=True).toarray()[:-1]
      eco = False
      if np.linalg.matrix_rank(C) < C.shape[0]:
        print "Reducing KCL"
        tmp = rref (C.copy())
        # Delete zero rows
        print "Before reduction ", C.shape
        M = np.delete (C, np.where((tmp==0).all(axis=1)),axis=0)
        print "Before reduction ", C.shape

      eqC = self.__buildKL__(C,"Current")

      return {"KVL":(V, eqV),"KCL":(C,eqC)}

    def __buildKL__ (self,M,kind):
      import sympy as sp

      if kind.lower() == "voltage":
        n = 0
      elif kind.lower() == "current":
        n = 1

      cc = ("v","i")
      nc = "net_"+cc[n]
      ic = "inp_"+cc[n]
      var = []
      for x,y,d in self.graph.edges(data=True):
        if d.has_key(nc):
          var.append(d[nc][0])
        elif d.has_key(ic):
          var.append(d[ic][0])
        elif d.has_key("source"):
          var.append(d["alg"][0])

      var = sp.symbols (var)

      eq = list()
      for row in M:
        ind = row.nonzero()[0]
        eq.append(np.sum([row[i]*var[i] for i in ind]).together())

      if n == 0:
        self.equations["KVL"] += [str(x) for x in eq]

      elif n == 1:
        if M.shape[0] == self.graph.number_of_nodes():
          print "Reduction failed, solving KCL explicitly. This may take some time ..."
          # Reduction failed but there is still too many equations
          # This is highly inefficient. TODO
          sol = sp.solve(eq,var)
          eq = [v-k for k,v in sol.iteritems()]

        self.equations["KCL"] += [str(x) for x in eq]

      return eq

    ###
    def buildDynamics(self):
      """Assemble the dynamics equations of the network and the dictionary
      with the parameters for IDA solver in SUNDIALS.
      """

      for e in self.graph.edges_iter(data=True):
        if not e[2].has_key("source"):
          for k,v in e[2]["eqs"].iteritems():
            if k in e[2]["dif"]:
              l = [v,True]
            else:
              l = [v,False]

            self.equations["Dyn"][k] = l

    ##### Getters ######
    def getVars(self,kind=[]):
      if not kind:
        return [v for varsk in self.var_order for v in self.vars[varsk]]
      if kind[0] == "a":
        return self.vars["alg"]
      elif kind[0] == "d":
        return self.vars["dif"]
      elif kind[0] == "inp_v":
        return self.vars["inp_v"]
      elif kind[0] == "inp_i":
        return self.vars["inp_i"]
      elif kind == "net_i":
        return self.vars["net_i"]
      elif kind == "net_v":
        return self.vars["net_v"]
      else:
        raise ValueError("No variable of kind "+kind)

    ###
    def getVarEdge(self,var):
      if type(var) != type(list()):
        var = [var]

      for e in self.graph.edges_iter(data=True):
        val = (e[0],e[1])
        types = e[2].keys()
        for typ in types:
          if any(map(lambda x: x in e[2][typ], var)) :
            return val, typ
      return None,None

    ###
    def getEquations(self):
      eqs = dict(self.equations["Dyn"].items())

      # Add Kirchhof laws
      # KVL
      var = set(self.vars["net_v"]) - set(eqs.keys())
#      var = [x for x in self.getVars() if x not in eqs.keys()]
      for v in self.equations["KVL"]:
         eqs[var.pop()] = [v,False]

      # KCL
      var = set(self.getVars()) - set(eqs.keys())
#      var = [x for x in self.getVars() if x not in eqs.keys()]
      for v in self.equations["KCL"]:
        eqs[var.pop()] = [v,False]

      return eqs

    ###
    def getInputs(self,p=[]):
      inp = dict()
      for e in self.graph.edges_iter(data=True):
        if e[2].has_key("source"):
          inp[e[2]["alg"][0]] = e[2]["source"]
      return inp

#########################################
#####       Helper functions        #####
def save4octave (net, path="./"):
  cmd,v  = getIndices(net)
  params = net.parameters["elems"]
  opts   = net.options["problem"]

  txt = "# Octave script. Automatically generated.\n" + \
        'x = load ("{}");'.format(net.options["outfile"]) + "\n" + \
        "t = x(1:end-1,1);\n" + \
        "x = x(1:end-1,2:end);\n" + \
        "# Problem parameters\n" + \
        "\n".join( \
           ["{var} = {val};".format(var=k,val=v) for k,v in opts.iteritems() ] \
           ) + \
        "# Memristance\n" + \
        "mem_R  = {};\n".format([y for x,y in params.iteritems() if x[1]=="R" and not x[2].isalpha()]) + \
        "mem_r  = {};\n".format([y for x,y in params.iteritems() if x[1]=="r"]) + \
        "mem_m  = {};\n".format([y for x,y in params.iteritems() if x[1]=="m"]) + \
        "mem_C  = {};\n".format([y for x,y in params.iteritems() if x[1]=="C"]) + \
        "M =@(x) mem_R - (mem_R-mem_r).*(tanh (x/2 + mem_C) + 1)/2;\n" + \
        cmd

  f = open(path + 'load_data_{}.m'.format(net.name), 'w')
  f.write(txt)
  f.close()

def fillParams(keys,defval={}):
  '''
  Returns dictionary with parameter values
  Use defval to pass uniform values to all parameters.
  '''

  params = dict().fromkeys(keys)

  for k in keys:
    idx = min([k.find(x) for x in "0123456789" if k.find(x)>=0])
    params[k] = defval[k[1:idx]]

  return params

def getIndices(net, octave=True):
  '''
  Returns indices of input readout , outputs, internal states, and ground currents.
  If octave=True (deault) it returns a string that can be executed in octave to generate
  the correspoinding indices.
  '''

  val = dict().fromkeys(["in","out","gnd","gndn","st"])
  val["in"]  = [i+1 for i,x in enumerate(net.getVars()) if x[1:3]=="Vi"]
  val["out"] = [i+1 for i,x in enumerate(net.getVars()) if x[1:3]=="Vo"]
  val["st"]  = [i+1 for i,x in enumerate(net.getVars()) if x[1]=="w"]
  val["gnd"] = [i+1 for i,x in enumerate(net.getVars()) if x[1:3]=="Ig"]
  val["gndn"] = [i+1 for i,x in enumerate(net.getVars()) if x[1:3]=="Vg"]

  edg = {"in":[],"out":[],"gnd":[],"st":[]}
  pos = []
  for x in net.getVars():
    if x[1:3]=="Vi":
      edg["in"].append(list(net.getVarEdge(x)[0]))
    elif x[1:3]=="Vo":
      edg["out"].append(list(net.getVarEdge(x)[0]))
      n = edg["out"][-1][0]
      pos.append([n, net.graph.nodes(data=True)[n][1]['pos']])
    elif x[1]=="w":
      edg["st"].append(list(net.getVarEdge(x)[0]))
    elif x[1:3]=="Ig":
      edg["gnd"].append(list(net.getVarEdge(x)[0]))

  pos2 = [None]*(max([x[0] for x in pos])+1)
  #import pdb; pdb.set_trace()
  for x in pos:
    pos2[x[0]] = list(x[1])
  pos = pos2
  del pos2

  if octave:
    cmd = 'ind = struct ("in",{i},"out",{o},"st",{s},"gnd",{g},"gndn",{gn});'.format( \
    i=val["in"],o=val["out"], \
    s=val["st"],g=val["gnd"],gn=val["gndn"])
    cmd = cmd + '\n' + \
    ('edges = struct ("in",reshape({i},2,{ni}),' + \
                  '"out",reshape({o},2,{no}),' + \
                  '"st",reshape({s},2,{ns}),' + \
                  '"gnd",reshape({g},2,{ng}));').format( \
    i=edg["in"],ni=len(edg["in"]), \
    o=edg["out"],no=len(edg["out"]), \
    s=edg["st"],ns=len(edg["st"]), \
    g=edg["gnd"],ng=len(edg["gnd"]))
    postxt = 'pos = [' + ';'.join(['{s}'.format(s=x) if x else '[NA,NA]'for x in pos]) + '];'
    cmd = cmd + '\n' + postxt

    return cmd, val

  return val

def trapesoidPulse (t,t_ini,t_end):

  beta=1e-3
  # Map time to [0,T] interval
  T   = t_end - t_ini;
  t_r = t - t_ini;

  if (t_r <= 0) :
    t_r = 0.0
  elif (t_r < beta):
    t_r = t_r/beta
  elif (t_r >= beta and t_r <= T - beta):
    t_r = 1.0
  elif (t_r > T - beta):
    t_r = (T - t_r)/beta
  elif (t_r > T):
    t_r = 0.0

  return t_r;

def trapesoidTrain(t,t_cycle,duty,t_max):

  if (t >= t_max):
    return 0.0

  t_end = t_cycle
  while (t_end < t):
    t_end += t_cycle

  return trapesoidPulse (t, t_end - duty*t_cycle, t_end);

def decoVar (basename,i):
  return ["$"+b+"{}".format(i) for b in basename]

def decoVars (basename,i):
  return ["$"+basename+"{}".format(i) for i in range(i)]
