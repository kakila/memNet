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

import sys
import os
sys.path.append('../')
import numpy as np
import networkx as nx

###########################################
#####  Topology generating functions  #####

def degreeGridGraph (n=2,m=2,d=2,l=1.0):

  HZ = np.cumsum(np.ones((n,m)),axis=1)-1
  VT = np.cumsum(np.ones((n,m)),axis=0)-1
  # Distance from all nodes to node i,j
  dm = lambda i,j: np.sqrt( (HZ-j)**2 + (VT-i)**2)

  G = nx.MultiGraph()
  G.add_nodes_from(range(n*m))
  for n in G:
    G.node[n]['pos']= (np.ravel(HZ)[n],np.ravel(VT)[n])

  nodes = G.nodes()
  while nodes:
    u = nodes.pop(np.random.randint(0,len(nodes)))
    i,j = np.unravel_index(u,(n,m))
    D = dm(i,j)
    order = np.argsort(D,axis=None)

    niter = 0
    while G.degree(u) < d and len(order)>0:
      niter = niter + 1
      # Update order to remove self and nodes with enough degree
      order = [x for x in order if (G.degree(x) < d and x != u)]
      for v in order:
        # Check if u has enough degree xor order is empty xor add an edge
        if G.degree(u) >= d:
          break
        elif len(order) == 1:
          G.add_edge(u,v)
        else:
          if np.random.random() < np.exp(-np.ravel(D)[v]/l):
              G.add_edge(u,v)

      if niter > np.log(1e3*(d*n*m)):
       print "Maximum iterations reached!"
       break

  return G

def waxmanGraph (n,alpha=1,beta=0.5,gamma=0.5):

  connected = False
  maxIter = 100
  j=0
  while not connected and j < maxIter:
    g = nx.waxman_graph (n,alpha=alpha,beta=beta)
    connected = nx.is_connected (g)
    j += 1

  if j >= maxIter:
    print "Generated {} graphs and none was connected.".format(j)
    raise ValueError

  nelm = g.number_of_edges()

  ## Swap grid randomly
  swap = np.max([np.round(gamma*nelm),1])
  ind = np.random.permutation(nelm)[0:swap]
  ngrid = g.edges()
  for i in ind:
    ngrid[i] = ngrid[i][::-1]

  # Find nodes that are the most appart
  pos = np.array(list(nx.get_node_attributes(g,'pos').values()))
  x = pos[:,0] - pos[:,0][np.newaxis].T
  y = pos[:,1] - pos[:,1][np.newaxis].T
  d = x**2 + y**2
  nodes = np.unravel_index(np.argmin(-d),d.shape)

  return ngrid,nodes

def gridGraph(rows=2,cols=2, kind="forward", swap=[]):
  from scipy.weave import inline
  n = rows*cols
  edges = np.zeros([n,n],dtype=np.int)
  code = r'''
         for (unsigned int j=0; j<n; ++j){
           for (unsigned int i=0; i<n; ++i){
             if ((j==i+1 && (i+1)%rows !=0) || j==i+rows) {
               edges[j+n*i] = 1;
             }
           }
          }
         '''
  inline(code,['n','rows','edges'],verbose=1)

  ind     = edges.nonzero()
  n_edges = len(ind[0]) # Number of edges

  if kind == "forward" and not swap:
    pass

  elif kind == "backward" or (swap and swap >= n_edges):
    edges = edges.T

  elif kind == "random" or swap:

    if swap:
      p = np.random.permutation(n_edges)[0:swap]
      for i in p:
        edges[ind[1][i],ind[0][i]]= 1
        edges[ind[0][i],ind[1][i]]= 0

    else:
      p = np.random.randint(0,2,n_edges)
      for i,v in enumerate(p):
        if v:
          edges[ind[1][i],ind[0][i]]= 1
          edges[ind[0][i],ind[1][i]]= 0
  else:
    raise ValueError("Unkown kind of grid:{}".format(kind))

  ind = edges.nonzero()
  edges = zip(ind[0],ind[1])
  return edges,n_edges

def generateGrid (n,m=None, memR=None, n_in=None,d=3,**kwargs):
  '''
  Generates the topology of the network.
  '''
  from circuit_elem import Resistor as res
  from circuit_elem import HP_TiO2 as hp

  if not memR:
    memR = hp
  if not m:
    m = n
  if not n_in:
    n_in = n

  ### Define Ground, and sides of the network
  Side = dict(E=[],N=[],W=[],S=[])
  Side["S"] = range(m)
  Side["W"] = range(0,n*m,m)
  Side["E"] = range(m-1,n*m,m)
  Side["N"] = range(m*(n-1),n*m,1)
  GND = -1;
  GND_n = Side["E"]
  IN  = range(n*m,n*m+n_in)

  # Generate grid
  G = degreeGridGraph(n=n,m=m,d=d,**kwargs)

  ## Erase connections between nodes in East side (grounded)
  notWorE = [i for i in xrange(n*m) if i not in Side["W"]+Side["E"]]
  if len(notWorE) == 0:
    notWorE = Side["E"][:]
  N = len(notWorE)

  to_add = []
  to_rem = []
  for i in Side["E"]:
    for j in [ii for ii in G.neighbors(i) if ii in Side["E"]]:
      to_rem.append((i,j))
      # Add a random edge to keep degree
      k = notWorE[np.random.randint(0,N)]
      to_add.append((i,k))

  G.remove_edges_from(to_rem)
  G.add_edges_from(to_add)

  ## Make sure we have a single connected component
  ## without adding edges between East nodes (grounded)
  if nx.number_connected_components(G) > 1:
    to_add=[]
    for x in nx.connected_components(G):
      xc = x
      if len(xc)>1:
          xc = [i for i in xc if i not in Side["E"]]
      to_add.append(xc[np.random.randint(0,len(xc))])
    G.add_path(to_add)

  ## Connect IN nodes with West side respecting degree
  for i in IN:
    to_add = set ()

    if len (Side["W"]) <= d: # go one by one
      for k in Side["W"]:
        to_add.add((i,k))
    else: # add without repetition
      while len(to_add) < d:
        k = Side["W"][np.random.randint(0,len(Side["W"]))]
        to_add.add((i,k))

    G.add_edges_from(list(to_add))

  for y,i in enumerate(IN):
    G.add_node(i,pos=(-1,y))

  '''
  ## Erase connections between nodes in West side
  notWorE = [i for i in xrange(n*m) if i not in Side["W"]+Side["E"]]
  if len(notWorE) == 0:
    notWorE = Side["E"][:]
  N = len(notWorE)

  to_add = []
  to_rem = []
  for i in Side["W"]:
    for j in [ii for ii in G.neighbors(i) if ii in Side["W"]]:
      to_rem.append((i,j))
      # Add a random edge to keep degree
      k = notWorE[np.random.randint(0,N)]
      to_add.append((i,k))

  G.remove_edges_from(to_rem)
  G.add_edges_from(to_add)
  '''

  grid = G.edges()
  nelm = len(grid)

  ## Swap grid randomly
  swap = np.max([np.round(0.5*nelm),1])
  ind  = np.random.randint(0,nelm,swap)
  for i in ind:
    grid[i] = grid[i][::-1]

  ## Generate input edges
  inp  = [(GND,x) for x in IN]

  ### Generate readouts edges
  OUT = [x for x in xrange (n*m) if x not in IN+GND_n]
  outp = [(x,GND) for x in OUT]
  nout = len(outp)
  ## input readout nodes
  outi = [(x,GND) for x in IN]

  ### Generate grounded edges
  gndp = [(x,GND) for x in GND_n]

  ### Store list and Type of each edge
  edges=dict().fromkeys(["net","input_readout","output","ground"])
  edges["net"]           = (grid, [memR()]*len(grid))
  edges["output"]        = (outp, [res(R="Ro",I="Io",V="Vo")]*len(outp))
  edges["input_readout"] = (outi, [res(R="Ro",I="Ii",V="Vi")]*len(outi))
  edges["ground"]        = (gndp, [res(R="Rg",I="Ig",V="Vg")]*len(gndp))

  return edges,inp,G

def simpleGrid (n=4,m=4,d=1, p=0.5, memR=None):

  from circuit_elem import Resistor as res
  from circuit_elem import HP_TiO2 as hp

  if not memR:
    memR = hp

  tmp = nx.grid_2d_graph(n, m)
  G = nx.DiGraph ()
  G.add_edges_from(tmp.edges())
  del tmp

  NODE_id = [x[0]+m*x[1] for x in G]
  ### Define Ground, and sides of the network
  Side = dict(E=set(),N=set(),W=set(),S=set())
  for x in G:
    k = x[0]+m*x[1]
    if x[0]==n-1:
      Side["S"].add(k)
    if x[0]==0:
      Side["N"].add(k)
    if x[1]==0:
      Side["W"].add(k)
    if x[1]==m-1:
      Side["E"].add(k)

  GND = -1;
  GND_n = Side["E"].intersection(Side["S"])
  IN    = Side["W"].intersection(Side["N"])
  OUT   = set(NODE_id).difference(IN.union(GND_n))
  NODE_id.append(GND)

  E = G.edges()
  N = G.nodes()
  pos = [{'pos':n} for n in N]
  pos.append({'pos':(-1,-1)})
  grid = [(NODE_id[N.index(e[0])],NODE_id[N.index(e[1])]) for e in E]
  nelm = len(grid)

  ## Swap grid randomly
  swap = np.max([np.round(p*nelm),1])
  ind  = np.random.randint(0,nelm,swap)
  for i in ind:
     grid[i] = grid[i][::-1]

  ## Increase degree
  s_grid = []
  for i in xrange(d):
    s_grid.append(grid)
  grid = [item for subgrid in s_grid for item in subgrid]
  nelm = len(grid)

  ## Generate input edges
  inp  = [(GND,x) for x in IN]

  ### Generate readouts edges
  outp = [(x,GND) for x in OUT]
  nout = len(outp)
  ## input readout nodes
  outi = [(x,GND) for x in IN]

  ### Generate grounded edges
  gndp = [(x,GND) for x in GND_n]

  G = nx.MultiDiGraph()
  G.add_edges_from(grid+inp+outp+outi+gndp)
  G.add_nodes_from(zip(NODE_id,pos))

  ### Store list and Type of each edge
  edges=dict().fromkeys(["net","input_readout","output","ground"])
  edges["net"]           = (grid, [memR()]*len(grid))
  edges["output"]        = (outp, [res(R="Ro",I="Io",V="Vo")]*len(outp))
  edges["input_readout"] = (outi, [res(R="Ro",I="Ii",V="Vi")]*len(outi))
  edges["ground"]        = (gndp, [res(R="Rg",I="Ig",V="Vg")]*len(gndp))

  return edges,inp,G

def symmetricGrid (n=4,m=4,d=1, p=0.5, memR=None):

  tmp = nx.grid_2d_graph(n, m)
  G = nx.DiGraph()
  G.add_edges_from(tmp.edges())
  del tmp

  NODE_id = [x[0]+m*x[1] for x in G]
  ### Define Ground, and sides of the network
  Side = dict(E=set(),N=set(),W=set(),S=set())
  for x in G:
    k = x[0]+m*x[1]
    if x[0]==n-1:
      Side["S"].add(k)
    if x[0]==0:
      Side["N"].add(k)
    if x[1]==0:
      Side["W"].add(k)
    if x[1]==m-1:
      Side["E"].add(k)

  GND = -1;
  GND_n = Side["E"].intersection(Side["S"])
  IN    = Side["W"].intersection(Side["N"])
  OUT   = set(NODE_id).difference(IN.union(GND_n))
  NODE_id.append(GND)

  E = G.edges()
  N = G.nodes()
  pos = [{'pos':n} for n in N]
  pos.append({'pos':(-1,-1)})
  grid = [(NODE_id[N.index(e[0])],NODE_id[N.index(e[1])]) for e in E]
  nelm = len(grid)

  ## Swap grid randomly
  swap = np.max([np.round(p*nelm),1])
  ind  = np.random.randint(0,nelm,swap)
  for i in ind:
     grid[i] = grid[i][::-1]

  ## Symmetrize grid
  s_grid = []
  for i in xrange(d):
    s_grid.append(grid)
    s_grid.append([x[::-1] for x in grid])
  grid = [item for subgrid in s_grid for item in subgrid]
  nelm = len(grid)

  ## Generate input edges
  inp  = [(GND,x) for x in IN]

  ### Generate readouts edges
  outp = [(x,GND) for x in OUT]
  nout = len(outp)
  ## input readout nodes
  outi = [(x,GND) for x in IN]

  ### Generate grounded edges
  gndp = [(x,GND) for x in GND_n]

  G = nx.MultiDiGraph()
  G.add_edges_from(grid+inp+outp+outi+gndp)
  G.add_nodes_from(zip(NODE_id,pos))

  return grid, outp, outi, gndp,inp,G
