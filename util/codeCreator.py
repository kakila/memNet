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


from time import strftime,localtime

from io import open
from string import Template


class codeCreator():

  def __init__(self,name,hpath="./include/",cpath="./src/"):

    self.license = '''/*
 * Copyright (C) {} - Juan Pablo Carbajal
 *
 * This progrm is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

 // This code was automatically generated from codeCreator.py
'''.format(strftime("%Y",localtime())) + "\n"

    self.name   = name
    self.path   = dict(h=hpath+name+".h",c=cpath+name+".c")
    self.header = Template('''// {name}.h
#ifndef _{nameC}_H
#define _{nameC}_H

$includes
// Variables
$variables
// Function prototypes
$functions

#endif // _{nameC}_H
'''.format(nameC=name.upper(),name=name))

    self.code = Template('''// {}.c
$includes
$variables
$functions
'''.format(name))

    self.strcode  = ""
    self.strheader= ""
    self.givens   = dict()

  def writeout (self,**kwargs):
    if kwargs.has_key("verbose"):
      v = kwargs.pop("verbose")

    if kwargs:
      self.set(**kwargs)

    if not self.updateCode(verbose=v,**kwargs):
      raise NameError("Could not write the code.")

    if not self.updateHead(verbose=v,**kwargs):
      raise NameError("Could not write the header.")

    self.writeHeader()
    self.writeCode()

  def updateHead(self,**kwargs):
    if kwargs.has_key("verbose"):
      v = kwargs.pop("verbose")

    if all (self.check_givens(verbose=v)):
      self.strheader= self.header.safe_substitute (self.givens)
      return True
    else:
      return False

  def updateCode(self,**kwargs):
    if kwargs.has_key("verbose"):
      v = kwargs.pop("verbose")

    if all (self.check_givens(verbose=v)):
      self.strcode= self.code.substitute (self.givens)
      return True
    else:
      return False

  def writeHeader (self):
    hfile = open (self.path["h"], 'w')
    hfile.write(unicode(self.strheader))
    hfile.close()

  def writeCode (self):
    cfile = open (self.path["c"], 'w')
    cfile.write(unicode(self.strcode))
    cfile.close()

  def check_givens (self,verbose=True):
    val = [False]*len(self.givens)
    for k,v in self.givens.iteritems():
      if type(v) is str:
        if v[0]=="$":
          if verbose:
            print k + " is not defined"
            continue

      if verbose:
        print k + ": {}".format(v)
      val[self.givens.keys().index(k)] = True

    return val

  def set(self,**kwargs):
    for k,v in kwargs.iteritems():
      if self.givens.has_key(k):
        method = getattr(self,"format_"+k, lambda x: x)
        self.givens[k] = method(v)
      else:
        raise ValueError ("The main file doesn't have a value named "+k)

class RHS (codeCreator):

  def __init__(self,name,**kwargs):
    codeCreator.__init__(self,name,**kwargs)

    self.strheader = self.header.substitute(\
          variables="extern realtype p_rhs[];\n" +\
                    "extern unsigned int Nprhs;\n",\
          functions= \
  "extern int rhs (realtype t, N_Vector x,N_Vector dx, N_Vector r, void *u);\n" +\
  "extern void IDfunc (N_Vector id);\n", \
          includes="#include <sundials/sundials_types.h>\n"
                                       )

    self.givens = dict (body = "$body",params = "$params", \
                        inputsfile = "$inputsfile")

    self.code = Template (\
          self.code.safe_substitute (\
          variables="realtype p_rhs[] = {" + self.givens["params"] + "};\n",
          functions= self.givens["body"],
          includes= "#include <math.h>\n" +\
                    '#include "helper.h"\n' +\
                    "#include <nvector/nvector_serial.h>\n"+\
                    "#include <sundials/sundials_types.h>\n"+\
                    '#include "' + self.givens["inputsfile"] +'"\n'
                                    )
                        )
    self.strcode=""

  def format_body (self,kw):
    """ Creates the string with the C code for the body of the right
    hand side.

    Input
    -----
    kw: A dict with the fileds "vrs","rhs", "params" and "input_name".
        "vars": a list of string with the name of the variables in the problem.
        "rhs": A dict with one key per variable in vars. the value of each key
               is a list with two elements, the first element is the string
               of the equation and the secind element is a boolean, that if
               True then the equation describes a differential equation.
        "params": a dict with the name of the parameters in its keys,
                  each one containing the value of the parameter.
        "input_name": The name of the input functions. If omitted it is asusmed that
                 there are no inputs.
    """
    var,rhs,par = kw["vars"],kw["rhs"],kw["params"]

    # Id of varibale sindicating differential or algebraic
    tmp = list ()
    for i,v in enumerate(var):
      if rhs.has_key(v):
        tmp.append("1.0" if rhs[v][1] else "0.0")

    strRhs_body ="void IDfunc (N_Vector id)\n" +\
                 "{\n\n"+\
                 "  realtype ID[]  = {" + ",".join(tmp) + "};\n"+\
                 "  realtype *data;\n"+\
                 "  unsigned int i;\n"+\
                 "  data = NV_DATA_S(id);\n"+\
                 "  for (i=0; i<{};++i) data[i] = ID[i];\n".format(len(tmp)) +\
                 "};\n\n"

    strRhs_body +=\
    "int rhs (realtype t, N_Vector x,N_Vector dx, N_Vector r, void *u)\n"+\
    "{\n\n" +\
    "  realtype *xval, *dxval, *rval;\n" +\
    "  xval  = NV_DATA_S(x);\n" + \
    "  dxval = NV_DATA_S(dx);\n" + \
    "  rval  = NV_DATA_S(r);\n\n" + \
    "  // Update inputs\n" + \
    "  inputs (t, x, dx);\n\n" + \
    "  //System equations\n"
    for i,v in enumerate(var):
      if rhs.has_key(v):
        tmp = Template (rhs[v][0].replace(" ", ""))
        subs_dict = dict([(u[1:],' xval[{}] '.format(j)) for j,u in enumerate(var)])

        strRhs_body += '  rval[{}] = '.format(i) + tmp.safe_substitute(subs_dict)
        if rhs[v][1]:
          strRhs_body += '- dxval[{}]'.format(i)

      strRhs_body +=';\n'
      subs_dict = dict(\
        [ ( x[1:],"p_rhs[{}]".format(i) ) for i,x in enumerate(par.keys())]
                      )
      strRhs_body = Template(strRhs_body).safe_substitute(subs_dict)

      if kw.has_key("input_name"):
        subs_dict = dict(\
          [ ( x[1:],"in[{}]".format(i) ) for i,x in enumerate(kw["input_name"])]
                        )
        strRhs_body = Template(strRhs_body).substitute(subs_dict)

    strRhs_body +="  return(0);\n};"

    return strRhs_body

  def format_params (self,par):
    strParam_values = "\n";
    for k,p in par.iteritems():
      strParam_values += " "*6 + "{:e}, {}\n".format(p,"//"+k)
    return strParam_values

  def format_inputsfile(self,f):
    ind = f.rfind("/")
    return f[ind+1:] if ind else f

  def updateHead(self,**kwargs):
    return True

class INPUT(codeCreator):

  def __init__(self,name,**kwargs):
    codeCreator.__init__(self,name,**kwargs)

    self.strheader = self.header.substitute(\
          variables="extern realtype p_in[];\n" +\
                    "extern unsigned int Npin;\n"+\
                    "extern realtype in[];\n" +\
                    "extern unsigned int Nin;\n",\
          functions= \
  "extern void inputs (realtype t, N_Vector x, N_Vector dx);\n",
          includes=
                   "#include <sundials/sundials_types.h>\n"
                                       )

    self.givens = dict(body="$body",params="$params", \
                       Ninputs = "$Ninputs")

    self.code = Template (\
          self.code.safe_substitute (\
          variables="realtype p_in[] = {" + self.givens["params"] + "};\n"+\
                    "unsigned int Npin = sizeof(p_in)/sizeof(p_in[0]);\n"+\
                    "unsigned int Nin = " + self.givens["Ninputs"] + ";\n"
                    "realtype in[" + self.givens["Ninputs"] + "];\n",
          functions= \
         "void inputs (realtype t, N_Vector x, N_Vector dx)\n"+\
         "{\n\n" +\
         "  //realtype *xval, *dxval;\n" +\
         "  //xval  = NV_DATA_S(x);\n" + \
         "  //dxval = NV_DATA_S(dx);\n" + \
          self.givens["body"] +"\n};",
          includes = "#include <math.h>\n"+\
                     "#include <nvector/nvector_serial.h>\n"+\
                     "#include <sundials/sundials_types.h>\n"+\
                     "// Helper functions\n"+\
                     '#include "helper.h"\n'
                                    )
                        )

  def format_body(self,kw):
    var,inputs,par = kw["vars"],kw["inputs"],kw["params"]
    strInput_body =\
    "  //Input functions\n"
    for k,v in inputs.iteritems():
      tmp = Template (v.replace(" ", ""))
      subs_dict = dict([(u[1:],' xval[{}] '.format(j)) for j,u in enumerate(var)])

      strInput_body += '  in[{}] = '.format(inputs.keys().index(k)) +\
                         tmp.safe_substitute(subs_dict)

      strInput_body +=';\n'

    subs_dict = dict(\
      [ ( x[1:],"p_in[{}]".format(i) ) for i,x in enumerate(par.keys())]
                    )
    strInput_body = Template(strInput_body).substitute(subs_dict)

    return strInput_body

  def format_params (self,par):
    strParam_values = "\n";
    for k,p in par.iteritems():
      strParam_values += " "*6 + "{:e}, {}\n".format(p,"//"+k)

    return strParam_values

  def updateHead(self,**kwargs):
    return True

class ICS(codeCreator):

  def __init__(self,name,**kwargs):
    codeCreator.__init__(self,name,**kwargs)

    self.strheader = self.header.substitute(\
          variables="",\
          functions= \
  "extern void ICSfunc (realtype *xval, realtype *dxval);\n",
          includes=
                   "#include <sundials/sundials_types.h>\n"
                                       )

    self.givens = dict(body="$body",params="$params")

    self.code = Template (\
          self.code.safe_substitute (\
          variables="",
          functions= \
         "void ICSfunc (realtype *xval, realtype *dxval)\n"+\
         "{\n\n" +\
          self.givens["body"] +"\n};",
          includes = "#include <math.h>\n"+\
                     "#include <nvector/nvector_serial.h>\n"+\
                     "#include <sundials/sundials_types.h>\n"+\
                     "// Helper functions\n"+\
                     '#include "helper.h"\n'
                                    )
                        )

  def updateHead(self,**kwargs):
    return True

  def format_body (self,kw):
    var,ics,par = kw["vars"],kw["ics"],kw["params"]
    str_body =\
    "  //Inital conditions\n" +\
'\n'.join(['  xval[{}]  = RCONST({:e});'.format(var.index(val),ics[val][0]) \
                                    for i,val in enumerate(ics.keys())]) + \
'\n' + \
'\n'.join(['  dxval[{}] = RCONST({:e});'.format(var.index(val),ics[val][1]) \
                                    for i,val in enumerate(ics.keys())])

    return str_body

class MAIN(codeCreator):

  def __init__(self,name="main",**kwargs):
    codeCreator.__init__(self,name,**kwargs)

    self.givens = dict(Neq="$Neq",rhsfile = "$rhsfile", \
                      tolerances = "$tolerances",outfile ="$outfile", \
                      wmode = "$wmode", tspan= "$tspan", period="$period")

    self.header = Template(self.header.safe_substitute(\
          variables="#define NEQ "+self.givens["Neq"] +"\n",\
          functions="",
          includes= \
          "#include <math.h>\n"+\
          "#include <ida/ida.h>\n"+\
          "#include <ida/ida_dense.h>\n"+\
          "#include <ida/ida_lapack.h>\n"+\
          "#include <nvector/nvector_serial.h>\n"+\
          "#include <sundials/sundials_math.h>\n"+\
          "#include <sundials/sundials_types.h>\n"+\
          "// Helper functions\n"+\
          '#include "helper.h"\n'+\
          "// RHS and INPUTS\n"+\
          '#include "' + self.givens["rhsfile"] +'"\n'
                                       ))

    self.code = Template (\
          self.code.safe_substitute (\
          variables="",
          functions= \
         "int main ()\n"+\
         "{\n\n" +\
         '''  void *mem;
  N_Vector yy, yp, id, avtol;
  realtype rtol, *yval, *ypval, *atval;
  realtype tstart, tend, tstep, tout, tret, tstop;
  int retval, i, nstop=1;
  realtype defatol;
  FILE *fout;

  // Initialize random number generator
  if ( getenv ("SEED") != NULL )
    srand (atoi (getenv ("SEED")));
  else
    srand (time(NULL));

  // Open files to write output
  fout = fopen("{fname}","{mode}");'''.format(fname=self.givens["outfile"],\
                                              mode=self.givens["wmode"]) +\
  '''
  if (fout == NULL) {''' +\
  '''printf("Could not open {fname} for writing.\\n");'''.format(fname=self.givens["outfile"],\
                                                                 mode=self.givens["wmode"]) +\
"""
  return (1);
  }

  mem   = NULL;
  yy    = NULL;
  yp    = NULL;
  avtol = NULL;
  yval  = NULL;
  ypval = NULL;
  atval = NULL;
  id    = NULL;

  // Allocate N-vectors.
  yy = N_VNew_Serial(NEQ);
  if (check_flag ((void *)yy, "N_VNew_Serial", 0)) return (1);

  yp = N_VNew_Serial(NEQ);
  if (check_flag ((void *)yp, "N_VNew_Serial", 0)) return (1);

  id  = N_VNew_Serial(NEQ);
  if(check_flag((void *)id, "N_VNew_Serial", 0)) return(1);

  avtol = N_VNew_Serial(NEQ);
  if (check_flag ((void *)avtol, "N_VNew_Serial", 0)) return (1);

  // Create and initialize  y, dydt, id.
  yval  = NV_DATA_S(yy);
  ypval = NV_DATA_S(yp);
  ICSfunc (yval, ypval);
  IDfunc (id);

  // Tolerances
""" + self.givens["tolerances"] +\
"""
  // Integration limits
""" + self.givens["tspan"] + \
"""
  /* Call IDACreate and IDAInit to initialize IDA memory */
  mem = IDACreate ();
  if (check_flag ((void *)mem, "IDACreate", 0)) return (1);

  retval = IDASetId(mem, id);
  if(check_flag(&retval, "IDASetId", 1)) return(1);

  retval = IDAInit (mem, rhs, tstart, yy, yp);
  if (check_flag (&retval, "IDAInit", 1)) return (1);

  /* Call IDASVtolerances to set tolerances */
  retval = IDASVtolerances (mem, rtol, avtol);
  if (check_flag(&retval, "IDASVtolerances", 1)) return (1);

  /* Free avtol */
  N_VDestroy_Serial (avtol);

  // Call IDADense and set up the linear solver.
  //  retval = IDADense (mem, NEQ);
  //  if (check_flag(&retval, "IDADense", 1)) return (1);

  // Call IDALapackDense and set up the lapack linear solver.
  retval = IDALapackDense (mem, NEQ);
  if (check_flag(&retval, "IDALapackDense", 1)) return (1);

  // Set up initial time step.
  retval = IDASetInitStep(mem, 1e-6);
  if (check_flag(&retval, "IDASetInitStep", 1)) return (1);

  // Set up maximum time step.
  retval = IDASetMaxStep(mem, 1e-1);
  if (check_flag(&retval, "IDASetMaxStep", 1)) return (1);

  // Set next stop time.
""" +\
"  tstop = {t};".format(t=self.givens["period"]) + \
"""
  retval = IDASetStopTime(mem, tstop);
  if (check_flag(&retval, "IDASetStopTime", 1)) return (1);

  // Main loop
  tout = tstart + tstep;

  // Correct Initial conditions
  retval = IDACalcIC(mem, IDA_YA_YDP_INIT, tout);
  if (check_flag(&retval, "IDACalcIC", 1)) return (1);

  while (TRUE) {

    Print2File(fout, NEQ, tret, yy);

    retval = IDASolve (mem, tout, &tret, yy, yp, IDA_NORMAL);
    if (check_flag (&retval, "IDASolve", 1)) return (1);
          
    if (retval == IDA_TSTOP_RETURN)
    {
      tout = tret;
      nstop += 1;
      // Set next stop time.
      i = IDASetStopTime(mem, nstop*tstop);
      if (check_flag(&i, "IDASetStopTime", 1)) return (1);
    }
   
    tout += tstep;

    if (tout > tend) break;
  }
  Print2File(fout, NEQ, tret, yy);
  // Print velocities, used for restart
  Print2File(fout, NEQ, tret, yp);

  PrintFinalStats (mem);

  // Free memory
  fclose(fout);
  IDAFree (&mem);
  N_VDestroy_Serial (yy);
  N_VDestroy_Serial (yp);
  N_VDestroy_Serial(id);
  return (0);

};
""",
          includes = '#include "{}"\n'.format(self.name+".h")
                                    )
                        )

  def format_tspan(self, t):
    if type(t) is not list:
      raise ValueError ("Time span must be a list fo real values.")

    strT = '  tstart = RCONST({:e});\n'.format(t[0]) + \
           '  tend   = RCONST({:e});\n'.format(t[1])
    dt = (t[1]-t[0])/10.0
    if len(t) == 3:
      dt = t[2]

    strT += '  tstep  = RCONST({:e});\n'.format(dt)
    return strT

  def format_tolerances(self,kw):
    v = kw["vars"]
    strTol = \
    '  rtol = RCONST({:e});\n'.format(1.0e-4 if not kw.has_key("rTol") else kw["rTol"])

    strTol += '  atval = NV_DATA_S(avtol);\n'
    strTol += '''
  defatol = RCONST(1.0e-6);
  for(i=0;i<{};++i)
    atval[i] = defatol;
'''.format(len(v))

    if kw.has_key("aTol"):
      strTol += \
      "".join(['  atval[{}] = RCONST({:e});\n'.format(v.index(val),kw["aTol"][val]) \
                                      for i,val in enumerate(kw["aTol"].keys())])

    return strTol

  def format_rhsfile(self,f):
    ind = f.rfind("/")
    return f[ind+1:] if ind else f

class Makefile ():

  def __init__(self,**kwargs):
    """ Creates a makefile with the names of the files to be compiled.

    Inuts must be all strings.
    """
    if map(type,kwargs.values()) != [str]*len(kwargs):
      raise ValueError("Makefile accepts only str arguments.")

    self.name = "Makefile_"+kwargs["main"]
    self.body = Template ("""
DIR = $cpath
CFLAGS= -Wall -Werror -fpic -c
LFLAGS= -shared
MFLAGS= -l$rhs -l$ics -l$inputs -lsundials_ida -lsundials_nvecserial -lm -lblas -llapack
MFILE= --file=$makefile

$(DIR)%.o: $(DIR)%.c
	gcc -I$hpath $(CFLAGS) $< -o $@

$main: $(DIR)helper.o $(DIR)$main.c
	gcc -I$hpath $(DIR)$main.c $(DIR)helper.o -L. -Wl,-rpath=. $(MFLAGS) -o $@

lib$inputs.so: $(DIR)$inputs.o
	gcc ${LFLAGS} $< -o $@

lib$rhs.so: $(DIR)$rhs.o
	gcc ${LFLAGS} $< -o $@

lib$ics.so: $(DIR)$ics.o
	gcc ${LFLAGS} $< -o $@

.PHONY: clean all libs

clean:
	 rm -f $(DIR)$rhs.o $(DIR)$inputs.o $(DIR)$ics.o

dist-clean:
	 make ${MFILE} clean
	 rm -f lib$inputs.so lib$rhs.so lib$ics.so $main

libs:
	 make ${MFILE} lib$inputs.so
	 make ${MFILE} lib$rhs.so
	 make ${MFILE} lib$ics.so

all:
	 make ${MFILE} libs
	 make ${MFILE} $main
""").safe_substitute(inputs=kwargs["inputs"],
                     rhs=kwargs["rhs"],
                     ics=kwargs["ics"],
                     main=kwargs["main"],
                     makefile=self.name,
                     hpath=kwargs["hpath"],
                     cpath=kwargs["cpath"])

  def writeout (self):
    mfile = open (self.name, 'w')
    mfile.write(unicode(self.body))
    mfile.close()
