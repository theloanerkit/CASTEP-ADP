import sys


def no_lattice_cart():
  string = "No lattice_cart block in cell file\n"
  string += "Reading atomic positions from a cell file requires lattice_cart block"
  print(string)
  sys.exit()

def file_not_found(fname,ext=""):
  string = f"No file found with filename {fname}{ext}\n"
  string += "Please check the seed matches your files."
  print(string)
  sys.exit()

def unexpected_format(line,fname):
  string = f"Unexpected format in {fname}.adp, line: {line}\n"
  string += "All lines should be of the form <parameter> : <value>"
  print(string)
  sys.exit()

def unknown_keyword(key,fname):
  string = f"Unknown keyword {key} found in {fname}.adp"
  print(string)
  sys.exit()

def invalid_parameter(key,value,optns=None,kind=None,lbound=None,ubound=None):
  string = ""
  if kind is None and lbound is None and ubound is None:
    string += f"Unexpected value {value} for {key}, expected "
    if optns is None:
      string += "true or false"
    else:
      substring = "/".join(optns)
      string += f"one of {substring}"
  elif lbound is not None:
    string += f"Unexpected value {value} for {key}, {key} >= {lbound}"
  elif ubound is not None:
    string += f"Unexpected value {value} for {key}, {key} <= {ubound}"
  else:
    n = ""
    if kind[0] in ["a","e","i","o","u"]: n = "n"
    string += f"{key} must be a{n} {kind}, got {value} instead"
  print(string)
  sys.exit()

class UnitError(Exception):
  def __init__(self,kind,unit,dim=None,loc=None):
    self.kind = kind
    self.loc = loc
    self.unit = unit
    self.dim = dim
    self.__suppress_context__=True
    sys.tracebacklimit=0

  def __str__(self):
    msg = ""
    if self.kind == "unknown unit":
      msg += f"Unknown unit {self.unit} in {self.loc}"
    if self.kind == "dimensionality":
      msg += f"Inconsistent dimensionality in units: expected {self.dim}, "
      msg += f"got {self.unit} with dimensionality {self.unit.dimensionality}"
    return msg

class IncompatibleCellError(Exception):
  def __init__(self,kind):
    self.kind = kind
    self.__suppress_context__=True
    sys.tracebacklimit=0

  def __str__(self):
    msg = ""
    if self.kind == "no lattice cart":
      msg += "No lattice_cart block in cell file\n"
      msg += "Reading atomic positions from a cell file requires lattice_cart block"
    if self.kind == "no positions block":
      msg += "No positions block in cell file"
    return msg

class NoLatticeCartError(IncompatibleCellError):
  def __init__(self):
    super().__init__(kind="no lattice cart")

class NoPositionsBlockError(IncompatibleCellError):
  def __init__(self):
    super().__init__(kind="no positions block")

class InvalidParamterError(Exception):
  def __init__(self,key,value,optns=None,kind=None,lbound=None,ubound=None):
    self.key = key
    self.value = value
    self.optns = optns
    self.kind = kind
    self.lbound = lbound
    self.ubound = ubound
    self.__suppress_context__=True
    sys.tracebacklimit=0

  def __str__(self):
    msg = ""
    if self.kind is None and self.lbound is None and self.ubound is None:
      msg += f"Unexpected value {self.value} for {self.key}, expected "
      if self.optns is None:
        msg += "true or false"
      else:
        substr = "/".join(self.optns)
        msg += f"one of {substr}"
    elif self.lbound is not None:
      msg += f"Unexpected value {self.value} for {self.key}, {self.key} >= {self.lbound}"
    elif self.ubound is not None:
      msg += f"Unexpected value {self.value} for {self.key}, {self.key} <= {self.ubound}"
    else:
      n = ""
      if self.kind[0] in ["a","e","i","o","u"]: n = "n"
      msg += f"{self.key} must be a{n} {self.kind}, got {self.value} instead"
    return msg

class ParseMDError(Exception):
  def __init__(self,kind,equ_timesteps=None):
    self.kind = kind
    self.equ_timesteps = equ_timesteps
    self.__suppress_context__=True
    sys.tracebacklimit=0

  def __str__(self):
    msg = "Error parsing MD file: "
    if self.kind == "no timesteps":
      msg += "no timesteps found\n"
      if self.equ_timesteps is None:
        msg += "  This is likely due to either an empty MD file, "
        msg += "or the number of equilibration timesteps set too high"
      elif self.equ_timesteps == 0:
        msg += "  This is likely due to an empty MD file"
      else:
        msg += "  This is likely due to the number of equilibration timesteps being set too high"
        msg += f" (equilibration timesteps set to {self.equ_timesteps})"
    return msg

class NoTimestepsError(ParseMDError):
  def __init__(self,equ_timesteps):
    super().__init__(kind="no timesteps",equ_timesteps=equ_timesteps)
