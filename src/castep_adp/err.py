import sys

def no_lattice_cart():
  string = f"No lattice_cart block in cell file\n"
  string += "Reading atomic positions from a cell file requires lattice_cart block"
  print(string)
  quit()

def file_not_found(fname,ext=""):
  string = f"No file found with filename {fname}{ext}\n"
  string += "Please check the seed matches your files."
  print(string)
  quit()

def unexpected_format(line,fname):
  string = f"Unexpected format in {fname}.adp, line: {line}\n"
  string += "All lines should be of the form <parameter> : <value>"
  print(string)
  quit()

def unknown_keyword(key,fname):
  string = f"Unknown keyword {key} found in {fname}.adp"
  print(string)
  quit()

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
  quit()

class InvalidParamter(Exception):
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
        string += "true or false"
      else:
        substr = "/".join(self.optns)
        string += f"one of {substr}"
    elif self.lbound is not None:
      msg += f"Unexpected value {self.value} for {self.key}, {self.key} >= {self.lbound}"
    elif self.ubound is not None:
      msg += f"Unexpected value {self.value} for {self.key}, {self.key} <= {self.ubound}"
    else:
      n = ""
      if self.kind[0] in ["a","e","i","o","u"]: n = "n"
      msg += f"{self.key} must be a{n} {self.kind}, got {self.value} instead"
    return msg