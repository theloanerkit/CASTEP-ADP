from . import err
from . import adp_constants
from . import obj

def check_int(obj,key,ubound=None,lbound=None):
  # check if parameter is an integer
  # optionally check if parameter value is within set bounds
  if not obj.valid:   # something else has gone wrong
    return
  
  # check if integer
  if obj.settings[key][0] in ("-","+"):
    is_int = obj.settings[key][1:].isdigit()
  else:
    is_int = obj.settings[key].isdigit()
  
  if not is_int:
    raise err.InvalidParamter(key,obj.settings[key],kind="integer")
  else:
    obj.settings[key] = int(obj.settings[key])

  if lbound is not None and obj.settings[key] < lbound:   # check lower bound
      obj.valid = False
      raise err.InvalidParamter(key,obj.settings[key],lbound=lbound)
  
  if ubound is not None and obj.settings[key] > ubound:   # check upper bound
      obj.valid = False
      raise err.InvalidParamter(key,obj.settings[key],ubound=ubound)
  
def check_float(obj,key,ubound=None,lbound=None):
  # check if parameter is a float
  # optionally check if parameter is within set bounds
  if not obj.valid:    # something else has gone wrong
     return
  
  # check if float
  if obj.settings[key][0] in ("-","+"):
    is_float = obj.settings[key][1:].replace(".","",1).isdigit()
  elif obj.settings[key][0] == ".":
    is_float = obj.settings[key][1:].isdigit()
  else:
    is_float = obj.settings[key].replace(".","",1).isdigit()

  if not is_float:
    raise err.InvalidParamter(key,obj.settings[key],kind="float")
  else:
    obj.settings[key] = float(obj.settings[key])

  if lbound is not None and obj.settings[key] < lbound:   # check lower bound
      obj.valid = False
      raise err.InvalidParamter(key,obj.settings[key],lbound=lbound)
  
  if ubound is not None and obj.settings[key] > ubound:   # check upper bound
      obj.valid = False
      raise err.InvalidParamter(key,obj.settings[key],ubound=ubound)

def check_boolean(obj,key):
  # check if parameter is a boolean
  if obj.valid and obj.settings[key] in ["t","true",True]:
    obj.settings[key] = True
  elif obj.valid and obj.settings[key] in ["f","false",False]:
    obj.settings[key] = False
  else:
    obj.valid = False
    raise err.InvalidParamter(key,obj.settings[key])
    
def check_keywords(obj,key,keywords,map=None):
  # check if parameter is within a set of keywords
  if obj.valid and obj.settings[key] not in keywords:
    obj.valid = False
    raise err.InvalidParamter(key,obj.settings[key],optns=keywords)
  if map is not None:
    if obj.settings[key] in map.keys():
      obj.settings[key] = map[obj.settings[key]]

def check_string_arr(obj,key,keywords=None):
    if obj.valid and keywords is None:
        obj.settings[key] = obj.settings[key].strip().split()
    elif obj.valid and keywords is not None:
        for val in obj.settings[key]:
            if val not in keywords:
                obj.valid = False
                # adp_err
    else:
        obj.valid = False

def check_float_unit(obj,key,units):
    if not obj.valid:
        return
    else:
        try:
            val = float(obj.settings[key].split()[0])
        except:
            obj.valid = False
            err.invalid_parameter(key,obj.settings[key],"float")

    u = obj.settings[key].split()[1]
    if u not in units:
        obj.valid = False
        # break with appropriate error here
    if u != "percent":
        u = adp_constants.units_dict[u]
    else:
        val = val/100
    obj.settings[key] = obj.Tolerance(val,u)
        
    
    

# first item in list is default value
# all further items are check functions to be called on user input
# and kwargs for those check functions
parameters = {
  "equilibration_timesteps"   :[0,[check_int,{"lbound":0}]],
  "calculate_adp"             :[True,[check_boolean,{}]],
  "write_adp"                 :[True,[check_boolean,{}]],
  "r_equilibrium"             :["finite",
                                [check_keywords,{"keywords":["finite","zero"]}]],
  "write_uij"                 :[False,[check_boolean,{}]],
  "calculate_ke"              :[False,[check_boolean,{}]],
  "write_ke"                  :[False,[check_boolean,{}]],
  "element_summary"           :[False,[check_boolean,{}]],
  "output_length"             :["angstrom",
                                [check_keywords,
                                 {"keywords":["angstrom","ang","bohr","nanometer","nm","atomic","atomic_length"],
                                  "map":{"atomic":"atomic_length"}}]],
  "output_energy"             :["electron_volt",
                                [check_keywords,
                                 {"keywords":["electron_volt","electron-volt","ev","hartree","ha","joule","j","atomic","atomic_energy"],
                                  "map":{"atomic":"atomic_energy"}}]],
  "write_jmol"                :[False,[check_boolean,{}]],
  "jmol_scale"                :[5,[check_float,{"lbound":0}]],
  "detect_environment"        :[None,[check_string_arr,{}]],
  "environment_tolerance_adp" :[None,[check_float_unit,{"units":["percent","angstrom"]}]],
  "environment_tolerance_ke"  :[None,[check_float_unit,{"units":["percent","electron-volt"]}]]
}

class Settings:
  def __init__(self,seed):
    self.seed = seed
    self.settings = {}
    self.user_def = {}
    self.output_units = {"length":[adp_constants.ureg.angstrom,"angstrom"],
                         "energy":[adp_constants.ureg.electron_volt,"electron-volt"]}
    self.valid = True
    self.initialise()   # set up settings and user_def dict
    self.parse()        # parse input file
    self.check()        # check values for parameters are valid and fix types
    self.set_units()

  def __repr__(self):
    string = ""
    return string

  def initialise(self):
    """sets up settings and user_def dictionary based on the values in
       parameters dictionary
    """
    for key in parameters.keys():
        self.settings[key] = parameters[key][0]
        self.user_def[key] = False

  def parse(self):
    """Parses the input file <self.seed>.adp
       For any parameter, the value specified will be stored as a string
       in self.settings and user_def will be set to True
    """
    try:
      # read in the settings file
      with open(f"{self.seed}.adp","r") as file:
          data = [line.strip() for line in file.readlines()]
    except:
      # settings file not found
      err.file_not_found(self.seed,".adp")

    for line in data:
      line = line.split("!")[0].strip()   # remove comments
      if len(line) == 0:      # ignore empty lines
        continue

      test = line.lower().split(":")      # all lowercase, : delimiter
      if len(test) == 2:
        k,v = test[0],test[1]
      else:
        err.unexpected_format(line,self.seed)

      k,v = k.strip(),v.strip()           # remove any extra whitespace from key/value

      if k not in self.settings.keys():
        self.valid = False
        err.unknown_keyword(k,self.seed)
      else:
        self.settings[k] = v
        self.user_def[k] = True

  def check(self):
    for key in self.settings.keys():
      if self.user_def[key]:
        for test in parameters[key][1:]:
          test[0](self,key,**test[1])

    if self.settings["write_adp"] and not self.settings["calculate_adp"]:
      self.settings["calculate_adp"] = True
    if self.settings["write_jmol"] and not self.settings["calculate_adp"]:
      self.settings["calculate_adp"] = True
    if self.settings["write_uij"] and not self.settings["calculate_adp"]:
      self.settings["calculate_adp"] = True
    if self.settings["write_ke"] and not self.settings["calculate_ke"]:
      self.settings["calculate_ke"] = True

  def set_units(self):
    units = {"angstrom"     :[adp_constants.ureg.angstrom,"angstrom"],
             "ang"          :[adp_constants.ureg.angstrom,"angstrom"],
             "bohr"         :[adp_constants.ureg.bohr,"bohr"],
             "nanometer"    :[adp_constants.ureg.nanometer,"nm"],
             "nm"           :[adp_constants.ureg.nanometer,"nm"],
             "atomic_length":[adp_constants.ureg.bohr,"bohr"],
             "electron_volt":[adp_constants.ureg.electron_volt,"electron-volt"],
             "electron-volt":[adp_constants.ureg.electron_volt,"electron-volt"],
             "ev"           :[adp_constants.ureg.electron_volt,"electron-volt"],
             "hartree"      :[adp_constants.ureg.hartree,"hartree"],
             "ha"           :[adp_constants.ureg.hartree,"hartree"],
             "atomic_energy":[adp_constants.ureg.hartree,"hartree"],
             "joule"        :[adp_constants.ureg.joule,"joule"],
             "j"            :[adp_constants.ureg.joule,"joule"]}
    self.output_units["length"] = units[self.settings["output_length"]]
    self.output_units["energy"] = units[self.settings["output_energy"]]
    for key in ["output_length","output_energy"]:
      if "atomic" in self.settings[key]:
        self.settings[key] = "atomic"


            
        
