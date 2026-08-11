import numpy as np

from . import err
from .adp_constants import units_dict, ureg


class MD:
  def __init__(self,timesteps,atoms,positions,velocities):
    self.timesteps = timesteps
    self.atoms = atoms
    self.positions = positions
    self.velocities = velocities
    self.apply_units()

  def __repr__(self):
    string = "MD object: \n"
    string += f"    number of timesteps: {self.timesteps}\n"
    string += f"    number of atoms: {len(self.atoms)}\n"
    string += f"    position units: {self.positions.units}\n"
    string += f"    velocity units: {self.velocities.units}"
    return string

  def apply_units(self):
    self.positions *= ureg.atomic_unit_of_length
    self.velocities *= (ureg.atomic_unit_of_length/ureg.atomic_unit_of_time)

class Cell:
  def __init__(self,lattice,positions):
    self.lattice = lattice
    self.positions = positions
    self.lattice_cart = None
    self.positions_abs = None
    self.positions_frac = None
    self.construct_cell()
    if self.positions_abs is None:
      self.get_positions_abs()

  def get_positions_abs(self):
    positions = {}
    for key in self.positions_frac.keys():
      pos = np.zeros(3,dtype=float)*ureg.angstrom
      for i in range(3):
        pos += self.positions_frac[key][i] * self.lattice_cart[i]
      positions[key] = pos
    self.positions_abs = positions

  def construct_cell(self):
    if len(self.lattice) == 0:
      raise err.NoLatticeCartError()
    if len(self.positions) == 0:
      raise err.NoPositionsBlockError()
    if "cart" in self.lattice[0].lower():
      self.construct_cart()
    elif "abc" in self.lattice[0].lower():
      raise err.NoLatticeCartError()
    if "abs" in self.positions[0].lower():
      self.construct_abs()
    elif "frac" in self.positions[0].lower():
      self.construct_frac()

  def construct_cart(self):
    start = 1
    unit = ureg.angstrom
    if len(self.lattice) == 6:
      if self.lattice[1] in units_dict.keys():
        unit = units_dict[self.lattice[1]][0]
      else:
        raise err.UnitError(kind="unknown unit",loc="cell file",unit=self.lattice[1])
      if unit.dimensionality != "[length]":
        raise err.UnitError(kind="dimensionality",unit=unit,dim="[length]")
      start += 1
    lattice = np.zeros((3,3),dtype=float)*ureg.angstrom
    for i in range(start,start+3):
      lattice[i-start] = np.asarray(self.lattice[i].split(),dtype=float)*unit
    self.lattice_cart = lattice

  def construct_abs(self):
    positions = {}
    elements = {}
    unit = ureg.angstrom
    for line in self.positions[1:-1]:
      words = line.split()
      if len(words)==1:
        if words[0] in units_dict.keys():
          unit = units_dict[words[0]][0]
        else:
          raise err.UnitError(kind="unknown unit",loc="cell file",unit=self.lattice[1])
        if unit.dimensionality != "[length]":
          raise err.UnitError(kind="dimensionality",unit=unit,dim="[length]")
        continue
      if words[0] not in elements.keys():
        elements[words[0]] = 0
      elements[words[0]] += 1
      name = f"{words[0]} {elements[words[0]]}"
      pos = np.asarray(words[1:4],dtype=float)*unit
      positions[name] = pos
    self.positions_abs = positions

  def construct_frac(self):
    positions = {}
    elements = {}
    for line in self.positions[1:-1]:
      words = line.split()
      if words[0] not in elements.keys():
        elements[words[0]] = 0
      elements[words[0]] += 1
      name = f"{words[0]} {elements[words[0]]}"
      pos = np.asarray(words[1:4],dtype=float)
      positions[name] = pos
    self.positions_frac = positions

class Tolerance:
  def __init__(self,value,unit):
    self.value = value
    self.unit = unit

  def within_tolerance(self,val1,val2):
    if isinstance(val1,(int,float)) and isinstance(val2,(int,float)):
      pass # both dimensionless -> fine (unless tolerance is not?)
    elif isinstance(val1,(int,float)) and not isinstance(val2,(int,float)):
      # val1: dimensionless, val2: has dimension
      raise err.UnitError(kind="dimensionality",dim="dimensionless",unit=val2.units)
    elif isinstance(val2,(int,float)) and not isinstance(val1,(int,float)):
      # val1: has dimension, val2: dimensionless
      raise err.UnitError(kind="dimensionality",dim="dimensionless",unit=val1.units)
    elif val1.units.dimensionality != val2.units.dimensionality:
      raise err.UnitError(kind="dimensionality",dim=val1.units.dimensionality,unit=val2.units)
    if self.unit == "percent":
      tol = self.value*val1
    else:
      tol = self.value*self.unit

      val1 = val1.to(tol.units)
      val2 = val2.to(tol.units)

    if abs(val2-val1) < tol:
      return True
    else:
      return False

class Atom:
  def __init__(self,name):
    self.name = name
    self.r_eq = None    # equilibrium positions
    self.adp = None     # vectors defining atomic displacement parameter
    self.uij = None     # ADP tensor
    self.ke = None      # KE tensor

    self.ke_vecs = None # vectors defining ke tensor
    self.adp_magnitudes,self.adp_sort = None,None#self.get_sorted_magnitudes()
    self.ke_magnitudes,self.ke_sort = None,None
    self.adp_ke_map = None

  def calc_magnitudes(self,adp,ke):
    if adp is not None:
      self.adp_magnitudes, self.adp_sort = self.get_sorted_magnitudes(self.adp)
    if ke is not None:
      self.ke_magnitudes, self.ke_sort = self.get_sorted_magnitudes(self.ke_vecs)
    if adp is not None and ke is not None:
      m = {0:-1,1:-1,2:-1}
      for i in range(3):
        m[self.adp_sort[i]] = self.ke_sort[i]
      self.adp_ke_map = m

  def get_sorted_magnitudes(self,vecs):
    magnitudes = [np.linalg.norm(v.magnitude) for v in vecs]
    idxs = np.argsort(magnitudes)
    magnitudes = magnitudes*vecs.units
    return magnitudes[idxs], idxs

class AtomEnvironment:
  def __init__(self,label,adp):
    self.name = label
    self.adp = adp
    self.magnitudes, self.sort_arr = self.get_sorted_magnitudes()

  def get_sorted_magnitudes(self):
    magnitudes = [np.linalg.norm(adp.magnitude) for adp in self.adp]
    idxs = np.argsort(magnitudes)
    magnitudes = magnitudes*self.adp.units
    return magnitudes[idxs], idxs
