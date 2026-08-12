import pytest
from castep_adp import constants, err, parse
import numpy as np

def check_md_obj(obj,timesteps,atoms):
  errors = 0
  errmsg = "\n"
  if obj.timesteps != timesteps:
    errors += 1
    errmsg += f"number of timesteps does not match: got {obj.timesteps}, expected {timesteps}\n"
  if obj.atoms != atoms:
    errors += 1
    errmsg += f"atoms do not match: got {obj.atoms}, expected {atoms}\n"
  if obj.positions.units != constants.ureg.bohr:
    errors += 1
    errmsg += f"units of position do not match, got {obj.positions.units}, expected {constants.ureg.bohr}\n"
  if obj.velocities.units != constants.ureg.bohr/constants.ureg.atomic_unit_of_time:
    errors += 1
    errmsg += f"units of velocity do not match, got {obj.velocities.units}, expected {constants.ureg.bohr/constants.ureg.atomic_unit_of_time}\n"
  if np.shape(obj.positions) != (timesteps,len(atoms),3):
    errors += 1
    errmsg += f"shape of positions array does not match, got {np.shape(obj.positions)}, expected {(timesteps,len(atoms),3)}\n"
  if np.shape(obj.velocities) != (timesteps,len(atoms),3):
    errors += 1
    errmsg += f"shape of velocities array does not match, got {np.shape(obj.velocities)}, expected {(timesteps,len(atoms),3)}\n"
  return errors, errmsg

def check_cell_obj(obj,cart,abs):
  errors = 0
  errmsg = "\n"
  if not np.all(obj.lattice_cart == cart):
    errors += 1
    errmsg += f"lattice does not match: \n  got      {obj.lattice_cart}\n  expected {cart}\n"
  for key in abs.keys():
    if key not in obj.positions_abs:
      errors += 1
      errmsg += f"atoms do not match: {key} not found\n"
    else:
      if not np.all(abs[key]==obj.positions_abs[key]):
        errors += 1
        errmsg += f"atoms positions do not match for {key}: \n  got      {obj.positions_abs[key]}\n  expected {abs[key]}\n"
  return errors, errmsg

def test_N2_md_no_equ():
  n2 = parse.parse_md("tests/data/test_N2")
  errors,errmsg = check_md_obj(n2,2,["N 1","N 2"])
  assert errors==0,errmsg

def test_N2_md():
  n2 = parse.parse_md("tests/data/test_N2",1)
  errors,errmsg = check_md_obj(n2,1,["N 1","N 2"])
  assert errors==0,errmsg

def test_N2_md_too_much_equ():
  with pytest.raises(err.ParseMDError):
    n2 = parse.parse_md("tests/data/test_N2",4)


def test_BTO_cart_frac():
  cart = np.diag([4.02,4.02,4.02],k=0)*constants.ureg.angstrom
  abs = {"O 1": np.asarray([2.01,2.01,0])*constants.ureg.angstrom,
         "O 2": np.asarray([2.01,0,2.01])*constants.ureg.angstrom,
         "O 3": np.asarray([0,2.01,2.01])*constants.ureg.angstrom,
         "Ti 1": np.asarray([2.01,2.01,2.01])*constants.ureg.angstrom,
         "Ba 1": np.asarray([0,0,0])*constants.ureg.angstrom}
  bto = parse.parse_cell("tests/data/test_BaTiO_cart_frac")
  errors,errmsg = check_cell_obj(bto,cart,abs)
  assert errors==0,errmsg

def test_BTO_cart_abs():
  cart = np.diag([4.02,4.02,4.02],k=0)*constants.ureg.angstrom
  abs = {"O 1": np.asarray([2.01,2.01,0])*constants.ureg.angstrom,
         "O 2": np.asarray([2.01,0,2.01])*constants.ureg.angstrom,
         "O 3": np.asarray([0,2.01,2.01])*constants.ureg.angstrom,
         "Ti 1": np.asarray([2.01,2.01,2.01])*constants.ureg.angstrom,
         "Ba 1": np.asarray([0,0,0])*constants.ureg.angstrom}
  bto = parse.parse_cell("tests/data/test_BaTiO_cart_abs")
  errors,errmsg = check_cell_obj(bto,cart,abs)
  assert errors==0,errmsg

def test_BTO_abc():
  with pytest.raises(err.IncompatibleCell):
    bto = parse.parse_cell("tests/data/test_BaTiO_abc")

def test_BTO_no_lattice():
  with pytest.raises(err.IncompatibleCell):
    bto = parse.parse_cell("tests/data/test_BaTiO_frac")

def test_BTO_no_positions():
  with pytest.raises(err.IncompatibleCell):
    bto = parse.parse_cell("tests/data/test_BaTiO_cart")

def test_BTO_cart_units_ang():
  cart = np.diag([4.02,4.02,4.02],k=0)*constants.ureg.angstrom
  abs = {"O 1": np.asarray([2.01,2.01,0])*constants.ureg.angstrom,
         "O 2": np.asarray([2.01,0,2.01])*constants.ureg.angstrom,
         "O 3": np.asarray([0,2.01,2.01])*constants.ureg.angstrom,
         "Ti 1": np.asarray([2.01,2.01,2.01])*constants.ureg.angstrom,
         "Ba 1": np.asarray([0,0,0])*constants.ureg.angstrom}
  bto = parse.parse_cell("tests/data/test_BaTiO_cart_ang")
  errors,errmsg = check_cell_obj(bto,cart,abs)
  assert errors==0,errmsg

def test_BTO_cart_units_bohr():
  cart = (np.diag([10,10,10],k=0)*constants.ureg.bohr).to("angstrom")
  abs = {"O 1": (np.asarray([5,5,0])*constants.ureg.bohr).to("angstrom"),
         "O 2": (np.asarray([5,0,5])*constants.ureg.bohr).to("angstrom"),
         "O 3": (np.asarray([0,5,5])*constants.ureg.bohr).to("angstrom"),
         "Ti 1": (np.asarray([5,5,5])*constants.ureg.bohr).to("angstrom"),
         "Ba 1": np.asarray([0,0,0])*constants.ureg.angstrom}
  bto = parse.parse_cell("tests/data/test_BaTiO_cart_bohr")
  errors,errmsg = check_cell_obj(bto,cart,abs)
  assert errors==0,errmsg

def test_BTO_abs_units_ang():
  cart = np.diag([4.02,4.02,4.02],k=0)*constants.ureg.angstrom
  abs = {"O 1": np.asarray([2.01,2.01,0])*constants.ureg.angstrom,
         "O 2": np.asarray([2.01,0,2.01])*constants.ureg.angstrom,
         "O 3": np.asarray([0,2.01,2.01])*constants.ureg.angstrom,
         "Ti 1": np.asarray([2.01,2.01,2.01])*constants.ureg.angstrom,
         "Ba 1": np.asarray([0,0,0])*constants.ureg.angstrom}
  bto = parse.parse_cell("tests/data/test_BaTiO_abs_ang")
  errors,errmsg = check_cell_obj(bto,cart,abs)
  assert errors==0,errmsg

def test_BTO_abs_units_nm():
  cart = np.diag([40.2,40.2,40.2],k=0)*constants.ureg.angstrom
  abs = {"O 1": np.asarray([20.1,20.1,0])*constants.ureg.angstrom,
         "O 2": np.asarray([20.1,0,20.1])*constants.ureg.angstrom,
         "O 3": np.asarray([0,20.1,20.1])*constants.ureg.angstrom,
         "Ti 1": np.asarray([20.1,20.1,20.1])*constants.ureg.angstrom,
         "Ba 1": np.asarray([0,0,0])*constants.ureg.angstrom}
  bto = parse.parse_cell("tests/data/test_BaTiO_abs_nm")
  errors,errmsg = check_cell_obj(bto,cart,abs)
  assert errors==0,errmsg

def test_BTO_cart_units_unknown():
  with pytest.raises(err.UnitError):
    bto = parse.parse_cell("tests/data/test_BaTiO_cart_furlong")

def test_BTO_cart_units_electronvolt():
  with pytest.raises(err.UnitError):
    bto = parse.parse_cell("tests/data/test_BaTiO_cart_electronvolt")
