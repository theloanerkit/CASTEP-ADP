import pytest
from castep_adp import settings,adp_constants,err

def check_dict(s,dict):
  errors = 0
  errmsg = ""
  for key in dict.keys():
    if s[key] != dict[key]:
      if errors == 0:
        errmsg += "Error(s) occured:\n"
      errors += 1
      errmsg += f"  setting {key} given default {s[key]}, intended default {dict[key]}"
  return errors, errmsg

def test_default_settings():
  params_dict = {
    "equilibration_timesteps"  : 0,
    "calculate_adp"            : True,
    "write_adp"                : True,
    "r_equilibrium"            : "finite",
    "write_uij"                : False,
    "calculate_ke"             : False,
    "write_ke"                 : False,
    "element_summary"          : False,
    "output_length"            : adp_constants.ureg.angstrom,
    "output_energy"            : adp_constants.ureg.electron_volt,
    "write_jmol"               : False,
    "jmol_scale"               : 5,
    "detect_environment"       : None,
    "environment_tolerance_adp": None,
    "environment_tolerance_ke" : None
  }
  s = settings.Settings("tests/data/test_default")
  errors, errmsg = check_dict(s.settings,params_dict)
  assert errors==0, errmsg

def test_adp_settings_1():
  params_dict = {
    "equilibration_timesteps"  : 0,
    "calculate_adp"            : True,
    "write_adp"                : True,
    "r_equilibrium"            : "finite",
    "write_uij"                : False,
    "calculate_ke"             : False,
    "write_ke"                 : False,
    "element_summary"          : False,
    "write_jmol"               : False,
    "jmol_scale"               : 5,
    "detect_environment"       : None,
    "environment_tolerance_adp": None,
    "environment_tolerance_ke" : None
  }
  s = settings.Settings("tests/data/test_adp_1")
  errors, errmsg = check_dict(s.settings,params_dict)
  assert errors==0, errmsg

def test_adp_settings_2():
  params_dict = {
    "equilibration_timesteps"  : 0,
    "calculate_adp"            : True,
    "write_adp"                : False,
    "r_equilibrium"            : "finite",
    "write_uij"                : True,
    "calculate_ke"             : False,
    "write_ke"                 : False,
    "element_summary"          : False,
    "write_jmol"               : False,
    "jmol_scale"               : 5,
    "detect_environment"       : None,
    "environment_tolerance_adp": None,
    "environment_tolerance_ke" : None
  }
  s = settings.Settings("tests/data/test_adp_2")
  errors, errmsg = check_dict(s.settings,params_dict)
  assert errors==0, errmsg

def test_ke_settings_1():
  params_dict = {
    "equilibration_timesteps"  : 0,
    "calculate_adp"            : True,
    "write_adp"                : True,
    "r_equilibrium"            : "finite",
    "write_uij"                : False,
    "calculate_ke"             : True,
    "write_ke"                 : True,
    "element_summary"          : False,
    "write_jmol"               : False,
    "jmol_scale"               : 5,
    "detect_environment"       : None,
    "environment_tolerance_adp": None,
    "environment_tolerance_ke" : None
  }
  s = settings.Settings("tests/data/test_ke_1")
  errors, errmsg = check_dict(s.settings,params_dict)
  assert errors==0, errmsg

def test_length_unit_angstrom():
  output_units = {"length":[adp_constants.ureg.angstrom,"angstrom"],
                  "energy":[adp_constants.ureg.electron_volt,"electron-volt"]}
  s1 = settings.Settings("tests/data/test_angstrom_1")
  errors, errmsg = check_dict(s1.output_units,output_units)
  assert errors==0, errmsg
  s2 = settings.Settings("tests/data/test_angstrom_2")
  errors, errmsg = check_dict(s2.output_units,output_units)
  assert errors==0, errmsg

def test_length_unit_bohr():
  output_units = {"length":[adp_constants.ureg.bohr,"bohr"],
                  "energy":[adp_constants.ureg.electron_volt,"electron-volt"]}
  s1 = settings.Settings("tests/data/test_bohr_1")
  errors, errmsg = check_dict(s1.output_units,output_units)
  assert errors==0, errmsg
  s2 = settings.Settings("tests/data/test_bohr_2")
  errors, errmsg = check_dict(s2.output_units,output_units)
  assert errors==0, errmsg
  s3 = settings.Settings("tests/data/test_bohr_3")
  errors, errmsg = check_dict(s3.output_units,output_units)
  assert errors==0, errmsg

def test_length_unit_nanometer():
  output_units = {"length":[adp_constants.ureg.nanometer,"nm"],
                  "energy":[adp_constants.ureg.electron_volt,"electron-volt"]}
  s1 = settings.Settings("tests/data/test_nanometer_1")
  errors, errmsg = check_dict(s1.output_units,output_units)
  assert errors==0, errmsg
  s2 = settings.Settings("tests/data/test_nanometer_2")
  errors, errmsg = check_dict(s2.output_units,output_units)
  assert errors==0, errmsg

def test_energy_unit_electron_volt():
  output_units = {"length":[adp_constants.ureg.angstrom,"angstrom"],
                  "energy":[adp_constants.ureg.electron_volt,"electron-volt"]}
  s1 = settings.Settings("tests/data/test_electron_volt_1")
  errors, errmsg = check_dict(s1.output_units,output_units)
  assert errors==0, errmsg
  s2 = settings.Settings("tests/data/test_electron_volt_2")
  errors, errmsg = check_dict(s2.output_units,output_units)
  assert errors==0, errmsg
  s3 = settings.Settings("tests/data/test_electron_volt_3")
  errors, errmsg = check_dict(s3.output_units,output_units)
  assert errors==0, errmsg

def test_energy_unit_hartree():
  output_units = {"length":[adp_constants.ureg.angstrom,"angstrom"],
                  "energy":[adp_constants.ureg.hartree,"hartree"]}
  s1 = settings.Settings("tests/data/test_hartree_1")
  errors, errmsg = check_dict(s1.output_units,output_units)
  assert errors==0, errmsg
  s2 = settings.Settings("tests/data/test_hartree_2")
  errors, errmsg = check_dict(s2.output_units,output_units)
  assert errors==0, errmsg
  s3 = settings.Settings("tests/data/test_hartree_3")
  errors, errmsg = check_dict(s3.output_units,output_units)
  assert errors==0, errmsg
  s4 = settings.Settings("tests/data/test_hartree_4")
  errors, errmsg = check_dict(s4.output_units,output_units)
  assert errors==0, errmsg

def test_energy_unit_hartree():
  output_units = {"length":[adp_constants.ureg.angstrom,"angstrom"],
                  "energy":[adp_constants.ureg.joule,"joule"]}
  s1 = settings.Settings("tests/data/test_joule_1")
  errors, errmsg = check_dict(s1.output_units,output_units)
  assert errors==0, errmsg
  s2 = settings.Settings("tests/data/test_joule_2")
  errors, errmsg = check_dict(s2.output_units,output_units)
  assert errors==0, errmsg

def test_eq_timesteps():
  params_dict = {
    "equilibration_timesteps"  : 100,
    "calculate_adp"            : True,
    "write_adp"                : True,
    "r_equilibrium"            : "finite",
    "write_uij"                : False,
    "calculate_ke"             : False,
    "write_ke"                 : False,
    "element_summary"          : False,
    "output_length"            : adp_constants.ureg.angstrom,
    "output_energy"            : adp_constants.ureg.electron_volt,
    "write_jmol"               : False,
    "jmol_scale"               : 5,
    "detect_environment"       : None,
    "environment_tolerance_adp": None,
    "environment_tolerance_ke" : None
  }
  s1 = settings.Settings("tests/data/test_eq_ts_1")
  errors, errmsg = check_dict(s1.settings,params_dict)
  assert errors==0, errmsg
  with pytest.raises(err.InvalidParamter):
    s2 = settings.Settings("tests/data/test_eq_ts_2")
  with pytest.raises(err.InvalidParamter):
    s3 = settings.Settings("tests/data/test_eq_ts_3")
  with pytest.raises(err.InvalidParamter):
    s4 = settings.Settings("tests/data/test_eq_ts_4")

def test_write_jmol():
  params_dict = {
    "equilibration_timesteps"  : 0,
    "calculate_adp"            : True,
    "write_adp"                : False,
    "r_equilibrium"            : "finite",
    "write_uij"                : False,
    "calculate_ke"             : False,
    "write_ke"                 : False,
    "element_summary"          : False,
    "output_length"            : adp_constants.ureg.angstrom,
    "output_energy"            : adp_constants.ureg.electron_volt,
    "write_jmol"               : True,
    "jmol_scale"               : 5,
    "detect_environment"       : None,
    "environment_tolerance_adp": None,
    "environment_tolerance_ke" : None
  }
  s1 = settings.Settings("tests/data/test_write_jmol_1")
  errors, errmsg = check_dict(s1.settings,params_dict)
  assert errors==0, errmsg
  with pytest.raises(err.InvalidParamter):
    s2 = settings.Settings("tests/data/test_write_jmol_2")

def test_r_equ():
  params_dict = {
    "equilibration_timesteps"  : 0,
    "calculate_adp"            : True,
    "write_adp"                : True,
    "r_equilibrium"            : "zero",
    "write_uij"                : False,
    "calculate_ke"             : False,
    "write_ke"                 : False,
    "element_summary"          : False,
    "output_length"            : adp_constants.ureg.angstrom,
    "output_energy"            : adp_constants.ureg.electron_volt,
    "write_jmol"               : False,
    "jmol_scale"               : 5,
    "detect_environment"       : None,
    "environment_tolerance_adp": None,
    "environment_tolerance_ke" : None
  }
  s1 = settings.Settings("tests/data/test_r_equ_1")
  errors, errmsg = check_dict(s1.settings,params_dict)
  assert errors==0, errmsg
  with pytest.raises(err.InvalidParamter):
    s2 = settings.Settings("tests/data/test_r_equ_2")

def test_detect_env():
  params_dict = {
    "equilibration_timesteps"  : 0,
    "calculate_adp"            : True,
    "write_adp"                : True,
    "r_equilibrium"            : "finite",
    "write_uij"                : False,
    "calculate_ke"             : False,
    "write_ke"                 : False,
    "element_summary"          : False,
    "output_length"            : adp_constants.ureg.angstrom,
    "output_energy"            : adp_constants.ureg.electron_volt,
    "write_jmol"               : False,
    "jmol_scale"               : 5,
    "detect_environment"       : ["h"],
    "environment_tolerance_adp": None,
    "environment_tolerance_ke" : None
  }
  s1 = settings.Settings("tests/data/test_detect_env_1")
  errors, errmsg = check_dict(s1.settings,params_dict)
  assert errors==0, errmsg
  params_dict["detect_environment"].append("o")
  s2 = settings.Settings("tests/data/test_detect_env_2")
  errors, errmsg = check_dict(s2.settings,params_dict)
  assert errors==0, errmsg