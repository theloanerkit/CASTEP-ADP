from castep_adp import settings,adp_constants

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
  errors = 0
  errmsg = ""
  for key in params_dict.keys():
    if s.settings[key] != params_dict[key]:
      if errors == 0:
        errmsg += "Error(s) occured:\n"
      errors += 1
      errmsg += f"  setting {key} given default {s.settings[key]}, intended default {params_dict[key]}"
  assert errors==0, errmsg
