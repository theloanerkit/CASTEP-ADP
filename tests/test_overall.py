from castep_adp import cli, obj, parse, settings, io

def test_module_imports():
  import castep_adp

def test_run():
  seed = "tests/data/test_N2_2"
  user_settings = settings.Settings(seed)
  md = parse.parse_md(seed,user_settings.settings["equilibration_timesteps"])
  atoms = {}
  for label in md.atoms:
    atoms[label] = obj.Atom(label)
  keys = list(atoms.keys())

  r_eq = cli.calc_r_eq_from_md(md)
  ke = cli.calc_ke_tensor(md)

  for i in range(len(keys)):
    atoms[keys[i]].r_eq = r_eq[i]
    atoms[keys[i]].ke = ke[i]

  io.write_out(seed,user_settings,atoms,False)
