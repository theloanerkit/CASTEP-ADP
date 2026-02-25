import argparse
import pint
import numpy as np
import scipy.linalg as sp
import os
# -----------------
import adp_constants
import adp_settings
import adp_err
import adp_io
import adp_parse
import adp_obj

# set up dictionary with masses of elements -> move to constants
masses = {}
for key in adp_constants.masses.keys():
    masses[key] = adp_constants.masses[key] * adp_constants.ureg.unified_atomic_mass_unit

def calc_r_eq_from_md(md):
  # initialise array
  r_eq = np.zeros((len(md.atoms),3))*md.positions[0,0,0].units

  # add up all positions
  for timestep in md.positions:
    r_eq += timestep

  # average over timesteps
  r_eq /= md.timesteps
  return r_eq

def reorder_positions(positions,atoms):
  ordered_pos = np.zeros((len(atoms),3),dtype=float)*adp_constants.ureg.angstrom
  for i in range(len(atoms)):
    ordered_pos[i] += positions[atoms[i]]
  return ordered_pos

def calc_covariance_matrix(md,r_eq):
  # initialise array
  cov_mat = np.zeros((len(md.atoms),3,3))*(md.positions[0,0,0].units**2)

  # add up |disp><disp| matrices
  for timestep in md.positions:
    disp = timestep - r_eq
    cov_mat += np.multiply(np.expand_dims(disp,2),np.expand_dims(disp,1))
    
  # average over timesteps
  cov_mat /= md.timesteps
  return cov_mat

def calc_ke_tensor(md):
  # initialise array
  v_tensor = np.zeros((len(md.atoms),3,3))*(md.velocities[0,0,0].units**2)
  ke_tensor = np.zeros(np.shape(v_tensor))*(v_tensor[0,0,0].units*masses["H"].units)

  # add up |v><v| matrices
  for timestep in md.velocities:
    v_tensor += np.multiply(np.expand_dims(timestep,2),np.expand_dims(timestep,1))
    
  # average over timesteps
  v_tensor /= md.timesteps

  # calculate ke tensor
  for i in range(len(md.atoms)):
    ke_tensor[i,:,:] = v_tensor[i,:,:] * 0.5 * masses[md.atoms[i].split()[0]]
  return ke_tensor

def evals_evecs(matrix,atoms,sqrt=False):
  # initialise array
  axes = np.zeros((len(atoms),3,3),dtype=float)

  # calculate eigenvalues and vectors for each atom
  for i in range(len(atoms)):
    evals,evecs = np.linalg.eig(matrix[i].magnitude)
    if sqrt:
      evals = evals**0.5
    for j in range(3):
      axes[i,j,:] += evals[j]*evecs[j,:]
  return axes

#def detect_environments_o(atom,data,tol):
#    atoms = []
#    #tol = adp_obj.Tolerance(0.2,"percent") # 20 percent
#    #lookup = 
#    for i in range(len(data["atoms"])):
#        element = data["atoms"][i].split()[0].lower()
#        if element == atom:
#            a = adp_obj.Atom_Environment(data["atoms"][i],data["adp"][i])
#            atoms.append(a)
#    environments = []
#    for a in atoms:
#        if len(environments) == 0:
#            env = set()
#            env.add(a)
#            environments.append(env)
#        else:
#            #print(f"checking environment of {a.name}")
#            in_any_env = False
#            for env in environments:
#                in_env = True
#                # check magnitudes:
#                for test in env:
#                    for i in range(3):
#                        mag = tol.within_tolerance(test.magnitudes[i],a.magnitudes[i])
#                        if not mag:
#                            in_env = False
#                            break
#                    if not in_env:
#                        break
#                if in_env:
#                    in_any_env = True
#                    env.add(a)
#                    #print(f"   found environment matching {a.name}")
#                    break
#            if not in_any_env:
#                env = set()
#                env.add(a)
#                environments.append(env)
#                #print(f"   not found an environment for {a.name}, creating a new one")
#    return environments
                
def detect_environments(elem,atoms,tol_adp,tol_ke):
  print("hello")
  environments = []
  selected_atoms = []     # atoms we are checking environment for based on input parameter
  for key in atoms.keys():
    #print(key,key.split()[0].lower,elem)
    if key.split()[0].lower() == elem:
      atoms[key].calc_magnitudes(tol_adp,tol_ke)
      selected_atoms.append(atoms[key])
  print(f"checking {len(selected_atoms)} atoms")
  for atom in selected_atoms:
    #print(f"checking atom {atom.name}")
    if len(environments) == 0:
      env = set()
      env.add(atom)
      environments.append(env)
    else:
      in_any_env = False
      for env in environments:
        in_env = True
        for other_atom in env:
          #print(f"other atom: {other_atom.name}")
          #print(f"  adps: {other_atom.adp_magnitudes}   {other_atom.adp_sort}")
          #print(f"  kes:  {other_atom.ke_magnitudes}   {other_atom.ke_sort}")
          adp_match = True
          ke_match = True
          rot_match = True
          for i in range(3):
            if tol_adp is not None:
              adp_match = tol_adp.within_tolerance(other_atom.adp_magnitudes[i],atom.adp_magnitudes[i])
            if tol_ke is not None:
              ke_match = tol_ke.within_tolerance(other_atom.ke_magnitudes[i],atom.ke_magnitudes[i])
            if tol_ke is not None and tol_adp is not None:
              rot_match = other_atom.adp_ke_map == atom.adp_ke_map
          if not (adp_match and ke_match and rot_match):
            in_env = False
            #print(adp_match,ke_match)
            break
        if in_env:
          #print(f"atom: {atom.name}")
          #print(f"  adps: {atom.adp_magnitudes}   {atom.adp_sort}")
          #print(f"  kes:  {atom.ke_magnitudes}   {atom.ke_sort}")
          in_any_env = True
          env.add(atom)
          break
      if not in_any_env:
        env = set()
        env.add(atom)
        environments.append(env)
    #input("next")
  for env in environments:
    print("environment: ")
    for atom in env:
      #m = get_map(atom.adp_sort,atom.ke_sort)
      #s = ""
      #for k in m.keys():
      #   s += f"{k} -> {m[k]}    "
      print(f"    atom: {atom.name}    {atom.adp_sort}, {atom.ke_sort}")
      print(f"        : {atom.adp_magnitudes}    {atom.ke_magnitudes}")
  return environments
  
def get_map(one,two):
  m = {0:-1,1:-1,2:-1}
  for i in range(3):
    m[one[i]] = two[i]
  return m

def get_parser() -> argparse.ArgumentParser:
    # for parsing command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("seed",
                        help="seedname for the files (.md/.phonon/.cell/.adp)")
    parser.add_argument("-d","--dryrun",
                        action="store_true")
    return parser


def main() -> None:
    args = get_parser().parse_args()

    print("loading settings")
    settings = adp_settings.Settings(args.seed)

    data = {"atoms":None,
            "r_eq" :None,
            "adp"  :None,
            "uij"  :None,
            "ke"   :None}
    
    atoms = {}

    if not args.dryrun:
        if not os.path.isfile(f"{args.seed}.md"):
            adp_err.file_not_found(args.seed,".md")

        print("parsing md file")
        md_obj = adp_parse.parse_md(args.seed,settings.settings["equilibration_timesteps"])
        data["atoms"] = md_obj.atoms
        for label in md_obj.atoms:
            atoms[label] = adp_obj.Atom(label)
        keys = list(atoms.keys())

        if settings.settings["detect_environment"] is not None:
            adp_settings.check_string_arr(settings,"detect_environment",md_obj.atoms)

        if settings.settings["r_equilibrium"] == "zero":
            print("reading r_eq from cell")
            if not os.path.isfile(f"{args.seed}.cell"):
                adp_err.file_not_found(args.seed,".cell")
            cell_obj = adp_parse.parse_cell(args.seed)
            for key in cell_obj.positions_abs.keys():
                atoms[key].r_eq = cell_obj.positions_abs[key]
            data["r_eq"] = reorder_positions(cell_obj.positions_abs,data["atoms"])

        elif settings.settings["r_equilibrium"] == "finite":
            print("calculating r_eq")
            r_eq = calc_r_eq_from_md(md_obj)
            data["r_eq"] = r_eq
            for i in range(len(keys)):
                atoms[keys[i]].r_eq = r_eq[i]

        if settings.settings["calculate_adp"]:
            print("calculating covariance matrix")
            cov_mat = calc_covariance_matrix(md_obj,data["r_eq"])
            data["uij"] = cov_mat
            

            print("calculating adp axes")
            unit = cov_mat.units**0.5
            axes = evals_evecs(cov_mat,md_obj.atoms,sqrt=True)
            axes *= unit
            data["adp"] = axes
            for i in range(len(keys)):
                atoms[keys[i]].uij = cov_mat[i]
                atoms[keys[i]].adp = axes[i]

        

        if settings.settings["calculate_ke"]:
            print("calculating ke")
            ke = calc_ke_tensor(md_obj)
            data["ke"] = ke
            unit = ke.units
            axes = evals_evecs(ke,md_obj.atoms)
            axes *= unit
            for i in range(len(keys)):
                atoms[keys[i]].ke = ke[i]
                atoms[keys[i]].ke_vecs = axes[i]

        environments = None
        if settings.settings["detect_environment"] is not None:
            print("searching for chemical environments")
            for atom in settings.settings["detect_environment"]:
                environments = detect_environments(atom,atoms,settings.settings["environment_tolerance_adp"],settings.settings["environment_tolerance_ke"])

        

    print("writing")
    adp_io.write_out(args.seed,settings,atoms,args.dryrun)
    if not args.dryrun:
        if settings.settings["write_jmol"]:
            print("writing jmol")
            adp_io.write_jmol_script(args.seed,
                                     data["adp"].to("angstrom").magnitude,
                                     data["atoms"],
                                     data["r_eq"].to("angstrom").magnitude,
                                     settings.settings["jmol_scale"])
        if settings.settings["detect_environment"] is not None:
            print("writing environments jmol")
            print(f"{len(environments)} environments found")
            adp_io.write_jmol_script(args.seed,
                                     data["adp"].to("angstrom").magnitude,
                                     data["atoms"],
                                     data["r_eq"].to("angstrom").magnitude,
                                     settings.settings["jmol_scale"],
                                     environments)
            with open(f"{args.seed}_env.out","w") as file:
              for i in range(len(environments)):
                file.write(f"environment {i}: ")
                for atom in environments[i]:
                  file.write(f"{atom.name}  ")
                file.write("\n")


if __name__=="__main__":
    main()
