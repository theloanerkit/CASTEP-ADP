import argparse
import os
import pint
import numpy as np
import scipy.linalg as sc
import subprocess

import time

ureg = pint.UnitRegistry()

parser = argparse.ArgumentParser()
parser.add_argument("-s","--seed",help="seedname for the files (.md/.phonon/.cell/.adp)")

jmol_colours = {"H":  "FFFFFF",
                "He": "D9FFFF",
                "Li": "CC80FF",
                "Be": "C2FF00",
                "B":  "FFB5B5",
                "C":  "909090",
                "N":  "3050F8",
                "O":  "FF0D0D",
                "F":  "90E050",
                "Ne": "B3E3F5",
                "Mo": "54B5B5"}

masses = {"H":  1.00794,
          "He": 4.0026,
          "Li": 6.941,
          "Be": 9.012187,
          "B":  10.811,
          "C":  12.0107,
          "N":  14.00674,
          "O":  15.9994,
          "F":  18.9984,
          "Ne": 20.1797,
          "Na": 22.98977,
          "Mg": 24.305,
          "Al": 26.98154,
          "Si": 28.0855,
          "P":  30.97376,
          "S":  32.066,
          "Cl": 35.4527,
          "Ar": 39.948,
          "K":  39.0983,
          "Ca": 40.078,
          "Mo": 95.94}

mkeys = masses.keys()
for k in mkeys:
    masses[k] = (masses[k]*ureg.a_u_mass)

#masses *= ureg.a_u_mass

class md_run:
    def __init__(self,timesteps,atoms,positions,velocities):
        self.timesteps = timesteps
        self.atoms = atoms
        self.positions = positions
        self.velocities = velocities
        self.apply_units()

    def apply_units(self):
        self.positions*=ureg.a_u_length
        self.velocities*=(ureg.a_u_length/ureg.a_u_time)

class cell:
    def __init__(self,lattice,positions):
        self.lattice = lattice
        self.positions = positions
        self.apply_units()

    def apply_units(self):
        self.lattice*=ureg.angstrom
        self.positions*=ureg.angstrom

class adp_settings:
    def __init__(self,seed):
        self.seed = seed
        self.settings = {"equilibration_timesteps":0,
                         "write_jmol":False,
                         "r_equilibrium":"finite", # finite, zero, time
                         "calculate_adp":True,
                         "jmol_scale":5,
                         "calculate_ke":False
                         } 
        self.valid = True
        self.parse()
        self.check()

    def parse(self):
        with open(f"{self.seed}.adp","r") as file:
            data = [line.strip() for line in file.readlines()]
        
        for line in data:
            try:
                line = line.split("!")[0]
                k,v = line.lower().split(":")
            except:
                print(f"unexpected format in {self.seed}.adp, line: {line}")
                break
            k = k.strip()
            v = v.strip()
            if k not in self.settings.keys():
                print(f"unknown paramter: {k}")
                break
            else:
                self.settings[k] = v

    def check(self):
        # checking equilibration_timesteps
        try:
            self.settings["equilibration_timesteps"] = int(self.settings["equilibration_timesteps"])
        except:
            self.valid = False
            print(f"equilibration_timesteps must be an integer, got {self.settings["equilibration_timesteps"]}")
        if self.valid and self.settings["equilibration_timesteps"] < 0:
            self.valid = False
            print(f"equilibration_timesteps must be non-negative, got {self.settings["equilibration_timesteps"]}")

        # checking write_jmol
        if self.valid and self.settings["write_jmol"] in ["t","true",True]:
            self.settings["write_jmol"] = True
        elif self.valid and self.settings["write_jmol"] in ["f","false",False]:
            self.settings["write_jmol"] = False
        else:
            self.valid = False
            print(f"unexpected value {self.settings["write_jmol"]} for write_jmol, expected true or false")

        # checking r_equilibrium
        if self.valid and self.settings["r_equilibrium"] not in ["finite", "zero"]: # add time later as an option
            self.valid = False
            print(f"unexpected value {self.settings["r_equilibrium"]} for r_equilibrium, expected one of finite/zero/time")

        # checking calculate_adp
        if self.valid and self.settings["calculate_adp"] in ["t","true",True]:
            self.settings["calculate_adp"] = True
        elif self.valid and self.settings["calculate_adp"] in ["f","false",False]:
            self.settings["calcuate_adp"] = False
        else:
            self.valid = False
            print(f"unexpected value {self.settings["calculate_adp"]} for calculate_adp, expected true or false")
        if self.valid and self.settings["write_jmol"] and not self.settings["calculate_adp"]:
            self.valid = False
            print(f"calculate_adp cannot be set to false if write_jmol is true")

        # checking calculate_ke
        if self.valid and self.settings["calculate_ke"] in ["t","true",True]:
            self.settings["calculate_ke"] = True
        elif self.valid and self.settings["calculate_ke"] in ["f","false",False]:
            self.settings["calculate_ke"] = False
        else:
            self.valid = False
            print(f"unexpected value {self.settings["calculate_ke"]} for calculate_ke, expected true or false")


def parse_md(fname,eq_timesteps):
    """parses a CASTEP .md file, returning a md_run object with data

    Args:
        fname (str): file name of the .md file

    Returns:
        md_run: object containing data from the md run
    """
    atoms = []
    # start and stop index for first md timestep
    start = int(subprocess.check_output(["grep","-m","1","-n","<-- R",f"{fname}.md"]).decode("utf-8").strip().split(":")[0])-1
    stop = int(subprocess.check_output(["grep","-m","1","-n","<-- V",f"{fname}.md"]).decode("utf-8").strip().split(":")[0])-1

    atoms = []
    # get atom labels from first md timestep
    with open(f"{fname}.md") as file:
        for i, line in enumerate(file):
            if i >= start and i < stop:
                atoms.append("-".join(line.split()[0:2]))
            if i >= stop:
                break
    # get positions
    process_pos = subprocess.Popen(["grep","<-- R",f"{fname}.md"],stdout=subprocess.PIPE)
    all_positions = np.fromiter(map(float,subprocess.check_output(["awk","{print $3,$4,$5}"],stdin=process_pos.stdout).decode("utf-8").split()),dtype=float)
    #all_positions = all_positions[eq_timesteps:]
    # get velocities
    process_vel = subprocess.Popen(["grep","<-- V",f"{fname}.md"],stdout=subprocess.PIPE)
    all_velocities = np.fromiter(map(float,subprocess.check_output(["awk","{print $3,$4,$5}"],stdin=process_vel.stdout).decode("utf-8").split()),dtype=float)
    #all_velocities = all_velocities[eq_timesteps:]
    # dimensions
    coords = 3
    atom_count = len(atoms)
    timesteps = int(len(all_positions)/(atom_count*coords))
    # reshape and make float
    all_positions = np.reshape(all_positions,(timesteps,atom_count,coords))
    all_velocities = np.reshape(all_velocities,(timesteps,atom_count,coords))

    timesteps -= eq_timesteps
    all_positions = all_positions[eq_timesteps:]
    all_velocities = all_velocities[eq_timesteps:]

    md = md_run(timesteps,atoms,all_positions,all_velocities)
    return md

def parse_cell(fname):
    """parses a CASTEP .cell file, returning a cell object with data

    Args:
        fname (str): file name of the .cell file

    Returns:
        cell: object containing cell file data
    """
    with open(f"{fname}.cell","r") as file:
        cell_file = [line.strip().lower() for line in file.readlines()]
    # need a lattice block and a positions block
    lattice_type = None
    position_type = None
    if "%block lattice_abc" in cell_file:
        lattice_type = "abc"
    elif "%block lattice_cart" in cell_file:
        lattice_type = "cart"
    if "%block positions_frac" in cell_file:
        position_type = "frac"
    elif "%block positions_abs" in cell_file:
        position_type = "abs"

    # get positions in cell file
    l_start = cell_file.index(f"%block lattice_{lattice_type}")+1
    l_stop = cell_file.index(f"%endblock lattice_{lattice_type}")
    p_start = cell_file.index(f"%block positions_{position_type}")+1
    p_stop = cell_file.index(f"%endblock positions_{position_type}")

    # md writes out in cartesians -> need this to be in cartesians too
    # read in the lattice
    if lattice_type == "cart":
        lattice = np.asarray([line.split() for line in cell_file[l_start:l_stop]],dtype=float)
    elif lattice_type == "abc":
        print("lattice_abc not yet supported")
        quit()
    # read in the positions
    position = np.asarray([line.split()[1:4] for line in cell_file[p_start:p_stop]],dtype=float)
    # if positions given in fractional, convert to absolute
    if position_type == "frac":
        for i in range(len(position)):
            position[i] = np.matmul(lattice,position[i])
    
    cell_info = cell(lattice, position)
    return cell_info

def calc_r_eq_from_md(md):
    """calculates r_eq as the average of positions from a .md file

    Args:
        md (md_run): object containing data from CASTEP .md file

    Returns:
        arr: array of atomic equilibrium positions (dimension [number of atoms, 3])
    """
    # set up positions array with units from the .md file (atomic units)
    r_eq = np.zeros((len(md.atoms),3))*md.positions[0,0,0].units
    # add up all positions
    for timestep in md.positions:
        r_eq += timestep
    # convert to angstroms
    r_eq = r_eq.to("angstrom")
    # average over timesteps
    r_eq /= md.timesteps
    return r_eq

def get_r_eq_from_cell(seed):
    cell = parse_cell(seed)
    return cell.positions

def calculate_covariance_matrix(md,r_eq):
    covariance_matrix = np.zeros((len(md.atoms),3,3))*(md.positions[0,0,0].units**2)
    for timestep in md.positions:
        disp = timestep-r_eq
        covariance_matrix += np.multiply(np.expand_dims(disp,2),np.expand_dims(disp,1))
    covariance_matrix /= md.timesteps
    return covariance_matrix

def evals_evecs(cov_mat,atoms):
    axes = []
    for i in range(len(atoms)):
        evals, evecs = sc.eig(cov_mat[i])
        temp = []
        for i in range(len(evecs)):
            temp.append(list(evals[i]*np.asarray(evecs[i])))
        axes.append(temp)
    return axes

def calculate_kinetic_energy_tensor(md):
    v_tensor = np.zeros((len(md.atoms),3,3))*(md.velocities[0,0,0].units**2)
    for timestep in md.velocities:
        v_tensor += np.multiply(np.expand_dims(timestep,2),np.expand_dims(timestep,1))
    v_tensor /= md.timesteps
    v_tensor *= 0.5
    ke_tensor = np.zeros(np.shape(v_tensor))*(v_tensor[0,0,0].units*masses["H"].units)
    for i in range(len(md.atoms)):
        ke_tensor[i,:,:] = v_tensor[i,:,:]*masses[md.atoms[i].split("-")[0]]
    ke_tensor = ke_tensor.to("electron_volt")
    return ke_tensor

            
def write_jmol_script(seed,axes,atoms,atom_pos,scale):
    with open(f"{seed}_ellipsoid.spt","w") as file:
        for i in range(len(axes)):
            file.write(f"ellipsoid ID {atoms[i].replace("-","")} AXES ")
            for j in range(len(axes[i])):
                file.write("{ ")
                for k in range(len(axes[i][j])):
                    file.write(f"{float(axes[i][j][k])} ")
                file.write("} ")
            file.write("center { ")
            for j in range(len(atom_pos[i])):
                file.write(f"{atom_pos[i][j].magnitude} ")
            file.write(f"}} scale {scale} ")
            file.write(f"color [x{jmol_colours[atoms[i].split("-")[0]]}]")
            file.write("\n")

if __name__ == "__main__":
    args = parser.parse_args()

    settings = adp_settings(args.seed)

    if not os.path.isfile(f"{args.seed}.md"):
        print(f"file not found {args.seed}.md")
        quit()

    print("parsing .md file")
    md = parse_md(args.seed,settings.settings["equilibration_timesteps"])

    r_eq = None
    if settings.settings["r_equilibrium"] == "zero":
        print("reading r_eq")
        r_eq = get_r_eq_from_cell(args.seed)
    elif settings.settings["r_equilibrium"] == "finite":
        print("calculating r_eq")
        r_eq = calc_r_eq_from_md(md)
    if r_eq is None:
        print("r_eq not calculated")
        quit()

    if settings.settings["calculate_adp"]:
        print("calculating covariance matrix")
        covariance_matrix = calculate_covariance_matrix(md,r_eq)
        print("calculating adp axes")
        axes = evals_evecs(covariance_matrix,md.atoms)

    if settings.settings["write_jmol"]:
        print("writing jmol")
        write_jmol_script(args.seed,axes,md.atoms,r_eq,settings.settings["jmol_scale"])

    if settings.settings["calculate_ke"]:
        print("calculating ke")
        ke_tensor = calculate_kinetic_energy_tensor(md)
        with open(f"{args.seed}_ke_tensor.dat","w") as file:
            for i in range(len(md.atoms)):
                file.write(f"{md.atoms[i]}\n")
                for row in ke_tensor[i]:
                    string = "    ".join(list(map(str,row.magnitude)))
                    file.write(f"{string}\n")
                file.write("\n\n")
    #print(ke_tensor)
#
    #cov_mat = calc_cov_matrix(md[settings.settings["equilibration_timesteps"]:],atoms,r_eqm)
    ##print(f"covariance matrix calc")
    #evecs = evals_evecs(cov_mat,atoms)
    #if settings.settings["write_jmol"]:
    #    write_jmol_script(args.seed,evecs,atoms,((r_eqm*ureg.a_u_length).to("angstrom")).magnitude)
    #print("writing jmol")
    #if write_jmol:
    #    pos = get_atom_pos(args.seed)
    #    #print(pos)
    #    write_jmol_script(args.seed,evecs,atoms,pos)
    #print(md_dict.keys())
    #print(md_dict["ions"],len(md_dict["ions"]))