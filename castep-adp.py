import argparse
import os
import pint
from copy import deepcopy as cp
import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as sc
import time
import subprocess
import struct

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
                "Ne": "B3E3F5"}

class md_frame:
    def __init__(self):
        self.time = None
        self.positions = []
        self.energy = None
        self.velocities = []

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

class adp_settings:
    def __init__(self,seed):
        self.seed = seed
        self.settings = {"equilibration_timesteps":0,
                         "write_jmol":False,
                         "r_equilibrium":"finite"} # finite, zero, time
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
        if self.valid and self.settings["write_jmol"] in ["t","true"]:
            self.settings["write_jmol"] = True
        elif self.valid and self.settings["write_jmol"] in ["f","false"]:
            self.settings["write_jmol"] = False
        else:
            self.valid = False
            print(f"unexpected value {self.settings["write_jmol"]} for write_jmol, expected true or false")

        # checking r_equilibrium
        if self.valid and self.settings["r_equilibrium"] not in ["finite", "zero"]: # add time later as an option
            self.valid = False
            print(f"unexpected value {self.settings["r_equilibrium"]} for r_equilibrium, expected one of finite/zero/time")

#def loading_bar(percent,msg):
#    string = f"{msg} [{"#"*percent}{" "*(50-percent)}]"
#    print(string,end="\r")

def parse_md(fname):
    """parses a CASTEP .md file, returning a md_run object with data

    Args:
        fname (str): file name of the .md file

    Returns:
        _md_run: object containing data from the md run
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
    # get velocities
    process_vel = subprocess.Popen(["grep","<-- V",f"{fname}.md"],stdout=subprocess.PIPE)
    all_velocities = np.fromiter(map(float,subprocess.check_output(["awk","{print $3,$4,$5}"],stdin=process_vel.stdout).decode("utf-8").split()),dtype=float)
    # dimensions
    coords = 3
    atom_count = len(atoms)
    timesteps = int(len(all_positions)/(atom_count*coords))
    # reshape and make float
    all_positions = np.reshape(all_positions,(timesteps,atom_count,coords))
    all_velocities = np.reshape(all_velocities,(timesteps,atom_count,coords))

    md = md_run(timesteps,atoms,all_positions,all_velocities)
    return md




def calc_r_eq_from_md(md):
    r_eq = np.zeros((len(md.atoms),3))*md.positions[0,0,0].units
    for timestep in md.positions:
        r_eq += timestep
    r_eq = r_eq.to("angstrom")
    r_eq /= md.timesteps
    return r_eq

def calc_cov_matrix(md,atoms,r_eqm):
    cov_mat = np.zeros((len(atoms),3,3))
    for frame in md:
        for i in range(len(atoms)):
            #cov_mat[i]
            disp1 = np.reshape((np.asarray(frame.positions[i]) - r_eqm[i]),(3,1))
            disp2 = np.reshape(disp1,(1,3))
            cov_mat[i] += np.matmul(disp1,disp2)
    cov_mat /= len(md)
    return cov_mat

def evals_evecs(cov_mat,atoms):
    axes = []
    for i in range(len(atoms)):
        evals, evecs = sc.eig(cov_mat[i])
        temp = []
        for i in range(len(evecs)):
            temp.append(list(evals[i]*np.asarray(evecs[i])))
        axes.append(temp)
    return axes

def get_atom_pos(seed):
    with open(f"{seed}.cell","r") as file:
        cell = [line.strip().lower() for line in file.readlines()]
    frac = False
    if "%block positions_abs" in cell:
        start = cell.index("%block positions_abs")+1
        stop = cell.index("%endblock positions_abs")
    else:
        start = cell.index("%block positions_frac")+1
        stop = cell.index("%endblock positions_frac")
        frac = True

    lattice = []
    if frac:
        start_l = cell.index("%block lattice_cart")+1
        stop_l = cell.index("%endblock lattice_cart")
        for line in cell[start_l:stop_l]:
            lattice.append(list(map(float,line.split())))
    positions = []
    for line in cell[start:stop]:
        pos = list(map(float,line.split()[1:]))
        positions.append(pos)
    temp = []
    if frac:
        for pos in positions:
            temp_pos = np.asarray([0.0,0.0,0.0])
            for i in range(3):
                temp_pos += pos[i]*np.asarray(lattice[i])
                print(temp_pos)
            print()
            temp.append(list(temp_pos))
        positions = temp
    return positions
            
def write_jmol_script(seed,axes,atoms,atom_pos):
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
                file.write(f"{atom_pos[i][j]} ")
            file.write("} scale 5 ")
            file.write(f"color [x{jmol_colours[atoms[i].split("-")[0]]}]")
            file.write("\n")

if __name__ == "__main__":
    args = parser.parse_args()

    settings = adp_settings(args.seed)

    if not os.path.isfile(f"{args.seed}.md"):
        print(f"file not found {args.seed}.md")
        quit()

    print("parsing .md file")
    s = time.time()
    md = parse_md(args.seed)
    print(f"parsing .md: {time.time()-s}")

    r_eq = None
    if settings.settings["r_equilibrium"] == "zero":
        print("reading r_eq")
        # read from cell file
        pass
    elif settings.settings["r_equilibrium"] == "finite":
        print("calculating r_eq")
        s = time.time()
        r_eq = calc_r_eq_from_md(md)
        print(f"calculating r_eq: {time.time()-s}")
    if r_eq is None:
        print("r_eq not calculated")
        quit()
    #if settings.settings["r_equilibrium"] == "zero":
    #    if not os.path.isfile(f"{args.seed}.cell"):
    #        print(f"file not found {args.seed}.cell, r_equilibrium at T=0K requires a .cell file")
    #        quit()
    #    print("r_equilibrium at T=0K is not yet implemented")
    #else:
    #    if settings.settings["r_equilibrium"] == "finite":
    #        r_eqm = calc_r_eqm(md[settings.settings["equilibration_timesteps"]:],atoms)
    #    elif settings.settings["r_equilibrium"] == "time":
    #        print("time-varying r_equilibrium is not yet implemented")
    #        quit()
#
    #print(r_eqm)
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