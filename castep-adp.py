import argparse
import os
import pint
from copy import deepcopy as cp
import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as sc
import colours as c

ureg = pint.UnitRegistry()

parser = argparse.ArgumentParser()
parser.add_argument("-s","--seed",help="seedname for the files (.md/.phonon/.cell/.adp)")

class md_frame:
    def __init__(self):
        self.time = None
        self.positions = []
        self.energy = None
        self.velocities = []

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
        if self.valid and self.settings["r_equilibrium"] not in ["finite", "zero", "time"]:
            self.valid = False
            print(f"unexpected value {self.settings["r_equilibrium"]} for r_equilibrium, expected one of finite/zero/time")

def parse_md(fname):
    # read in .md file
    with open(f"{fname}.md","r") as file:
        data = [line.strip() for line in file.readlines()]
    
    # skip header
    start = data.index("END header")
    data = data[start:]

    get_atoms = True
    atoms = []
    frame = None
    frames = []

    for i in range(len(data)):
        if len(data[i]) == 0: # empty line -> new md timestep
            if frame is not None: # we have done at least one timestep
                get_atoms = False
                frames.append(frame)
            frame = md_frame()
        if "<-- R" in data[i]: # position line
            line = data[i].split()
            if get_atoms: # we haven't recorded the atoms yet
                atoms.append("-".join(line[0:2]))
            frame.positions.append(list(map(float,line[2:5]))) # add position to frame position array
        if "<-- V" in data[i]: # velocity line
            line = data[i].split()
            frame.velocities.append(list(map(float,line[2:5]))) # add velocity to frame velocity array
        if "<-- E" in data[i]: # energy line
            frame.time = float(data[i-1]) # time is on the line before energy
            frame.energy = float(data[i].split()[0])
    return frames, atoms

#def parse_files(seed):
#    print("parsing md file")
#    frames, atoms = parse_md(seed)
#
#    #with open(f"{seed}.md","r") as file:
#    #    data = [line.strip() for line in file.readlines()]
#    #atoms = []
#    #get_atoms = True
#    #positions = []
#    #velocities = []
#    #frames = []
#    #frame = None
#    #start = data.index("END header")
#    #data = data[start:]
#    #for i in range(len(data)):
#    #    if len(data[i]) == 0:
#    #        if frame is not None:
#    #            get_atoms = False
#    #            frame.positions = cp(positions)
#    #            frame.velocities = cp(velocities)
#    #            frames.append(frame)
#    #            positions = []
#    #            velocities = []
#    #        frame = md_frame()
#    #    if "<-- R" in data[i]:
#    #        line = data[i].split()
#    #        if get_atoms:
#    #            atoms.append("".join(line[0:2]))
#    #        pos = list(map(float,line[2:5]))
#    #        positions.append(pos)
#    #    if "<-- V" in data[i]:
#    #        line = data[i].split()
#    #        vel = list(map(float,line[2:5]))
#    #        velocities.append(vel)
#    #    if "<-- E" in data[i]:
#    #        frame.time = (float(data[i-1])*ureg.a_u_time).to("s").magnitude
#    #        frame.energy = (float(data[i].split()[0])*ureg.a_u_energy).to("electron_volt").magnitude
#    return frames, atoms

def calc_r_eqm(md,atoms):
    r_eqm = np.zeros((len(atoms),3))
    for frame in md:
        for i in range(len(atoms)):
            r_eqm[i] += frame.positions[i]
    r_eqm /= len(md)

    return r_eqm

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
            file.write(f"color [x{c.jmol_colours[atoms[i].split("-")[0]]}]")
            file.write("\n")

if __name__ == "__main__":
    args = parser.parse_args()

    settings = adp_settings(args.seed)
    if not os.path.isfile(f"{args.seed}.md"):
        print(f"file not found {args.seed}.md")
        quit()

    md, atoms = parse_md(args.seed)
    if settings.settings["r_equilibrium"] == "zero":
        if not os.path.isfile(f"{args.seed}.cell"):
            print(f"file not found {args.seed}.cell, r_equilibrium at T=0K requires a .cell file")
            quit()
        print("r_equilibrium at T=0K is not yet implemented")
    else:
        if settings.settings["r_equilibrium"] == "finite":
            r_eqm = calc_r_eqm(md[settings.settings["equilibration_timesteps"]:],atoms)
        elif settings.settings["r_equilibrium"] == "time":
            print("time-varying r_equilibrium is not yet implemented")
            quit()

    print(r_eqm)

    cov_mat = calc_cov_matrix(md[settings.settings["equilibration_timesteps"]:],atoms,r_eqm)
    #print(f"covariance matrix calc")
    evecs = evals_evecs(cov_mat,atoms)
    if settings.settings["write_jmol"]:
        write_jmol_script(args.seed,evecs,atoms,((r_eqm*ureg.a_u_length).to("angstrom")).magnitude)
    #print("writing jmol")
    #if write_jmol:
    #    pos = get_atom_pos(args.seed)
    #    #print(pos)
    #    write_jmol_script(args.seed,evecs,atoms,pos)
    #print(md_dict.keys())
    #print(md_dict["ions"],len(md_dict["ions"]))