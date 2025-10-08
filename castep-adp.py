import argparse
import os
import pint
from copy import deepcopy as cp
import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as sc

ureg = pint.UnitRegistry()

parser = argparse.ArgumentParser()
parser.add_argument("-s","--seed")

class md_frame:
    def __init__(self):
        self.time = None
        self.positions = None
        self.energy = None

def parse_files(seed):
    print("parsing md file")
    with open(f"{seed}.md","r") as file:
        data = [line.strip() for line in file.readlines()]
    atoms = []
    get_atoms = True
    positions = []
    frames = []
    frame = None
    start = data.index("END header")
    data = data[start:]
    for i in range(len(data)):
        if len(data[i]) == 0:
            if frame is not None:
                get_atoms = False
                frame.positions = cp(positions)
                frames.append(frame)
                positions = []
            frame = md_frame()
        if "<-- R" in data[i]:
            line = data[i].split()
            if get_atoms:
                atoms.append("".join(line[0:2]))
            pos = list(map(float,line[2:5]))
            positions.append(pos)
        if "<-- E" in data[i]:
            frame.time = (float(data[i-1])*ureg.a_u_time).to("s").magnitude
            frame.energy = (float(data[i].split()[0])*ureg.a_u_energy).to("electron_volt").magnitude
    return frames, atoms

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
            file.write(f"ellipsoid ID {atoms[i]} AXES ")
            for j in range(len(axes[i])):
                file.write("{ ")
                for k in range(len(axes[i][j])):
                    file.write(f"{float(axes[i][j][k])} ")
                file.write("} ")
            file.write("center { ")
            for j in range(len(atom_pos[i])):
                file.write(f"{atom_pos[i][j]} ")
            file.write("} scale 5")

            file.write("\n")

if __name__ == "__main__":
    args = parser.parse_args()
    print(f"seed = {args.seed}")
    write_jmol = True
    if not os.path.isfile(f"{args.seed}.md"):
        print(f"no file found: [{args.seed}.md]")
        quit()
    if not os.path.isfile(f"{args.seed}.cell"):
        print(f"cell file not found - no jmol output")
        write_jmol = False
    equ_steps = int(input("number of equilibration steps: "))
    md,atoms = parse_files(args.seed)
    t = []
    e = []
    for frame in md:
        t.append(frame.time)
        e.append(frame.energy)
    plt.plot(t,e)
    plt.show()
    r_eqm = calc_r_eqm(md[equ_steps:],atoms)
    #print(f"r_eqm = {r_eqm}")
    cov_mat = calc_cov_matrix(md[equ_steps:],atoms,r_eqm)
    print(f"covariance matrix calc")
    evecs = evals_evecs(cov_mat,atoms)
    print("writing jmol")
    if write_jmol:
        pos = get_atom_pos(args.seed)
        #print(pos)
        write_jmol_script(args.seed,evecs,atoms,pos)
    #print(md_dict.keys())
    #print(md_dict["ions"],len(md_dict["ions"]))