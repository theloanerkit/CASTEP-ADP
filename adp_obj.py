import numpy as np
import adp_err

class MD:
    def __init__(self,timesteps,atoms,positions,velocities,ureg):
        self.timesteps = timesteps
        self.atoms = atoms
        self.positions = positions
        self.velocities = velocities
        self.apply_units(ureg)

    def __repr__(self):
        string = "MD object: \n"
        string += f"    number of timesteps: {self.timesteps}\n"
        string += f"    number of atoms: {len(self.atoms)}\n"
        string += f"    position units: {self.positions.units}\n"
        string += f"    velocity units: {self.velocities.units}"
        return string

    def apply_units(self,ureg):
        self.positions *= ureg.atomic_unit_of_length
        self.velocities *= (ureg.atomic_unit_of_length/ureg.atomic_unit_of_time)

class Cell:
    def __init__(self,lattice,positions,ureg):
        self.lattice = lattice
        self.positions = positions
        self.lattice_cart = None
        self.positions_abs = None
        self.positions_frac = None
        self.construct_cell(ureg)
        if self.positions_abs is None:
            self.get_positions_abs(ureg)

    def get_positions_abs(self,ureg):
        positions = {}
        for key in self.positions_frac.keys():
            pos = np.zeros(3,dtype=float)*ureg.angstrom
            for i in range(3):
                pos += self.positions_frac[key][i] * self.lattice_cart[i]
            positions[key] = pos
        self.positions_abs = positions

    def construct_cell(self,ureg):
        if "cart" in self.lattice[0].lower():
            self.construct_cart(ureg)
        elif "abc" in self.lattice[0].lower():
            adp_err.no_lattice_cart()
        if "abs" in self.positions[0].lower():
            self.construct_abs(ureg)
        elif "frac" in self.positions[0].lower():
            self.construct_frac()

    def construct_cart(self,ureg):
        start = 1
        if len(self.lattice) == 6:
            # units included, but assuming angstrom for now
            start += 1
        lattice = np.zeros((3,3),dtype=float)*ureg.angstrom
        for i in range(start,start+3):
            lattice[i-start] = np.asarray(self.lattice[i].split(),dtype=float)*ureg.angstrom
        self.lattice_cart = lattice

    def construct_abs(self,ureg):
        positions = {}
        elements = {}
        for line in self.positions[1:-1]:
            line = line.split()
            if len(line)==1:
                # units line, ignore (still assuming angstroms)
                continue
            if line[0] not in elements.keys():
                elements[line[0]] = 0
            elements[line[0]] += 1
            name = f"{line[0]} {elements[line[0]]}"
            pos = np.asarray(line[1:4],dtype=float)*ureg.angstrom
            positions[name] = pos
        self.positions_abs = positions

    def construct_frac(self):
        positions = {}
        elements = {}
        for line in self.positions[1:-1]:
            line = line.split()
            if len(line)==1:
                # units line, ignore (still assuming angstroms)
                continue
            if line[0] not in elements.keys():
                elements[line[0]] = 0
            elements[line[0]] += 1
            name = f"{line[0]} {elements[line[0]]}"
            pos = np.asarray(line[1:4],dtype=float)
            positions[name] = pos
        self.positions_frac = positions