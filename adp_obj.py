import numpy as np
import adp_err
import adp_constants

class MD:
    def __init__(self,timesteps,atoms,positions,velocities):
        self.timesteps = timesteps
        self.atoms = atoms
        self.positions = positions
        self.velocities = velocities
        self.apply_units()

    def __repr__(self):
        string = "MD object: \n"
        string += f"    number of timesteps: {self.timesteps}\n"
        string += f"    number of atoms: {len(self.atoms)}\n"
        string += f"    position units: {self.positions.units}\n"
        string += f"    velocity units: {self.velocities.units}"
        return string

    def apply_units(self):
        self.positions *= adp_constants.ureg.atomic_unit_of_length
        self.velocities *= (adp_constants.ureg.atomic_unit_of_length/adp_constants.ureg.atomic_unit_of_time)

class Cell:
    def __init__(self,lattice,positions):
        self.lattice = lattice
        self.positions = positions
        self.lattice_cart = None
        self.positions_abs = None
        self.positions_frac = None
        self.construct_cell()
        if self.positions_abs is None:
            self.get_positions_abs()

    def get_positions_abs(self):
        positions = {}
        for key in self.positions_frac.keys():
            pos = np.zeros(3,dtype=float)*adp_constants.ureg.angstrom
            for i in range(3):
                pos += self.positions_frac[key][i] * self.lattice_cart[i]
            positions[key] = pos
        self.positions_abs = positions

    def construct_cell(self):
        if "cart" in self.lattice[0].lower():
            self.construct_cart()
        elif "abc" in self.lattice[0].lower():
            adp_err.no_lattice_cart()
        if "abs" in self.positions[0].lower():
            self.construct_abs()
        elif "frac" in self.positions[0].lower():
            self.construct_frac()

    def construct_cart(self):
        start = 1
        if len(self.lattice) == 6:
            # units included, but assuming angstrom for now
            start += 1
        lattice = np.zeros((3,3),dtype=float)*adp_constants.ureg.angstrom
        for i in range(start,start+3):
            lattice[i-start] = np.asarray(self.lattice[i].split(),dtype=float)*adp_constants.ureg.angstrom
        self.lattice_cart = lattice

    def construct_abs(self):
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
            pos = np.asarray(line[1:4],dtype=float)*adp_constants.ureg.angstrom
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