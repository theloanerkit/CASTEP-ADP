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

# for handling units
ureg = pint.UnitRegistry()

# for parsing command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("seed",
                    help="seedname for the files (.md/.phonon/.cell/.adp)")
parser.add_argument("-d","--dryrun",
                    action="store_true")

# set up dictionary with masses of elements
masses = {}
for key in adp_constants.masses.keys():
    masses[key] = adp_constants.masses[key] * ureg.unified_atomic_mass_unit

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
    ordered_pos = np.zeros((len(atoms),3),dtype=float)*ureg.angstrom
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



if __name__=="__main__":
    args = parser.parse_args()

    print("loading settings")
    settings = adp_settings.Settings(args.seed,ureg)

    data = {"atoms":None,
            "r_eq" :None,
            "adp"  :None,
            "uij"  :None,
            "ke"   :None}

    if not args.dryrun:
        if not os.path.isfile(f"{args.seed}.md"):
            adp_err.file_not_found(args.seed,".md")

        print("parsing md file")
        md_obj = adp_parse.parse_md(args.seed,settings.settings["equilibration_timesteps"],ureg)
        data["atoms"] = md_obj.atoms

        if settings.settings["r_equilibrium"] == "zero":
            print("reading r_eq from cell")
            if not os.path.isfile(f"{args.seed}.cell"):
                adp_err.file_not_found(args.seed,".cell")
            cell_obj = adp_parse.parse_cell(args.seed,ureg)
            data["r_eq"] = reorder_positions(cell_obj.positions_abs,data["atoms"])

        elif settings.settings["r_equilibrium"] == "finite":
            print("calculating r_eq")
            data["r_eq"] = calc_r_eq_from_md(md_obj)

        if settings.settings["calculate_adp"]:
            print("calculating covariance matrix")
            cov_mat = calc_covariance_matrix(md_obj,data["r_eq"])
            data["uij"] = cov_mat

            print("calculating adp axes")
            unit = cov_mat.units**0.5
            axes = evals_evecs(cov_mat,md_obj.atoms,sqrt=True)
            axes *= unit
            data["adp"] = axes

        if settings.settings["calculate_ke"]:
            print("calculating ke")
            data["ke"] = calc_ke_tensor(md_obj)

    print("writing")
    adp_io.write_out(args.seed,settings,data,args.dryrun)

    if settings.settings["write_jmol"]:
        print("writing jmol")
        adp_io.write_jmol_script(args.seed,
                                 data["adp"].to("angstrom").magnitude,
                                 data["atoms"],
                                 data["r_eq"].to("angstrom").magnitude,
                                 settings.settings["jmol_scale"])

