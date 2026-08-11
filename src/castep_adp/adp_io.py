from . import adp_constants
import numpy as np

col = 12
col_num = 24
col_long = 32

def extend_str(string,width):
    string = " "+string
    string = string + " "*(width-len(string))
    return string

def vec_to_string(vec):
    string = []
    for v in vec:
        try:
            num = v.magnitude
        except AttributeError:
            num = v
        num = f"{num}"
        if num[0] != "-":
            num = " "+num
        string.append(num)
    return string

def write_columns(cols,width):
    string = ""
    for c in cols:
        string += extend_str(c,width)
    return string

def write_jmol_script(seed,axes,atoms,atom_pos,scale,environments=None):
    env_dict = {}
    other_dict = {}
    env = ""
    if environments is not None:
        env = "_env"
        for i in range(len(environments)):
            for elem in environments[i]:
                env_dict[elem.name] = i
        idx = 0
        for j in range(len(atoms)):
            if atoms[j] not in env_dict.keys() and atoms[j].split()[0] not in other_dict.keys():
                other_dict[atoms[j].split()[0]] = idx
                idx += 1
    with open(f"{seed}_ellipsoid{env}.spt","w") as file:
        for i in range(len(axes)):
            tmp = atoms[i].replace(" ","")
            file.write(f"ellipsoid ID {tmp} AXES ")
            for j in range(len(axes[i])):
                file.write("{ ")
                for k in range(len(axes[i][j])):
                    file.write(f"{float(axes[i][j][k])} ")
                file.write("} ")
            file.write("center { ")
            for j in range(len(atom_pos[i])):
                file.write(f"{atom_pos[i][j]} ")
            file.write(f"}} scale {scale} ")
            if environments is None:
                col = adp_constants.jmol_colours[atoms[i].split()[0]]
            else:
                if atoms[i].split()[0] in other_dict.keys():
                    col = adp_constants.jmol_colours[atoms[i].split()[0]]
                else:
                    col = adp_constants.environment_colours[env_dict[atoms[i]]]
            file.write(f"color [x{col}]")
            file.write("\n")

def write_header(file,keys,settings):
    length_unit = settings.output_units["length"][1]
    energy_unit = settings.output_units["energy"][1]
    file.write("begin HEADER\n\n")
    file.write(f"version={adp_constants.VERSION}\n\n")
    file.write("--- Parameters ---\n")
    for key in keys:
        if settings.user_def[key]:
            file.write(extend_str("<userdef>",col))
        else:
            file.write(extend_str("<default>",col))
        file.write(extend_str(key,col_long))
        file.write(f"{settings.settings[key]}\n")
    file.write("\n")
    file.write("--- Outputs ---\n")
    # write out r_eq
    string = extend_str("<r_eq>",col)
    string += write_columns(["x","y","z"],col_num)
    string += length_unit
    file.write(f"{string}\n\n")

    if settings.settings["write_adp"]:
        string = extend_str("<adp>",col)
        string += write_columns(["x","y","z"],col_num)
        string += length_unit
        file.write(f"{string}\n\n")

    if settings.settings["write_uij"]:
        string = extend_str("<uij>",col)
        string1 = write_columns(["xx","xy","xz"],col_num)
        string2 = write_columns(["yx","yy","yz"],col_num)
        string3 = write_columns(["zx","zy","zz"],col_num)
        unit = length_unit
        file.write(f"{string}{string1}{unit}^2\n")
        file.write(f"{string}{string2}{unit}^2\n")
        file.write(f"{string}{string3}{unit}^2\n")

    if settings.settings["write_ke"]:
        string = extend_str("<ke>",col)
        string1 = write_columns(["xx","xy","xz"],col_num)
        string2 = write_columns(["yx","yy","yz"],col_num)
        string3 = write_columns(["zx","zy","zz"],col_num)
        unit = energy_unit
        file.write(f"{string}{string1}{unit}\n")
        file.write(f"{string}{string2}{unit}\n")
        file.write(f"{string}{string3}{unit}\n")

    file.write("\n")
    file.write("end HEADER\n\n")

def write_summary(file,settings,atoms):
    energy_unit = settings.output_units["energy"][0]
    elements = set()
    keys = list(atoms.keys())
    for atom in keys:
        elements.add(atom.split()[0])
    elements = list(elements)
    #print(elements)
    file.write("begin SUMMARY\n\n")
    for elem in elements:
        string = ""
        string += write_columns(["<atom>",elem],col)
        string += "\n"
        file.write(string)
        if settings.settings["write_ke"]:
            ke_tensor = np.zeros((3,3),dtype=float)*atoms[keys[0]].ke.units
            for atom in keys:
                if atom.split()[0] == elem:
                    ke_tensor += atoms[atom].ke
            for vec in ke_tensor:
                string = extend_str("<ke>",col)
                #unit = settings.output_unit["energy"][0]
                str_vec = vec_to_string(vec.to(energy_unit))
                string += write_columns(str_vec,col_num)
                file.write(f"{string}\n")
        file.write("\n")
    file.write("end SUMMARY\n\n")

def write_out(seed,settings,atoms,dryrun):
    length_unit = settings.output_units["length"][0]
    energy_unit = settings.output_units["energy"][0]
    keys = list(settings.settings.keys())
    keys.sort()
    with open(f"{seed}.out","w") as file:
        write_header(file,keys,settings)
        if not dryrun:
            if settings.settings["element_summary"]:
                write_summary(file,settings,atoms)
            for atom in atoms.keys():
                string = extend_str("<atom>",col)
                string += extend_str(atom,col)
                file.write(f"{string}\n")

                if atoms[atom].r_eq is not None:
                    string = extend_str("<r_eq>",col)

                    r_eq = vec_to_string(atoms[atom].r_eq.to(length_unit))
                    string += write_columns(r_eq,col_num)
                    file.write(f"{string}\n")

                if settings.settings["write_adp"]:
                    for vec in atoms[atom].adp:
                        string = extend_str("<adp>",col)
                        str_vec = vec_to_string(vec.to(length_unit))
                        string += write_columns(str_vec,col_num)
                        file.write(f"{string}\n")

                if settings.settings["write_uij"]:
                    for vec in atoms[atom].uij:
                        string = extend_str("<uij>",col)
                        str_vec = vec_to_string(vec.to(length_unit**2))
                        string += write_columns(str_vec,col_num)
                        file.write(f"{string}\n")

                if settings.settings["write_ke"]:
                    for vec in atoms[atom].ke:
                        string = extend_str("<ke>",col)
                        str_vec = vec_to_string(vec.to(energy_unit))
                        string += write_columns(str_vec,col_num)
                        file.write(f"{string}\n")

                file.write("\n")
