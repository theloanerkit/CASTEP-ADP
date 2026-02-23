def no_lattice_cart():
    string = f"No lattice_cart block in cell file\n"
    string += "Reading atomic positions from a cell file requires lattice_cart block"
    print(string)
    quit()

def file_not_found(fname,ext=""):
    string = f"No file found with filename {fname}{ext}\n"
    string += "Please check the seed matches your files."
    print(string)
    quit()

def unexpected_format(line,fname):
    string = f"Unexpected format in {fname}.adp, line: {line}\n"
    string += "All lines should be of the form <parameter> : <value>"
    print(string)
    quit()

def unknown_keyword(key,fname):
    string = f"Unknown keyword {key} found in {fname}.adp"
    print(string)
    quit()

def invalid_parameter(key,value,optns=None,kind=None,lbound=None,ubound=None):
    string = ""
    if kind is None and lbound is None and ubound is None:
        string += f"Unexpected value {value} for {key}, expected "
        if optns is None:
            string += "true or false"
        else:
            substring = "/".join(optns)
            string += f"one of {substring}"
    elif lbound is not None:
        string += f"Unexpected value {value} for {key}, {key} >= {lbound}"
    elif ubound is not None:
        string += f"Unexpected value {value} for {key}, {key} <= {ubound}"
    else:
        n = ""
        if kind[0] in ["a","e","i","o","u"]: n = "n"
        string += f"{key} must be a{n} {kind}, got {value} instead"
    print(string)
    quit()