import re
import numpy as np
from . import obj
from . import err

def parse_md(seed,equ_timesteps=0):
  # load md file
  with open(f"{seed}.md") as file:
    md_file = [line.strip() for line in file.readlines()]

  # get the start index of each of the timestep blocks
  # searches for "<-- E" and subtracts 1 to get starting index
  idxs = np.asarray(
    [i for i in range(len(md_file)) if re.match(".*<-- E$",md_file[i]) is not None]
  )-1

  # create list of atom names
  offset = 0
  in_block = False
  atoms = []
  start,stop = idxs[0],idxs[1]      # +6 gets to the start of <-- R block
  for line in md_file[start:stop]:
    if "<-- R" not in line and in_block:  # finished looking at positions, so we have all the atoms
      break
    elif "<-- R" not in line and not in_block:
      offset += 1
    elif "<-- R" in line:
      in_block = True
      line = line.split()
      atoms.append(f"{line[0]} {line[1]}")    # add atom name and number to list

  idxs = idxs[equ_timesteps+1:]       # indexs that we care about (first block is at t=0)

  if len(idxs) == 0:
    raise err.NoTimesteps(equ_timesteps)

  # set up positions array
  positions = np.zeros((len(idxs),len(atoms),3),dtype=float)
  for i in range(len(idxs)):
    start = idxs[i] + offset
    stop = start + len(atoms)
    for j in range(start,stop):
      positions[i,j-start,:] = np.asarray(md_file[j].split()[2:5],dtype=float)

  # set up velocities array
  velocities = np.zeros((len(idxs),len(atoms),3),dtype=float)
  for i in range(len(idxs)):
    start = idxs[i] + offset + len(atoms)
    stop = start + len(atoms)
    for j in range(start,stop):
      velocities[i,j-start,:] = np.asarray(md_file[j].split()[2:5],dtype=float)

  # create MD object
  md = obj.MD(len(idxs),atoms,positions,velocities)

  return md

def parse_cell(seed):
  with open(f"{seed}.cell") as file:
    cell_file = [line.strip() for line in file.readlines()]

  lattice_block = []
  positions_block = []

  l_block = False
  p_block = False

  for line in cell_file:
    if "%block" in line.lower():
      if "positions" in line.lower():
        p_block = True
      elif "lattice" in line.lower():
        l_block = True
    if l_block:
      lattice_block.append(line)
    elif p_block:
      positions_block.append(line)
    if "%endblock" in line.lower():
      if "positions" in line.lower():
        p_block = False
      elif "lattice" in line.lower():
        l_block = False

  cell = obj.Cell(lattice_block,positions_block)

  return cell
