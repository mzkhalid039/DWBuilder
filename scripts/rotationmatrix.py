#!/usr/bin/env python
"""
 Overview
 --------
 This script automates the creation and manipulation of atomic structures
 of a bulk phase. The script performs the following operations:
 - Reads atomic structure data from an input file representing a bulk phase.
 - Determines space group number, international symbol, and lattice type for the phase using pymatgen.
 - Allows users to define lattice directions for the phase to create and manipulate slabs using ASE.

"""

import os
import ase
import ase.io
from ase.build import cut, rotate
import numpy as np
from colorama import init, Fore, Style
from numpy.linalg import norm
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from ase.neighborlist import NeighborList

init(autoreset=True)

def read_structure(file_name):
    return ase.io.read(file_name)

def get_symmetry_info(file_name, log):
    struct = Structure.from_file(file_name)
    sym = SpacegroupAnalyzer(struct)
    data = sym.get_symmetry_dataset()
    log.append(f"Space group number: {data['number']}")
    log.append(f"International symbol: {data['international']}")
    log.append(f"Lattice type: {sym.get_lattice_type()}")
    print(f"Space group number: {data['number']}")
    print(f"International symbol: {data['international']}")
    print(f"Lattice type: {sym.get_lattice_type()}")

def get_lattice_direction(prompt, log):
    direction_input = input(prompt)
    log.append(f"{prompt.strip()}: {direction_input}")
    return [float(value.strip()) for value in direction_input.split(',')]

def create_slab(D, a, b, c):
    slab = cut(D, a=a, b=b, c=c)
    rotate(slab, slab.cell[0], (0, 1, 0), slab.cell[1], (1, 0, 0))
    return slab

def get_unique_path(directory):
    if not os.path.exists(directory):
        return directory
    base, extension = os.path.splitext(directory)
    counter = 1
    while os.path.exists(directory):
        directory = f"{base}_{counter}"
        counter += 1
    return directory

def save_slab(slab, directory, file_name):
    os.makedirs(directory, exist_ok=True)
    slab.write(os.path.join(directory, file_name), sort=True, vasp5=True)

def remove_close_atoms(atoms, cutoff):
    nl = NeighborList([cutoff / 2.0] * len(atoms), self_interaction=False, bothways=True)
    nl.update(atoms)
    indices_to_remove = set()
    for i in range(len(atoms)):
        indices, offsets = nl.get_neighbors(i)
        for idx in indices:
            if i < idx:  # To avoid double counting
                distance = atoms.get_distance(i, idx)
                if distance < cutoff:
                    indices_to_remove.add(idx)
    atoms = atoms[[atom.index not in indices_to_remove for atom in atoms]]
    return atoms

def calculate_strain(slab, log):
    a, b, c = slab.cell

    strain_a = (norm(a) - norm(a)) / norm(a) * 100
    strain_b = (norm(b) - norm(b)) / norm(b) * 100
    strain_c = (norm(c) - norm(c)) / norm(c) * 100

    log.append(f"Strain along a (%): {strain_a}")
    log.append(f"Strain along b (%): {strain_b}")
    log.append(f"Strain along c (%): {strain_c}")

    print('Strain along a (%):', strain_a)
    print('Strain along b (%):', strain_b)
    print('Strain along c (%):', strain_c)

def main():
    init(autoreset=True)
    log = []

    filename = input("Enter the bulk phase file name (with extension): ")
    log.append(f"Bulk phase file: {filename}")

    D = read_structure(filename)

    print("Bulk phase:")
    get_symmetry_info(filename, log)

    current_path = os.getcwd()
    directory = get_unique_path(os.path.join(current_path, 'RotationMatrix'))
    log.append(f"Output directory: {directory}")

    a = get_lattice_direction(Fore.GREEN + Style.BRIGHT + "Enter three comma-separated values for lattice direction a: " + Style.RESET_ALL, log)
    b = get_lattice_direction(Fore.GREEN + Style.BRIGHT + "Enter three comma-separated values for lattice direction b: " + Style.RESET_ALL, log)
    c = get_lattice_direction(Fore.GREEN + Style.BRIGHT + "Enter three comma-separated values for lattice direction c: " + Style.RESET_ALL, log)

    slab = create_slab(D, a, b, c)
    save_slab(slab, directory, 'bulk.vasp')

    cutoff = float(input("Enter the cutoff distance in angstroms to remove close atoms: "))
    log.append(f"Cutoff distance: {cutoff}")
    slab = remove_close_atoms(slab, cutoff)
    save_slab(slab, directory, 'bulk_cleaned.vasp')

    calculate_strain(slab, log)


    with open(os.path.join(directory, "LOGFILE.txt"), "w") as logfile:
        logfile.write("\n".join(log))
    print(f"Log file written to {os.path.join(directory, 'LOGFILE.txt')}")

if __name__ == "__main__":
    main()
