#!/usr/bin/env python
"""
 Overview
 --------
 This script automates the creation and manipulation of atomic interface structures
 between two different bulk phases. The script performs the following operations:
 - Reads atomic structure data from two input files representing bulk phase 1 and bulk phase 2.
 - Determines space group number, international symbol, and lattice type for each phase using pymatgen.
 - Allows users to define lattice directions for both phases to create and manipulate slabs using ASE.
 - Stacks the two phases to create an interface structure based on the user's chosen stacking direction.
 - Calculates and displays lattice and angular strain between the two phases.
 - Provides a warning about potential interface/domain wall artifacts that may require manual adjustment.

 Dependencies:
 - ase
 - numpy
 - colorama
 - pymatgen

 Usage:
 - Ensure the required dependencies are installed and the script is executed in a Python environment.
 - Follow the prompts to input file names, lattice directions, and other user inputs.
 - The script creates and saves manipulated structures in a "HIS" directory.
 - Review the warnings for potential artifacts in the generated interface/domain wall structures.
"""

import os
import shutil
import ase
import ase.io
from ase.build import cut, rotate, stack
import numpy as np
from colorama import init, Fore, Style
from numpy.linalg import norm
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from ase.neighborlist import NeighborList


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


def calculate_strain(slab1, slab2, log):
    a1, b1, c1 = slab1.cell
    a2, b2, c2 = slab2.cell

    strain_a = (norm(a1) - norm(a2)) / norm(a2) * 100
    strain_b = (norm(b1) - norm(b2)) / norm(b2) * 100
    strain_c = (norm(c1) - norm(c2)) / norm(c2) * 100

    theta_a = np.arccos(np.dot(a1, a2) / (norm(a1) * norm(a2)))
    theta_b = np.arccos(np.dot(b1, b2) / (norm(b1) * norm(b2)))
    theta_c = np.arccos(np.dot(c1, c2) / (norm(c1) * norm(c2)))

    log.append(f"Strain along a (%): {strain_a}")
    log.append(f"Strain along b (%): {strain_b}")
    log.append(f"Strain along c (%): {strain_c}")

    log.append(f"Angular strain along a (radians): {theta_a}")
    log.append(f"Angular strain along b (radians): {theta_b}")
    log.append(f"Angular strain along c (radians): {theta_c}")

    print('Strain along a (%):', strain_a)
    print('Strain along b (%):', strain_b)
    print('Strain along c (%):', strain_c)

    print('Angular strain along a (radians):', theta_a)
    print('Angular strain along b (radians):', theta_b)
    print('Angular strain along c (radians):', theta_c)


def main():
    init(autoreset=True)
    log = []

    filename1 = input("Enter the bulk phase 1 file name (with extension): ")
    filename2 = input("Enter the bulk phase 2 file name (with extension): ")
    log.append(f"Bulk phase 1 file: {filename1}")
    log.append(f"Bulk phase 2 file: {filename2}")

    D1 = read_structure(filename1)
    D2 = read_structure(filename2)

    print("Bulk phase 1:")
    get_symmetry_info(filename1, log)
    print("Bulk phase 2:")
    get_symmetry_info(filename2, log)

    current_path = os.getcwd()
    directory = get_unique_path(os.path.join(current_path, 'HIS'))
    log.append(f"Output directory: {directory}")

    a = get_lattice_direction(Fore.GREEN + Style.BRIGHT + "Enter three comma-separated values for bulk 1 lattice direction a: " + Style.RESET_ALL, log)
    b = get_lattice_direction(Fore.GREEN + Style.BRIGHT + "Enter three comma-separated values for bulk 1 lattice direction b: " + Style.RESET_ALL, log)
    c = get_lattice_direction(Fore.GREEN + Style.BRIGHT + "Enter three comma-separated values for bulk 1 lattice direction c: " + Style.RESET_ALL, log)

    a1 = get_lattice_direction(Fore.BLUE + Style.BRIGHT + "Enter three comma-separated values for bulk 2 lattice direction a: " + Style.RESET_ALL, log)
    b1 = get_lattice_direction(Fore.BLUE + Style.BRIGHT + "Enter three comma-separated values for bulk 2 lattice direction b: " + Style.RESET_ALL, log)
    c1 = get_lattice_direction(Fore.BLUE + Style.BRIGHT + "Enter three comma-separated values for bulk 2 lattice direction c: " + Style.RESET_ALL, log)

    slab1 = create_slab(D1, a, b, c)
    save_slab(slab1, directory, 'bulk1.vasp')

    slab2 = create_slab(D2, a1, b1, c1)
    save_slab(slab2, directory, 'bulk2.vasp')

    interface_direction = int(input("Enter the stacking direction (0 for a, 1 for b, and 2 for c): "))
    log.append(f"Stacking direction: {interface_direction}")
    slab = stack(slab1, slab2, axis=interface_direction, maxstrain=None)
    save_slab(slab, directory, 'interface.vasp')

    cutoff = float(input("Enter the cutoff distance in angstroms to remove close atoms: "))
    log.append(f"Cutoff distance: {cutoff}")
    slab = remove_close_atoms(slab, cutoff)
    save_slab(slab, directory, 'interface_cleaned.vasp')

    calculate_strain(slab1, slab2, log)

    print("\033[1;31;40mWarning: This code could generate interface/domain wall artifacts (i.e., Oxygen atoms, duplicate atoms etc,) at the interface, thus it requires manual adjustment.\033[0m")
    log.append("Warning: This code could generate interface/domain wall artifacts (i.e., Oxygen atoms, duplicate atoms etc,) at the interface, thus it requires manual adjustment.")

    with open(os.path.join(directory, "LOGFILE.txt"), "w") as logfile:
        logfile.write("\n".join(log))
    print(f"Log file written to {os.path.join(directory, 'LOGFILE.txt')}")


if __name__ == "__main__":
    main()
