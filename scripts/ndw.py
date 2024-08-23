#!/usr/bin/env python

"""
Script: ndw.py
Description:
    This script constructs a neutral domain wall (NDW) structure by stacking two crystal structures (P1 and P2) 
    along the b-axis in a supercell. The final structure is written to a VASP POSCAR file. The script also logs 
    the process in a logfile.

Usage:
    python ndw.py <supercell_size> <P1_filename> <P2_filename> <output_file>

    - <supercell_size>: The size of the supercell along the b-axis.
    - <P1_filename>: Path to the POSCAR file of the first crystal structure (P1).
    - <P2_filename>: Path to the POSCAR file of the second crystal structure (P2).
    - <output_file>: The name of the output POSCAR file containing the final NDW structure.

Functions:
    - main(supercell, P1_filename, P2_filename, output_file, log): Main function that drives the script.
    - import_poscar(file_path, log): Reads and parses a VASP POSCAR file, returning lattice vectors, atom types, atom counts, and atomic coordinates.
    - write_poscar(file_path, lattice_vectors, atom_types, atom_counts, coordinates, log): Writes the final supercell structure to a POSCAR file.

Logging:
    - The script logs key steps and any errors encountered to "LOGFILE.txt".
    - The command executed is also logged for reference.

Notes:
    - The script automatically adjusts atomic coordinates to ensure they are within the unit cell.
    - The supercell is constructed by stacking P1 in the first half and P2 in the second half along the b-axis.
"""


import numpy as np
import sys

def main(supercell, P1_filename, P2_filename, output_file, log):
    log.append(f"Running NDW with supercell size: {supercell}")
    print(f"Running NDW with supercell size: {supercell}")

    try:
        lattice_vectors_P1, atom_types_P1, atom_counts_P1, P1_coord = import_poscar(P1_filename, log)
        _, _, _, P2_coord = import_poscar(P2_filename, log)
    except FileNotFoundError as e:
        log.append(f"Error: {e}")
        print(f"Error: {e}")
        return

    P1_coord = np.where(P1_coord > 0.97, P1_coord - 1, P1_coord)
    P2_coord = np.where(P2_coord > 0.97, P2_coord - 1, P2_coord)

    A, B, O = atom_counts_P1

    P1_coord_A = P1_coord[:A]
    P1_coord_B = P1_coord[A:A+B]
    P1_coord_O = P1_coord[A+B:]
    P2_coord_A = P2_coord[:A]
    P2_coord_B = P2_coord[A:A+B]
    P2_coord_O = P2_coord[A+B:]

    coordinates_A = []
    coordinates_B = []
    coordinates_O = []
    for ii in range(1, supercell + 1):
        if ii <= supercell / 2:
            coordinates_A.append(P1_coord_A[:, :3] / [1, supercell, 1] + [0, (ii - 1) / supercell, 0])
            coordinates_B.append(P1_coord_B[:, :3] / [1, supercell, 1] + [0, (ii - 1) / supercell, 0])
            coordinates_O.append(P1_coord_O[:, :3] / [1, supercell, 1] + [0, (ii - 1) / supercell, 0])
        else:
            coordinates_A.append(P2_coord_A[:, :3] / [1, supercell, 1] + [0, (ii - 1) / supercell, 0])
            coordinates_B.append(P2_coord_B[:, :3] / [1, supercell, 1] + [0, (ii - 1) / supercell, 0])
            coordinates_O.append(P2_coord_O[:, :3] / [1, supercell, 1] + [0, (ii - 1) / supercell, 0])

    coordinates_A = np.vstack(coordinates_A)
    coordinates_B = np.vstack(coordinates_B)
    coordinates_O = np.vstack(coordinates_O)

    coordinates = np.vstack([coordinates_A, coordinates_B, coordinates_O])

    lattice_vectors = lattice_vectors_P1.copy()
    lattice_vectors[1] *= supercell

    write_poscar(output_file, lattice_vectors, atom_types_P1, [A*supercell, B*supercell, O*supercell], coordinates, log)

    log.append(f"{output_file} has been written.")
    print(f"{output_file} has been written.")

def import_poscar(file_path, log):
    log.append(f"Importing POSCAR from {file_path}")
    with open(file_path, 'r') as file:
        lines = file.readlines()

    lattice_vectors = np.array([list(map(float, line.split())) for line in lines[2:5]])
    atom_types = lines[5].split()
    atom_counts = list(map(int, lines[6].split()))
    coords = np.array([list(map(float, line.split())) for line in lines[8:8+sum(atom_counts)]])
    log.append(f"Imported POSCAR from {file_path}")
    return lattice_vectors, atom_types, atom_counts, coords

def write_poscar(file_path, lattice_vectors, atom_types, atom_counts, coordinates, log):
    log.append(f"Writing POSCAR to {file_path}")
    with open(file_path, 'w') as file:
        file.write("Generated Supercell\n")
        file.write("1.0\n")
        for vec in lattice_vectors:
            file.write(f" {vec[0]:20.10f} {vec[1]:20.10f} {vec[2]:20.10f}\n")
        file.write(" ".join(atom_types) + "\n")
        file.write(" ".join(map(str, atom_counts)) + "\n")
        file.write("Direct\n")
        for coord in coordinates:
            file.write(f" {coord[0]:20.10f} {coord[1]:20.10f} {coord[2]:20.10f}\n")
    log.append(f"POSCAR written to {file_path}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python ndw.py <supercell_size> <P1_filename> <P2_filename> <output_file>")
        sys.exit(1)
    supercell_size = int(sys.argv[1])
    P1_filename = sys.argv[2]
    P2_filename = sys.argv[3]
    output_file = sys.argv[4]
    
    log = []
    log.append(f"Command: python ndw.py {supercell_size} {P1_filename} {P2_filename} {output_file}")
    main(supercell_size, P1_filename, P2_filename, output_file, log)
    
    with open("LOGFILE.txt", "w") as logfile:
        logfile.write("\n".join(log))
