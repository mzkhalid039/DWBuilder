#!/usr/bin/env python

"""
Script: supercell.py
Description:
    This script generates a supercell from a given VASP POSCAR file. The user specifies the supercell dimensions 
    along the a, b, and c crystallographic directions. The resulting supercell structure is then saved to a new 
    VASP POSCAR file.

Key Features:
    - Prompts the user for input: file name of the VASP POSCAR file and supercell dimensions.
    - Validates the existence of the input file in the current directory.
    - Constructs the supercell using the specified dimensions.
    - Saves the generated supercell to a new POSCAR file in the current directory.

Usage:
    - Run the script by executing `python supercell.py` and follow the on-screen prompts.

Dependencies:
    - ASE (Atomic Simulation Environment): For handling atomic structures, reading/writing POSCAR files, and creating supercells.
    - Numpy: For matrix operations during supercell construction.
    - OS: For handling file paths and checking file existence.

Functions:
    - get_user_input(): Prompts the user for the input file name and supercell dimensions, and validates the file's existence.
    - build_supercell(structure, a, b, c): Constructs the supercell based on the user-provided dimensions.
    - main(): The main function that coordinates user input, supercell construction, and file output.

Output:
    - The script outputs the supercell structure to a new POSCAR file in the current directory, with a modified file name indicating the supercell.
"""

import numpy as np
from ase import io
from ase.build import make_supercell
from ase.io import read, write
import os

def get_user_input():
    while True:
        file_name = input("Please input the name of the VASP format file (e.g., POSCAR): ").strip()
        file_path = os.path.join(os.getcwd(), file_name)
        if os.path.isfile(file_path):
            break
        print(f"File '{file_name}' not found in the current directory. Please provide a valid file name.")
    
    a = int(input("Enter the supercell size along the a direction: "))
    b = int(input("Enter the supercell size along the b direction: "))
    c = int(input("Enter the supercell size along the c direction: "))
    
    return file_name, file_path, a, b, c

def build_supercell(structure, a, b, c):
    # Define the supercell matrix
    supercell_matrix = np.diag([a, b, c])
    
    # Build the supercell
    supercell = make_supercell(structure, supercell_matrix)
    
    return supercell

def main():
    file_name, file_path, a, b, c = get_user_input()
    structure = read(file_path, format='vasp')
    supercell = build_supercell(structure, a, b, c)
    
    # Create the output file name
    base_name, ext = os.path.splitext(file_name)
    output_file_name = f"{base_name}_supercell{ext}"
    
    # Output the supercell to a VASP POSCAR file
    write(output_file_name, supercell)
    print(f"Supercell structure written to '{output_file_name}' file.")

if __name__ == "__main__":
    main()
