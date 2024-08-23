#!/usr/bin/env python

"""
Script: vasp2cif.py
Description:
    This script facilitates the conversion of atomic structure files from VASP POSCAR format to other formats 
    such as CIF, XYZ, or Quantum ESPRESSO input files. It interacts with the user to obtain the necessary 
    file name and desired output format, then performs the conversion using the ASE (Atomic Simulation Environment) library.

Key Features:
    - User-friendly: Guides the user through file selection and format conversion options.
    - Supports multiple output formats: CIF (.cif), XYZ (.xyz), and Quantum ESPRESSO (.in).
    - Automatically generates output files with appropriate extensions based on user choice.

Usage:
    - Run the script by executing `python vasp2cif.py` and follow the prompts.
    - The script will ask for the name of the VASP POSCAR file and the desired conversion format.
    - The converted file is saved in the same directory as the original file with the appropriate extension.

Dependencies:
    - ASE (Atomic Simulation Environment): Used for reading the POSCAR file and writing the converted output files.
    - OS: For handling file paths and checking file existence.

Functions:
    - get_user_input(): Interacts with the user to get the input file name and desired output format.
    - convert_files(structure, base_name, format_choice): Converts the structure to the selected format and saves it to a file.
    - main(): The main function that orchestrates user input, file reading, and file conversion.

Execution:
    - The script runs in a loop, ensuring valid user input for both the file name and conversion format.
    - It then reads the POSCAR file and performs the conversion based on the user's choice.
    - The converted file is saved in the current working directory.

Example:
    $ python vasp2cif.py
    Please input the name of the VASP POSCAR file (e.g., POSCAR): POSCAR
    Select the format to convert to:
    1. CIF (.cif)
    2. XYZ (.xyz)
    3. Quantum ESPRESSO (.in)
    Enter the number corresponding to your choice (1, 2, or 3): 1
    Structure written to 'POSCAR.cif' file.

Output:
    - The script outputs the converted file in the format chosen by the user, saved in the same directory as the input file.
"""

from ase.io import read, write
import os

def get_user_input():
    while True:
        file_name = input("Please input the name of the VASP POSCAR file (e.g., POSCAR): ").strip()
        file_path = os.path.join(os.getcwd(), file_name)
        if os.path.isfile(file_path):
            break
        print(f"File '{file_name}' not found in the current directory. Please provide a valid file name.")
    
    while True:
        print("Select the format to convert to:")
        print("1. CIF (.cif)")
        print("2. XYZ (.xyz)")
        print("3. Quantum ESPRESSO (.in)")
        format_choice = input("Enter the number corresponding to your choice (1, 2, or 3): ").strip()
        if format_choice in ['1', '2', '3']:
            break
        print("Invalid choice. Please enter 1, 2, or 3.")
    
    return file_name, file_path, format_choice

def convert_files(structure, base_name, format_choice):
    if format_choice == '1':
        output_file = f"{base_name}.cif"
        write(output_file, structure)
        print(f"Structure written to '{output_file}' file.")
    elif format_choice == '2':
        output_file = f"{base_name}.xyz"
        write(output_file, structure)
        print(f"Structure written to '{output_file}' file.")
    elif format_choice == '3':
        output_file = f"{base_name}.in"
        write(output_file, structure, format='espresso-in')
        print(f"Structure written to '{output_file}' file.")

def main():
    file_name, file_path, format_choice = get_user_input()
    structure = read(file_path, format='vasp')
    
    # Create the base name for output files
    base_name, _ = os.path.splitext(file_name)
    
    convert_files(structure, base_name, format_choice)

if __name__ == "__main__":
    main()

