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
