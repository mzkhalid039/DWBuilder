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
