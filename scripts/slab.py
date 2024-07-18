import os
import ase.io
from colorama import init, Fore, Style

init(autoreset=True)

class VacuumLayerManager:
    def __init__(self, filename):
        self.filename = filename
        self.structure = ase.io.read(filename)
        self.current_path = os.getcwd()
        self.log = []

    def prompt_vacuum_layer(self):
        vacuum_size = float(input("Enter the size of the vacuum layer in angstroms: "))
        direction = input("Enter the direction to add the vacuum layer (a, b, or c): ").lower()
        if direction not in ['a', 'b', 'c']:
            print("Invalid input. Defaulting to 'c'.")
            direction = 'c'
        self.log.append(f"Vacuum layer size: {vacuum_size} angstroms, direction: {direction}")
        return vacuum_size, direction

    def add_vacuum_layer(self, vacuum_size, direction):
        cell = self.structure.get_cell()
        new_cell = cell.copy()
        if direction == 'a':
            new_cell[0, 0] += vacuum_size
        elif direction == 'b':
            new_cell[1, 1] += vacuum_size
        elif direction == 'c':
            new_cell[2, 2] += vacuum_size
        self.structure.set_cell(new_cell, scale_atoms=False)
        
        vacuum_filename = os.path.join(self.current_path, 'structure_with_vacuum.vasp')
        self.structure.write(vacuum_filename, sort=True, vasp5=True)
        self.log.append(f"Structure with vacuum saved to {vacuum_filename}")
        print(f"Structure with vacuum saved to {vacuum_filename}")

    def run(self):
        vacuum_size, direction = self.prompt_vacuum_layer()
        self.add_vacuum_layer(vacuum_size, direction)

        with open("LOGFILE.txt", "w") as logfile:
            logfile.write("\n".join(self.log))
        print("Operation completed! Log file written to LOGFILE.txt")

if __name__ == "__main__":
    filename = input("Enter the input file name (with extension): ")
    manager = VacuumLayerManager(filename)
    manager.run()
