import os
import ase.io
from ase.build import cut, rotate
from ase.neighborlist import NeighborList
from numpy.linalg import norm
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from enum import Enum
from colorama import init, Fore, Style

init(autoreset=True)

class DomainWallAngle(Enum):
    R180 = 'R180'
    R71 = 'R71'
    R109 = 'R109'
    T180 = 'T180'
    T90 = 'T90'
    FDW = 'FDW'
    O120_HH_TT = '1'
    O120_HT = '2'
    O180 = '3'
    O90 = '4'
    ALL = 'ALL'

def get_angles_dict(domain_size):
    return {
        "R3m": {
            DomainWallAngle.R180: [
                {'a': [1, 1, 0], 'b': [0, 0, 1], 'c': [domain_size, -domain_size, 0], 'stack_axis': 2},
                {'a': [-1, -1, 0], 'b': [0, 0, -1], 'c': [domain_size, -domain_size, 0], 'stack_axis': 2}
            ],
            DomainWallAngle.R71: [
                {'a': [domain_size, domain_size, 0], 'b': [0, 0, 1], 'c': [1, -1, 0], 'stack_axis': 0},
                {'a': [domain_size, domain_size, 0], 'b': [0, 0, -1], 'c': [-1, 1, 0], 'stack_axis': 0}
            ],
            DomainWallAngle.R109: [
                {'a': [domain_size, 0.0, 0], 'b': [0, 1, 0], 'c': [0, 0, 1], 'stack_axis': 0},
                {'a': [domain_size, 0.0, 0], 'b': [0, -1, 0], 'c': [0, 0, -1], 'stack_axis': 0}
            ]
        },
        "R3c": {
            DomainWallAngle.R109: [
                {'a': [0.0, 1, 0], 'b': [1, 0, 0], 'c': [0, 0, domain_size], 'stack_axis': 2},
                {'a': [0.0, -1, 0], 'b': [-1, 0, 0], 'c': [0, 0, domain_size], 'stack_axis': 2}
            ],
            DomainWallAngle.R71: [
                {'a': [1, -1, 0], 'b': [0, 0, 1], 'c': [domain_size, domain_size, 0], 'stack_axis': 2},
                {'a': [-1, 1, 0], 'b': [0, 0, -1], 'c': [domain_size, domain_size, 0], 'stack_axis': 2}
            ],
            DomainWallAngle.R180: [
                {'a': [1, 1, 0], 'b': [0, 0, -1], 'c': [-domain_size, domain_size, 0], 'stack_axis': 2},
                {'a': [-1, -1, 0], 'b': [0, 0, 1], 'c': [-domain_size, domain_size, 0], 'stack_axis': 2}
            ]
        },
        "P4mm": {
            DomainWallAngle.T180: [
                {'a': [1.01, 0, 0], 'b': [0, domain_size, 0], 'c': [0, 0, 1.01], 'stack_axis': 1},
                {'a': [-1.01, 0, 0], 'b': [0, domain_size, 0], 'c': [0, 0, -1.01], 'stack_axis': 1}
            ],
            DomainWallAngle.T90: [
                {'a': [0, 1, 0], 'b': [-1, 0, 1], 'c': [domain_size, 0, domain_size], 'stack_axis': 2},
                {'a': [0, -1, 0], 'b': [1, 0, -1], 'c': [domain_size, 0, domain_size], 'stack_axis': 2}
            ]
        },
        "Pnma": {
            DomainWallAngle.FDW: [
                {'a': [1.0, 1.0, 0], 'b': [0.0, 0.0, 1], 'c': [domain_size, -domain_size, 0], 'stack_axis': 2},
                {'a': [-1, 1, 0], 'b': [0.0, 0.0, 1], 'c': [domain_size, domain_size, 0], 'stack_axis': 2}
            ]
        },
        "Pmc2_1": {
            DomainWallAngle.O120_HH_TT: [
                {'a': [-1, -2, 0], 'b': [2.0*domain_size, -1.0*domain_size, 0], 'c': [0, 0, 1], 'stack_axis': 1},
                {'a': [-1, 2, 0], 'b': [-2.0*domain_size, -1.0*domain_size, 0], 'c': [0, 0, 1], 'stack_axis': 1}
            ],
            DomainWallAngle.O120_HT: [
                {'a': [-1*domain_size, -2*domain_size, 0], 'b': [2.0, -1.0, 0], 'c': [0, 0, 1], 'stack_axis': 0},
                {'a': [-1*domain_size, 2*domain_size, 0], 'b': [-2.0, -1.0, 0], 'c': [0, 0, 1], 'stack_axis': 0}
            ],
            DomainWallAngle.O180: [
                {'a': [1, 0, 0], 'b': [0, domain_size, 0], 'c': [0, 0, 1], 'stack_axis': 1},
                {'a': [-1, 0, 0], 'b': [0, -domain_size, 0], 'c': [0, 0, 1], 'stack_axis': 1}
            ],
            DomainWallAngle.O90: [
                {'a': [0, 1*domain_size, 0], 'b': [-1, 0, 0], 'c': [0, 0, 1], 'stack_axis': 0},
                {'a': [-1*domain_size, 0, 0], 'b': [0, -1, 0], 'c': [0, 0, 1], 'stack_axis': 0}
            ]
        }
    }

def print_symmetry_info(sym, log):
    data = sym.get_symmetry_dataset()
    log.append(f"Space group number: {data['number']}")
    log.append(f"International symbol: {data['international']}")
    log.append(f"Lattice type: {sym.get_lattice_type()}")
    print(f"Space group number: {data['number']}")
    print(f"International symbol: {data['international']}")
    print(f"Lattice type: {sym.get_lattice_type()}")

def cut_and_rotate(D1, a, b, c):
    slab = cut(D1, a=a, b=b, c=c)
    rotate(slab, slab.cell[0], (0, 1, 0), slab.cell[1], (1, 0, 0))
    return slab

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

class DomainWallManager:
    def __init__(self, filename):
        self.filename = filename
        self.D1 = ase.io.read(filename)
        self.struct = Structure.from_file(filename)
        self.sym = SpacegroupAnalyzer(self.struct)
        self.current_path = os.getcwd()
        self.log = []

    def get_system_selection(self):
        print("Select the system:")
        print("1 - R3c")
        print("2 - R3m")
        print("3 - P4mm")
        print("4 - Pmc2_1")
        print("5 - Pnma")
        print("6 - P6_3cm")
        print("7 - Manual Orientation relationships")
        return int(input())

    def prompt_lattice_directions(self):
        a_input = input(Fore.GREEN + Style.BRIGHT + "Enter three comma-separated values for lattice direction a: " + Style.RESET_ALL)
        a = [float(value.strip()) for value in a_input.split(',')]
        self.log.append(f"Lattice direction a: {a}")

        b_input = input(Fore.GREEN + Style.BRIGHT + "Enter three comma-separated values for lattice direction b: " + Style.RESET_ALL)
        b = [float(value.strip()) for value in b_input.split(',')]
        self.log.append(f"Lattice direction b: {b}")

        c_input = input(Fore.GREEN + Style.BRIGHT + "Enter three comma-separated values for lattice direction c: " + Style.RESET_ALL)
        c = [float(value.strip()) for value in c_input.split(',')]
        self.log.append(f"Lattice direction c: {c}")

        return a, b, c

    def prompt_domain_wall_angle(self):
        valid_angles = self.get_valid_angles()
        print("Select the domain wall angle:")
        for angle in valid_angles:
            print(f"{angle.value} - {angle.name.replace('_', ' ')}")
        return DomainWallAngle(input())

    def get_valid_angles(self):
        if self.system_selection == 1:  # R3c
            return [DomainWallAngle.R180, DomainWallAngle.R71, DomainWallAngle.R109, DomainWallAngle.ALL]
        elif self.system_selection == 2:  # R3m
            return [DomainWallAngle.R180, DomainWallAngle.R71, DomainWallAngle.R109, DomainWallAngle.ALL]
        elif self.system_selection == 3:  # P4mm
            return [DomainWallAngle.T180, DomainWallAngle.T90, DomainWallAngle.ALL]
        elif self.system_selection == 4:  # Pmc2_1
            return [DomainWallAngle.O120_HH_TT, DomainWallAngle.O120_HT, DomainWallAngle.O180, DomainWallAngle.O90, DomainWallAngle.ALL]
        elif self.system_selection == 5:  # Pnma
            return [DomainWallAngle.FDW]
        elif self.system_selection == 6:  # P6_3cm
            return ["NDW", "CDW"]
        else:
            return list(DomainWallAngle)

    def convert_structure(self):
        cif_filename = os.path.join(self.current_path, 'structure.cif')
        xyz_filename = os.path.join(self.current_path, 'structure.xyz')

        self.D1.write(cif_filename)
        self.D1.write(xyz_filename)

        self.log.append(f"Structure converted to CIF format: {cif_filename}")
        self.log.append(f"Structure converted to XYZ format: {xyz_filename}")

    def process_slab(self, angles, angle_type, space_group, cutoff, domain_size, supercell_size=None):
        directory_name = os.path.join(self.current_path, f"{space_group}_{angle_type}")
        os.makedirs(directory_name, exist_ok=True)

        for i, angle in enumerate(angles):
            slab = cut_and_rotate(self.D1, **{k: v for k, v in angle.items() if k in ('a', 'b', 'c')})
            slab = remove_close_atoms(slab, cutoff)
            slab_filename = os.path.join(directory_name, f"{angle_type}_domain{i+1}.vasp")
            slab.write(slab_filename, sort=True, vasp5=True)
            self.log.append(f"Domain slab {i+1} saved to {slab_filename}")

            if supercell_size:
                supercell = slab.repeat(supercell_size)
                supercell = remove_close_atoms(supercell, cutoff)
                supercell_filename = os.path.join(directory_name, f"{angle_type}_domain{i+1}_supercell.vasp")
                supercell.write(supercell_filename, sort=True, vasp5=True)
                self.log.append(f"Supercell domain slab {i+1} saved to {supercell_filename}")

    def run(self):
        print_symmetry_info(self.sym, self.log)
        self.system_selection = self.get_system_selection()
        self.log.append(f"System selection: {self.system_selection}")

        if self.system_selection == 7:
            print("Manual Orientation relationships selected.")
            self.create_slab()
        else:
            space_groups = ["R3c", "R3m", "P4mm", "Pmc2_1", "Pnma", "P6_3cm"]
            space_group = space_groups[self.system_selection - 1]
            angle_selection = self.prompt_domain_wall_angle()
            domain_size = int(input("Enter the domain wall size (in number of unit cells): "))
            cutoff = float(input("Enter the cutoff distance in angstroms to remove close atoms: "))
            supercell_size = None
            if input("Do you want to create a supercell? (y/n): ").lower() == 'y':
                supercell_size = [int(x) for x in input("Enter the supercell size along a, b, and c directions (comma-separated): ").split(',')]
            angles_dict = get_angles_dict(domain_size)[space_group]
            if angle_selection == DomainWallAngle.ALL:
                for angle_type, angles in angles_dict.items():
                    self.process_slab(angles, angle_type.value, space_group, cutoff, domain_size, supercell_size)
            else:
                self.process_slab(angles_dict[angle_selection], angle_selection.value, space_group, cutoff, domain_size, supercell_size)

        if input("Do you want to convert the structure to CIF and XYZ formats? (y/n): ").lower() == 'y':
            self.convert_structure()

        with open("LOGFILE.txt", "w") as logfile:
            logfile.write("\n".join(self.log))
        print("Operation completed! Log file written to LOGFILE.txt")

if __name__ == "__main__":
    filename = input("Enter the input file name (with extension): ")
    manager = DomainWallManager(filename)
    manager.run()
