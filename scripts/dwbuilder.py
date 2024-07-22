#!/usr/bin/env python
"""
 Overview
 --------
 This script processes a given input file representing a crystal structure 
 and manipulates it to produce various domain wall structures based on 
 different space group symmetries. 

 The script uses `ase` and `pymatgen` libraries to manipulate atomic structures,
 determine symmetry properties, and generate domain wall configurations. It allows
 the user to choose from a range of predefined space groups and domain wall angles.
 
 The main features include:
 - Automatic space group identification and lattice type determination using `pymatgen`.
 - User input for selecting predefined systems and domain wall angles.
 - Creation of various domain wall structures (109°, 71°, 180°) and supercells.
 - Calculations of lattice strain between different domains.
 - Warnings about potential domain wall artifacts that may arise during the 
   process, advising manual adjustment.
 
 Usage:
 - Ensure required libraries (`ase`, `pymatgen`, `numpy`, and `colorama`) are 
   installed.
 - Run the script and follow prompts for user input.
 - The script generates domain wall structures and supercells in the respective 
   directories for the specified system and angle selections.

 Dependencies:
 - ase
 - pymatgen
 - numpy
 - colorama

 For detailed usage instructions, refer to the prompts in the script.
"""

import os
import ase.io
from ase.build import cut, rotate, stack
from ase.neighborlist import NeighborList
from numpy.linalg import norm
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from enum import Enum
from colorama import init, Fore, Style
import numpy as np
import subprocess

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
                {'a': [0, -1*domain_size, 0], 'b': [0, -1, 0], 'c': [0, 0, 1], 'stack_axis': 0}
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

def calculate_lattice_strain(slab1, slab2, angle_type, log):
    a1, b1, c1 = slab1.cell
    a2, b2, c2 = slab2.cell
    strain_a = (norm(a1) - norm(a2)) / norm(a2) * 100
    strain_b = (norm(b1) - norm(b2)) / norm(b2) * 100
    strain_c = (norm(c1) - norm(c2)) / norm(c2) * 100
    log.append(f"Lattice strain for {angle_type} DW:")
    log.append(f"strain along a (%): {strain_a:.2f}")
    log.append(f"strain along b (%): {strain_b:.2f}")
    log.append(f"strain along c (%): {strain_c:.2f}")
    print(f"Lattice strain for {angle_type} DW:")
    print(f"strain along a (%): {strain_a:.2f}")
    print(f"strain along b (%): {strain_b:.2f}")
    print(f"strain along c (%): {strain_c:.2f}")

def cut_and_rotate(D1, a, b, c):
    slab = cut(D1, a=a, b=b, c=c)
    rotate(slab, slab.cell[0], (0, 1, 0), slab.cell[1], (1, 0, 0))
    return slab

def get_unique_path(directory, filename):
    base, extension = os.path.splitext(filename)
    counter = 1
    unique_filename = filename
    while os.path.exists(os.path.join(directory, unique_filename)):
        unique_filename = f"{base}_{counter}{extension}"
        counter += 1
    return unique_filename

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

def process_slab(D1, angles, angle_type, space_group, cutoff, log, supercell_size=None):
    directory_name = get_unique_path('', f"{space_group}_{angle_type}")
    os.makedirs(directory_name, exist_ok=True)

    stack_axis = angles[0]['stack_axis']
    slabs = [cut_and_rotate(D1, **{k: v for k, v in angle.items() if k in ('a', 'b', 'c')}) for angle in angles]
    stacked_slab = stack(*slabs, axis=stack_axis, maxstrain=None)
    stacked_slab = remove_close_atoms(stacked_slab, cutoff)

    stacked_filename = get_unique_path(directory_name, f"{angle_type}_stacked.vasp")
    stacked_slab.write(os.path.join(directory_name, stacked_filename), sort=True, vasp5=True)
    log.append(f"Stacked slab saved to {stacked_filename}")

    for i, slab in enumerate(slabs):
        slab = remove_close_atoms(slab, cutoff)
        slab_filename = get_unique_path(directory_name, f"{angle_type}_domain{i+1}.vasp")
        slab.write(os.path.join(directory_name, slab_filename), sort=True, vasp5=True)
        log.append(f"Domain slab {i+1} saved to {slab_filename}")

    calculate_lattice_strain(*slabs, angle_type, log)

    if supercell_size:
        supercell = stacked_slab.repeat(supercell_size)
        supercell = remove_close_atoms(supercell, cutoff)
        supercell_filename = get_unique_path(directory_name, f"{angle_type}_supercell.vasp")
        supercell.write(os.path.join(directory_name, supercell_filename), sort=True, vasp5=True)
        log.append(f"Supercell saved to {supercell_filename}")

class DomainWallManager:
    def __init__(self, filename):
        self.filename = filename
        self.D1 = ase.io.read(filename)
        self.struct = Structure.from_file(filename)
        self.sym = SpacegroupAnalyzer(self.struct)
        self.current_path = os.getcwd()
        self.log = []

    def get_system_selection(self):
        data = self.sym.get_symmetry_dataset()
        predefined_systems = {
            "R3c": 1,
            "R3m": 2,
            "P4mm": 3,
            "Pmc2_1": 4,
            "Pnma": 5,
            "P6_3cm": 6
        }
        return predefined_systems.get(data["international"], 7)

    def manual_system_selection(self):
        print("Select the system:")
        print("1 - R3c")
        print("2 - R3m")
        print("3 - P4mm")
        print("4 - Pmc2_1")
        print("5 - Pnma")
        print("6 - P6_3cm")
        print("7 - Manual Orientation relationships")
        return int(input())

    def prompt_domain_wall_angle(self):
        valid_angles = self.get_valid_angles()
        print("Select the domain wall angle:")
        for angle in valid_angles:
            print(f"{angle.value} - {angle.name.replace('_', ' ')}")
        return DomainWallAngle(input())

    def prompt_polar_axis(self):
        print("Select the polar axis (a, b, or c): ")
        axis = input()
        if axis not in ['a', 'b', 'c']:
            print("Invalid input. Defaulting to 'c'.")
            axis = 'c'
        return axis

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

    def prompt_domain_size(self):
        size = float(input("Enter the domain wall size (in number of unit cells): "))
        self.log.append(f"Domain wall size: {size} unit cells")
        return size

    def prompt_supercell_size(self):
        size_a = int(input("Enter the size of the supercell along the a direction: "))
        size_b = int(input("Enter the size of the supercell along the b direction: "))
        size_c = int(input("Enter the size of the supercell along the c direction: "))
        self.log.append(f"Supercell size: a={size_a}, b={size_b}, c={size_c}")
        return [size_a, size_b, size_c]

    def prompt_cutoff_distance(self):
        cutoff = float(input("Enter the cutoff distance in angstroms to remove close atoms: "))
        self.log.append(f"Cutoff distance: {cutoff} angstroms")
        return cutoff

    def process_R3m(self, angle_selection, domain_size, cutoff, supercell_size):
        angles_dict = get_angles_dict(domain_size)["R3m"]
        space_group = "R3m"
        if angle_selection == DomainWallAngle.ALL:
            for angle_type, angles in angles_dict.items():
                process_slab(self.D1, angles, angle_type.value, space_group, cutoff, self.log, supercell_size)
        else:
            process_slab(self.D1, angles_dict[angle_selection], angle_selection.value, space_group, cutoff, self.log, supercell_size)

    def process_R3c(self, angle_selection, domain_size, cutoff, supercell_size):
        pseudo_cubic = cut(self.D1, a=[1.0, 1.0, -1], b=[-1, 1, 1], c=[1, -1, 1])
        rotate(pseudo_cubic, pseudo_cubic.cell[0], (0, 1, 0), pseudo_cubic.cell[1], (1, 0, 0))
        pseudo_cubic.write('Pseudo_cubic_R3c.vasp', sort=True, vasp5=True)
        BFO = ase.io.read('Pseudo_cubic_R3c.vasp')

        angles_dict = get_angles_dict(domain_size)["R3c"]
        space_group = "R3c"
        if angle_selection == DomainWallAngle.ALL:
            for angle_type, angles in angles_dict.items():
                process_slab(BFO, angles, f'{angle_type.value}', space_group, cutoff, self.log, supercell_size)
        else:
            process_slab(BFO, angles_dict[angle_selection], f'{angle_selection.value}', space_group, cutoff, self.log, supercell_size)

    def process_P4mm(self, angle_selection, domain_size, cutoff, supercell_size):
        angles_dict = get_angles_dict(domain_size)["P4mm"]
        space_group = "P4mm"
        if angle_selection == DomainWallAngle.ALL:
            for angle_type, angles in angles_dict.items():
                process_slab(self.D1, angles, angle_type.value, space_group, cutoff, self.log, supercell_size)
        else:
            process_slab(self.D1, angles_dict[angle_selection], angle_selection.value, space_group, cutoff, self.log, supercell_size)

    def process_Pnma(self, domain_size, cutoff, supercell_size):
        angles_dict = get_angles_dict(domain_size)["Pnma"]
        space_group = "Pnma"
        process_slab(self.D1, angles_dict[DomainWallAngle.FDW], DomainWallAngle.FDW.value, space_group, cutoff, self.log, supercell_size)

    def process_Pmc2_1(self, angle_selection, domain_size, cutoff, supercell_size, polar_axis):
        if polar_axis == 'a':
            slab1 = cut(self.D1, a=[1, 0, 0], b=[0.0, 1.0, 0], c=[0, 0, 1])
            rotate(slab1, slab1.cell[0], (0, 1, 0), slab1.cell[1], (1, 0, 0))
        elif polar_axis == 'b':
            slab1 = cut(self.D1, a=[0, 1, 0], b=[0.0, 0.0, 1], c=[1, 0, 0])
            rotate(slab1, slab1.cell[0], (0, 1, 0), slab1.cell[1], (1, 0, 0))
        elif polar_axis == 'c':
            slab1 = cut(self.D1, a=[0, 0, 1], b=[1.0, 0.0, 0], c=[0, 1, 0])
            rotate(slab1, slab1.cell[0], (0, 1, 0), slab1.cell[1], (1, 0, 0))
        else:
            print("Default c axis was selected")
            slab1 = cut(self.D1, a=[0, 0, 1], b=[1.0, 0.0, 0], c=[0, 1, 0])
            rotate(slab1, slab1.cell[0], (0, 1, 0), slab1.cell[1], (1, 0, 0))
        
        slab1.write('rot.vasp', sort=True, vasp5=True)

        Rotation = ase.io.read("rot.vasp")
        angles_dict = get_angles_dict(domain_size)["Pmc2_1"]
        space_group = "Pmc2_1"

        if angle_selection == DomainWallAngle.ALL:
            for angle_type, angles in angles_dict.items():
                process_slab(Rotation, angles, angle_type.value, space_group, cutoff, self.log, supercell_size)
        else:
            process_slab(Rotation, angles_dict[angle_selection], angle_selection.value, space_group, cutoff, self.log, supercell_size)

    def process_P6_3cm(self):
        while True:
            print("Select the domain wall type:")
            print("1 - NDW")
            print("2 - CDW")
            selection = input()
            if selection in ['1', '2']:
                script = 'ndw.py' if selection == '1' else 'cdw.py'
                command = input(f"Enter the arguments for {script} (Usage: <supercell_size> <P1_filename> <P2_filename> <output_file>): ")
                self.log.append(f"Executing {script} with arguments: {command}")
                subprocess.run(['python', script] + command.split())
                break
            else:
                print("Invalid selection. Please enter 1 or 2.")

    def create_domain_wall(self):
        filename2 = input("Enter the input domain2 file name (with extension): ")
        D2 = ase.io.read(filename2)
        self.log.append(f"Read domain2 file: {filename2}")

        current_path = os.getcwd()
        a_input = input(Fore.GREEN + Style.BRIGHT + "Enter three comma-separated values for domain 1 lattice direction a: " + Style.RESET_ALL)
        a = [float(value.strip()) for value in a_input.split(',')]
        self.log.append(f"Domain 1 lattice direction a: {a}")

        b_input = input(Fore.GREEN + Style.BRIGHT + "Enter three comma-separated values for domain 1 lattice direction b: " + Style.RESET_ALL)
        b = [float(value.strip()) for value in b_input.split(',')]
        self.log.append(f"Domain 1 lattice direction b: {b}")

        c_input = input(Fore.GREEN + Style.BRIGHT + "Enter three comma-separated values for domain 1 lattice direction c: " + Style.RESET_ALL)
        c = [float(value.strip()) for value in c_input.split(',')]
        self.log.append(f"Domain 1 lattice direction c: {c}")

        a1_input = input(Fore.BLUE + Style.BRIGHT + "Enter three comma-separated values for domain 2 lattice direction a: " + Style.RESET_ALL)
        a1 = [float(value.strip()) for value in a1_input.split(',')]
        self.log.append(f"Domain 2 lattice direction a: {a1}")

        b1_input = input(Fore.BLUE + Style.BRIGHT + "Enter three comma-separated values for domain 2 lattice direction b: " + Style.RESET_ALL)
        b1 = [float(value.strip()) for value in b1_input.split(',')]
        self.log.append(f"Domain 2 lattice direction b: {b1}")

        c1_input = input(Fore.BLUE + Style.BRIGHT + "Enter three comma-separated values for domain 2 lattice direction c: " + Style.RESET_ALL)
        c1 = [float(value.strip()) for value in c1_input.split(',')]
        self.log.append(f"Domain 2 lattice direction c: {c1}")

        slab1 = cut(self.D1, a=a, b=b, c=c)
        rotate(slab1, slab1.cell[0], (0, 1, 0), slab1.cell[1], (1, 0, 0))
        os.makedirs(os.path.join(current_path, 'Manual'))
        slab1_filename = os.path.join(current_path, 'Manual', 'Manual_domain1.vasp')
        slab1.write(slab1_filename, sort=True, vasp5=True)
        self.log.append(f"Manual domain1 slab saved to {slab1_filename}")

        slab2 = cut(D2, a=a1, b=b1, c=c1)
        rotate(slab2, slab2.cell[0], (0, 1, 0), slab2.cell[1], (1, 0, 0))
        slab2_filename = os.path.join(current_path, 'Manual', 'Manual_domain2.vasp')
        slab2.write(slab2_filename, sort=True, vasp5=True)
        self.log.append(f"Manual domain2 slab saved to {slab2_filename}")

        interface_direction = int(input("Enter the stacking direction (0 for a, 1 for b, and 2 for c): "))
        self.log.append(f"Stacking direction: {interface_direction}")

        slab = stack(slab1, slab2, axis=interface_direction, maxstrain=None)
        interface_filename = os.path.join(current_path, 'Manual', 'Manual_interface.vasp')
        slab.write(interface_filename, sort=True, vasp5=True)
        self.log.append(f"Manual interface slab saved to {interface_filename}")

    def run(self):
        print_symmetry_info(self.sym, self.log)
        self.system_selection = self.get_system_selection()
        self.log.append(f"System selection: {self.system_selection}")
        if self.system_selection == 7:
            self.system_selection = self.manual_system_selection()
            if self.system_selection == 7:
                self.create_domain_wall()
                return
        if self.system_selection == 6:
            self.process_P6_3cm()
            return
        if self.system_selection == 4:
            self.polar_axis = self.prompt_polar_axis()
            self.log.append(f"Polar axis: {self.polar_axis}")
        self.angle_selection = self.prompt_domain_wall_angle()
        self.log.append(f"Domain wall angle: {self.angle_selection}")
        self.domain_size = self.prompt_domain_size()
        self.supercell_size = self.prompt_supercell_size()
        self.cutoff = self.prompt_cutoff_distance()
        
        if self.system_selection == 1:
            self.process_R3c(self.angle_selection, self.domain_size, self.cutoff, self.supercell_size)
        elif self.system_selection == 2:
            self.process_R3m(self.angle_selection, self.domain_size, self.cutoff, self.supercell_size)
        elif self.system_selection == 3:
            self.process_P4mm(self.angle_selection, self.domain_size, self.cutoff, self.supercell_size)
        elif self.system_selection == 4:
            self.process_Pmc2_1(self.angle_selection, self.domain_size, self.cutoff, self.supercell_size, self.polar_axis)
        elif self.system_selection == 5:
            self.process_Pnma(self.domain_size, self.cutoff, self.supercell_size)

        self.log.append("Domain wall structures and supercells created successfully!")
        with open("LOGFILE.txt", "w") as logfile:
            logfile.write("\n".join(self.log))
        print("Domain wall structures and supercells created successfully! Log file written to LOGFILE.txt")

if __name__ == "__main__":
    filename = input("Enter the input file name (with extension): ")
    manager = DomainWallManager(filename)
    manager.run()
