import os
import sys
import ase
import ase.io
from ase.build import cut, rotate, stack
import numpy as np
from colorama import init, Fore, Style
from numpy.linalg import norm
from colorama import Fore, Style
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


filename1 = input("Enter the bulk phase 1file name (with extension): ")
filename2 = input("Enter the bulk phase 2 file name (with extension): ")


# Read the input file
D1 = ase.io.read(filename1)
D2 = ase.io.read(filename2)


struct = Structure.from_file(filename1)
sym = SpacegroupAnalyzer(struct)
data = sym.get_symmetry_dataset()

print("Bulk phase 1 Space group number: {}".format(data["number"]))
print("Bulk phase 1 International symbol: {}".format(data["international"]))
print(f"Bulk phase 1 Lattice type: {sym.get_lattice_type()}")

struct = Structure.from_file(filename2)
sym = SpacegroupAnalyzer(struct)
data = sym.get_symmetry_dataset()

print("Bulk phase 2 Space group number: {}".format(data["number"]))
print("Bulk phase 2 International symbol: {}".format(data["international"]))
print(f"Bulk phase 2 Lattice type: {sym.get_lattice_type()}")


# Manually develop interfaces/DW if you have ORs

# Find the path of the current directory
current_path = os.getcwd()
    
# initialize colorama
init()

#Creating planes

# define lattice directions for slab1
a_input = input(Fore.GREEN + Style.BRIGHT + "Enter three comma-separated values for bulk 1 lattice direction a: " + Style.RESET_ALL)
a = [float(value.strip()) for value in a_input.split(',')]

b_input = input(Fore.GREEN + Style.BRIGHT + "Enter three comma-separated values for bulk 1 lattice direction b: " + Style.RESET_ALL)
b = [float(value.strip()) for value in b_input.split(',')]

c_input = input(Fore.GREEN + Style.BRIGHT + "Enter three comma-separated values for bulk 1 lattice direction c: " + Style.RESET_ALL)
c = [float(value.strip()) for value in c_input.split(',')]

 # define lattice directions for slab2
a_input = input(Fore.BLUE + Style.BRIGHT + "Enter three comma-separated values for bulk 2 lattice direction a: " + Style.RESET_ALL)
a1 = [float(value.strip()) for value in a_input.split(',')]

b_input = input(Fore.BLUE + Style.BRIGHT + "Enter three comma-separated values for bulk 2 lattice direction b: " + Style.RESET_ALL)
b1 = [float(value.strip()) for value in b_input.split(',')]

c_input = input(Fore.BLUE + Style.BRIGHT + "Enter three comma-separated values for bulk 2 lattice direction c: " + Style.RESET_ALL)
c1 = [float(value.strip()) for value in c_input.split(',')]

# create and manipulate slab1 model
slab1 = cut(D1, a=a, b=b, c=c)
rotate(slab1, slab1.cell[0], (0,1,0), slab1.cell[1], (1,0,0))
os.makedirs(os.path.join(current_path, 'HIS'))
slab1.write(os.path.join(current_path, 'HIS', 'bulk1.vasp'), sort=True, vasp5=True)

    #view(slab1)

# create and manipulate slab2 model
slab2 = cut(D2, a=a1, b=b1, c=c1)
rotate(slab2, slab2.cell[0], (0,1,0), slab2.cell[1], (1,0,0))
slab2.write(os.path.join(current_path, 'HIS', 'bulk2.vasp'), sort=True, vasp5=True)

#view(slab2)

interface_direction = int(input("Enter the stacking direction 0 for a, 1 for b and 2 for c): "))

slab = stack(slab1, slab2, axis=interface_direction, maxstrain=None)
slab.write(os.path.join(current_path, 'HIS', 'interface.vasp'), sort=True, vasp5=True)
    

slab.write('H-interface.vasp', sort=True, vasp5=True)


# Calculate lattice strain

a1, b1, c1 = slab1.cell
a2, b2, c2 = slab2.cell

print ('strain along a (%):', (norm(a1) - norm(a2)) / norm(a2) * 100.)
print ('strain along b (%):', (norm(b1) - norm(b2)) / norm(b2) * 100.)
print ('strain along c (%):', (norm(c1) - norm(c2)) / norm(c2) * 100.)

# Calculate angular strain
theta_a = np.arccos(np.dot(a1, a2) / (np.linalg.norm(a1) * np.linalg.norm(a2))) 
theta_b = np.arccos(np.dot(b1, b2) / (np.linalg.norm(b1) * np.linalg.norm(b2))) 
theta_c = np.arccos(np.dot(c1, c2) / (np.linalg.norm(c1) * np.linalg.norm(c2)))
print('Angular strain along a (radians):', theta_a) 
print('Angular strain along b (radians):', theta_b) 
print('Angular strain along c (radians):', theta_c)

print("\033[1;31;40mWarning: This code could generate interface/domain wall artifacts (i.e., Oxygen atoms, duplicate atoms etc,) at the interface, thus it requires manual adjustment.\033[0m")

