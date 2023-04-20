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


# Find the path of the current directory
current_path = os.getcwd()

filename = input("Enter the input file name (with extension): ")
D1 = ase.io.read(filename)

struct = Structure.from_file(filename)
sym = SpacegroupAnalyzer(struct)
data = sym.get_symmetry_dataset()

print("Space group number: {}".format(data["number"]))
print("International symbol: {}".format(data["international"]))
print(f"Lattice type: {sym.get_lattice_type()}")

# Conditionally select the system based on international symbol
if data["international"] == "R3c":
    system_selection = 1
elif data["international"] == "R3m":
    system_selection = 2
elif data["international"] == "P4mm":
    system_selection = 3
elif data["international"] == "P6_3cm":
    system_selection = 4
elif data["international"] == "Pnma":
    system_selection = 5

else:
    print("Warning: The international symbol does not match any pre-defined systems.")
    print("Please manually select the system.")
    system_selection = 6

# Prompt the user to select the system if not automatically selected
if system_selection == 6:
    print("Select the system:")
    print("1 - R3c")
    print("2 - R3m")
    print("3 - P4mm")
    print("4 - P6_3cm")
    print("5 - Pnma")
    print("6 - Create slab")

    system_selection = int(input())

# BFO system
if system_selection == 1:
    # Prompt the user to enter the input file name
    filename = input("Enter the input file name (with extension): ")
    # Read the input file
    BFO1 = ase.io.read(filename)

    # Cut and rotate the structure to obtain a pseudo-cubic structure
    pseudo_cubic = cut(BFO1, a=[1.0, 1.0, -1], b=[-1, 1, 1], c=[1, -1, 1])
    rotate(pseudo_cubic, pseudo_cubic.cell[0], (0,1,0), pseudo_cubic.cell[1], (1,0,0))

    # Write the pseudo-cubic structure to a VASP file
    pseudo_cubic.write('Pseudo_cubic_BFO.vasp', sort=True, vasp5=True)

    # Read the pseudo-cubic structure
    BFO = ase.io.read('Pseudo_cubic_BFO.vasp')

    # Prompt the user to select the domain wall angle
    print("Select the domain wall angle:")
    print("1 - 109 degrees")
    print("2 - 71 degrees")
    print("3 - 180 degrees")
    print("4 - All angles")
    angle_selection = int(input())

    # Prompt the user to enter the domain wall size
    domain_size = int(input("Enter the domain wall size (in number of unit cells): "))

    # Cut and rotate the slabs based on the selected domain wall angle(s)
    if angle_selection == 1 or angle_selection == 4:
        slab1 = cut(BFO, a=[0.0, 1.01, 0], b=[1.01, 0, 0], c=[0, 0, domain_size])
        rotate(slab1, slab1.cell[0], (0,1,0), slab1.cell[1], (1,0,0))
        slab2 = cut(BFO, a=[0.0, -1.01, 0], b=[-1.01, 0, 0], c=[0, 0, domain_size])
        rotate(slab2, slab2.cell[0], (0,1,0), slab2.cell[1], (1,0,0))
        os.makedirs('D_109')
        slab1.write('D_109/109_domain1.vasp', sort=True, vasp5=True)
        slab2.write('D_109/109_domain2.vasp', sort=True, vasp5=True)
        print("Domain 1 Orientation (109 DW): a = [0, 1, 0), b = [1, 0, 0], c = (0, 0, 1)")
        print("Domain 2 Orientation (109 DW): a = [0, -1, 0), b = [-1, 0, 0], c = (0, 0, 1)")

    if angle_selection == 2 or angle_selection == 4:
        slab1 = cut(BFO, a=[1.0, -1.0, 0], b=[0, 0, 1.0], c=[domain_size, domain_size, 0])
        rotate(slab1, slab1.cell[0], (0,1,0), slab1.cell[1], (1,0,0))
        slab2 = cut(BFO, a=[-1.0, 1.0, 0], b=[0, 0, -1.0], c=[domain_size, domain_size, 0])
        rotate(slab2, slab2.cell[0], (0,1,0), slab2.cell[1], (1,0,0))
        os.makedirs(os.path.join(current_path, 'D_71'))
        slab1.write(os.path.join(current_path, 'D_71', '71_domain1.vasp'), sort=True, vasp5=True)
        slab2.write(os.path.join(current_path, 'D_71', '71_domain2.vasp'), sort=True, vasp5=True)
        print("Domain 1 Orientation (71 DW): a = [1, -1, 0), b = [0, 0, 1], c = (1, 1, 0)")
        print("Domain 2 Orientation (71 DW): a = [-1, 1, 0), b = [0, 0, -1], c = (1, 1, 0)")


    if angle_selection == 3 or angle_selection == 4:
        slab1 = cut(BFO, a=[1.0, 1.0, 0], b=[0, 0, -1.0], c=[-domain_size, domain_size, 0])
        rotate(slab1, slab1.cell[0], (0,1,0), slab1.cell[1], (1,0,0))
        slab2 = cut(BFO, a=[-1.0, -1.0, 0], b=[0, 0, 1.0], c=[-domain_size, domain_size, 0])
        rotate(slab2, slab2.cell[0], (0,1,0), slab2.cell[1], (1,0,0))
        os.makedirs(os.path.join(current_path, 'D_180'))    
        slab1.write(os.path.join(current_path,'D_180', '180_domain1.vasp'), sort=True, vasp5=True)
        slab2.write(os.path.join(current_path,'D_180', '180_domain2.vasp'), sort=True, vasp5=True)
        print("Domain 1 Orientation (180 DW): a = [1, 1, 0), b = [0, 0, -1], c = (-1, 1, 0)")
        print("Domain2 Orientation (180 DW): a = [-1, -1, 0), b = [0, 0, -1], c = (-1, 1, 0)")

        

    #Print confirmation message to user

print("Domain  structures created successfully!")

 # BaTiO3 system
if system_selection == 2:
    # Prompt the user to enter the input file name
    filename = input("Enter the input file name (with extension): ")
    # Read the input file
    BTO1 = ase.io.read(filename)

    # Prompt the user to select the domain wall angle
    print("Select the domain wall angle:")
    print("1 - R180 degrees")
    print("2 - R71 degrees")
    print("3 - R109 degrees")
    print("4 - All angles")
    angle_selection = int(input())

    # Prompt the user to enter the lattice vector sizes
    print("Note: R109 and  R71 domain walls lie along the a-direction and R180 along c-direction")

    # Prompt the user to enter the lattice vector sizes
    domain_size = int(input("Enter the domain wall size (in number of unit cells): "))

    # Cut and rotate the slabs based on the selected domain wall angle(s)
    if angle_selection == 1 or angle_selection == 4:
        slab1 = cut(BTO1, a=[1, 1, 0], b=[0, 0, 1], c=[domain_size, -domain_size, 0])
        rotate(slab1, slab1.cell[0], (0,1,0), slab1.cell[1], (1,0,0))
        slab2 = cut(BTO1, a=[-1, -1, 0], b=[0, 0, -1], c=[domain_size, -domain_size, 0])
        rotate(slab2, slab2.cell[0], (0,1,0), slab2.cell[1], (1,0,0))
        os.makedirs(os.path.join(current_path, 'R180'))
        slab1.write(os.path.join(current_path, 'R180', 'R180_domain1.vasp'), sort=True, vasp5=True)
        slab2.write(os.path.join(current_path, 'R180', 'R180_domain2.vasp'), sort=True, vasp5=True)
        print("Domain 1 Orientation (R180 DW): a = [1, 1, 0), b = [0, 0, 1], c = (1, -1, 0)")
        print("Domain 2 Orientation (R180 DW): a = [-1, -1, 0), b = [0, 0, -1], c = (1, -1, 0)")

    
    if angle_selection == 2 or angle_selection == 4:
        slab1 = cut(BTO1, a=[domain_size, domain_size, 0], b=[0, 0, 1], c=[1, -1, 0])
        rotate(slab1, slab1.cell[0], (0,1,0), slab1.cell[1], (1,0,0))
        slab2 = cut(BTO1, a=[domain_size, domain_size, 0], b=[0, 0, -1], c=[-1, 1, 0])
        rotate(slab2, slab2.cell[0], (0,1,0), slab2.cell[1], (1,0,0))
        os.makedirs(os.path.join(current_path, 'R71'))
        slab1.write(os.path.join(current_path, 'R71', 'R71_domain1.vasp'), sort=True, vasp5=True)
        slab2.write(os.path.join(current_path, 'R71', 'R71_domain2.vasp'), sort=True, vasp5=True)
        print("Domain 1 Orientation (R71 DW): a = (1, 1, 0), b = [0, 0, 1], c = [1, -1, 0]")
        print("Domain 2 Orientation (R71 DW): a = (1, 1, 0), b = [0, 0, -1], c = [-1, 1, 0]")



    if angle_selection == 3 or angle_selection == 4:
        slab1 = cut(BTO1, a=[domain_size, 0.0, 0], b=[0, 1, 0], c=[0, 0, 1])
        rotate(slab1, slab1.cell[0], (0,1,0), slab1.cell[1], (1,0,0))
        slab2 = cut(BTO1, a=[domain_size, 0.0, 0], b=[0, -1, 0], c=[0, 0, -1])
        rotate(slab2, slab2.cell[0], (0,1,0), slab2.cell[1], (1,0,0))
        os.makedirs(os.path.join(current_path, 'R109'))
        slab1.write(os.path.join(current_path, 'R109', 'R109_domain1.vasp'), sort=True, vasp5=True)
        slab2.write(os.path.join(current_path, 'R109', 'R109_domain2.vasp'), sort=True, vasp5=True)
        print("Domain 1 Orientation (R109 DW): a = (1, 0, 0), b = [0, 1, 0], c = [0, 0, 1]")
        print("Domain 2 Orientation (R109 DW): a = (1, 0, 0), b = [0, -1, 0], c = [0, 0, -1]")



    #Print confirmation message to user

    print("Domain structures created successfully!")



 # PTO system
if system_selection == 3:
    # Prompt the user to enter the input file name
    filename = input("Enter the input file name (with extension): ")
    # Read the input file
    PT = ase.io.read(filename)
    
    # Prompt the user to select the domain wall angle
    print("Select the domain wall angle:")
    print("1 - T180 degrees")
    print("2 - T90 degrees")
    print("3 - All angles")
    angle_selection = int(input())
    # Prompt the user to enter the lattice vector sizes
    domain_size = int(input("Enter the domain wall size (in number of unit cells): "))

    # Cut and rotate the slabs based on the selected domain wall angle(s)
    if angle_selection == 1 or angle_selection == 3:
        slab1 = cut(PT, a=[1.01, 0, 0], b=[0, domain_size, 0], c=[0, 0, 1.01])
        rotate(slab1, slab1.cell[0], (0,1,0), slab1.cell[1], (1,0,0))
        slab2 = cut(PT, a=[-1.01, 0, 0], b=[0, domain_size, 0], c=[0, 0, -1.01])
        rotate(slab2, slab2.cell[0], (0,1,0), slab2.cell[1], (1,0,0))
        os.makedirs(os.path.join(current_path, 'T180'))
        slab1.write(os.path.join(current_path, 'T180', 'T180_domain1.vasp'), sort=True, vasp5=True)
        slab2.write(os.path.join(current_path, 'T180', 'T180_domain2.vasp'), sort=True, vasp5=True)
        print("Domain 1 Orientation (T180 DW): a = [1, 0, 0], b = (0, 1, 0), c = [0, 0, 1]")
        print("Domain 2 Orientation (T180 DW): a = [-1, 0, 0], b = (0, 1, 0), c = [0, 0, -1]")

    
    if angle_selection == 2 or angle_selection == 3:
        slab1 = cut(PT, a=[0, 1, 0], b=[-1, 0, 1], c=[domain_size, 0, domain_size])
        rotate(slab1, slab1.cell[0], (0,1,0), slab1.cell[1], (1,0,0))
        slab2 = cut(PT, a=[0, -1, 0], b=[1, 0, -1], c=[domain_size, 0, domain_size])
        rotate(slab2, slab2.cell[0], (0,1,0), slab2.cell[1], (1,0,0))
        slab90 = stack(slab1, slab2, axis=2, maxstrain=None)
        os.makedirs(os.path.join(current_path, 'T90'))
        slab1.write(os.path.join(current_path, 'T90', 'T90_domain1.vasp'), sort=True, vasp5=True)
        slab2.write(os.path.join(current_path, 'T90', 'T90_domain2.vasp'), sort=True, vasp5=True)
        print("Domain 1 Orientation (T90 DW): a = [0, 1, 0], b = [-1, 0, 1], c = (1, 0, 1)")
        print("Domain 2 Orientation (T90 DW): a = [0, -1, 0], b = [1, 0, -1], c = (1, 0, 1)")



    #Print confirmation message to user

    print("Domain wall structures created successfully!")



    
    
# YMO system
if system_selection == 4:
    # Prompt the user to enter the input file name
    filename1 = input("Enter the input domain1 file name (with extension): ")
    filename2 = input("Enter the input domain2 file name (with extension): ")

    # Read the input file
    YMO1 = ase.io.read(filename1)
    YMO2 = ase.io.read(filename2)


    # Prompt the user to select the domain wall angle
    print("Select the domain wall angle:")
    print("1 - Neutral domain wall")
    print("2 - Charged domain wall")
    print("3 - both")
    angle_selection = int(input())

    # Find the path of the current directory
    current_path = os.getcwd()



    domain_size = float(input("Enter the domain wall size (in number of unit cells): "))

    # Cut and rotate the slabs based on the selected domain wall angle(s)
    if angle_selection == 1 or angle_selection == 3:
        slab1 = cut(YMO1, a=[-domain_size, domain_size, 0], b=[0.0, 0.0, 1.01], c=[-1.0, -1.0, 0])
        rotate(slab1, slab1.cell[0], (0,1,0), slab1.cell[1], (1,0,0))
        slab2 = cut(YMO2, a=[-domain_size, domain_size, 0], b=[0.0, 0.0, 1.01], c=[-1.0, -1.0, 0])
        rotate(slab2, slab2.cell[0], (0,1,0), slab2.cell[1], (1,0,0))
        os.makedirs(os.path.join(current_path, 'NDW'))
        slab1.write(os.path.join(current_path, 'NDW', 'domain1.vasp'), sort=True, vasp5=True)
        slab2.write(os.path.join(current_path, 'NDW', 'domain2.vasp'), sort=True, vasp5=True)
        print("Domain Orientation: a = (-1, 1, 0), b = [0, 0, 1], c = [-1, -1, 0]")

    if angle_selection == 2 or angle_selection == 3:
        slab1 = cut(YMO1, a=[-1, 1, 0], b=[0.0, 0.0, domain_size+0.01], c=[-1.0, -1.0, 0])
        rotate(slab1, slab1.cell[0], (0,1,0), slab1.cell[1], (1,0,0))
        slab2 = cut(YMO2, a=[-1, 1, 0], b=[0.0, 0.0, domain_size+0.01], c=[-1.0, -1.0, 0])
        rotate(slab2, slab2.cell[0], (0,1,0), slab2.cell[1], (1,0,0))
        os.makedirs(os.path.join(current_path, 'CDW'))
        slab1.write(os.path.join(current_path, 'CDW', 'CDW_domain1.vasp'), sort=True, vasp5=True)
        slab2.write(os.path.join(current_path, 'CDW', 'CDW_domain2.vasp'), sort=True, vasp5=True)
        print("Domain Orientation: a = [-1, 1, 0], b = (0, 0, 1), c = [-1, -1, 0]")

        print("\033[1;31;40mWarning: This code could generate domain wall artifacts (i.e., missing Oxygen atoms, duplicate atoms) at the domain wall, thus it requires manual adjustment.\033[0m")

    #Print confirmation message to user

    print("Domain structures created successfully!")

 # CTO system
if system_selection == 5:
    # Prompt the user to enter the input file name


    # Prompt the user to enter the lattice vector sizes
    print("Pnma space group has ferroelastic domain walls")

    # Prompt the user to enter the lattice vector sizes
    domain_size = int(input("Enter the domain wall size (in number of unit cells): "))

    # Cut and rotate the slabs based on the selected domain wall angle(s)
    slab1 = cut(D1, a=[1, 1, 0], b=[0.0, 0.0, 1], c=[domain_size, -domain_size, 0])
    rotate(slab1, slab1.cell[0], (0,1,0), slab1.cell[1], (1,0,0))
    slab2 = cut(D1, a=[-1, 1, 0], b=[0.0, 0.0, 1], c=[domain_size, domain_size, 0])
    rotate(slab2, slab2.cell[0], (0,1,0), slab2.cell[1], (1,0,0))
    os.makedirs(os.path.join(current_path, 'Pnma'))
    slab1.write(os.path.join(current_path, 'Pnma', 'D1.vasp'), sort=True, vasp5=True)
    slab2.write(os.path.join(current_path, 'Pnma', 'D2.vasp'), sort=True, vasp5=True)
    print("Domain 1 Orientation (R180 DW): a = [1, 1, 0), b = [0, 0, 1], c = (1, -1, 0)")
    print("Domain 2 Orientation (R180 DW): a = [-1, 1, 0), b = [0, 0, 1], c = (1, 1, 0)")
    print("\033[1;31;40mWarning: This code might create domain artifacts, like missing or duplicate atoms. Therefore, manual adjustment may be needed.\033[0m")


# Manually develop interfaces/DW if you have ORs
if system_selection == 6:
    def create_domain_wall(filename1, filename2, domain_size):

    # Find the path of the current directory
     current_path = os.getcwd()


    # Prompt the user to select the domain wall angle
    print("Select the relevant option:")
    print("1 - Unitcell")
    print("2 - Manual OR single slab")
    print("3 - Manual OR two slabs")

    angle_selection = int(input())



    
    # initialize colorama
    init()
    
            #Creating planes
    if angle_selection == 1:
        
        filename = input("Enter the input file name (with extension): ")
        D1 = ase.io.read(filename1)

        # define lattice directions for slab
        size_a = int(input("Enter the supercell size along a-direction (in number of unit cells): "))
        size_b = int(input("Enter the supercell size along b-direction (in number of unit cells): "))
        size_c = int(input("Enter the supercell size along c-direction (in number of unit cells): "))

        # create and manipulate slab1 model
        slab1 = cut(D1, a=[size_a, 0, 0], b=[0.0, size_b, 0], c=[0, 0, size_c])
        rotate(slab1, slab1.cell[0], (0,1,0), slab1.cell[1], (1,0,0))
        vacuum_thickness = float(input("Enter the vacuum thickness in Angstroms:"))  # in Angstroms
        vacuum_direction = int(input("Enter the vacuum direction 0 for a, 1 for b and 2 for c): "))
        slab1.center(vacuum=vacuum_thickness, axis=vacuum_direction)
        os.makedirs(os.path.join(current_path, 'Slab'))
        slab1.write(os.path.join(current_path, 'Slab', 'supercell.vasp'), sort=True, vasp5=True)
        print("Slab created successfully!")


        #Creating planes
    if angle_selection == 2:
        
        filename = input("Enter the input file name (with extension): ")
        D1 = ase.io.read(filename1)


        # define lattice directions for slab1
        a_input = input(Fore.GREEN + Style.BRIGHT + "Enter three comma-separated values for domain 1 lattice direction a: " + Style.RESET_ALL)
        a = [float(value.strip()) for value in a_input.split(',')]

        b_input = input(Fore.GREEN + Style.BRIGHT + "Enter three comma-separated values for domain 1 lattice direction b: " + Style.RESET_ALL)
        b = [float(value.strip()) for value in b_input.split(',')]

        c_input = input(Fore.GREEN + Style.BRIGHT + "Enter three comma-separated values for domain 1 lattice direction c: " + Style.RESET_ALL)
        c = [float(value.strip()) for value in c_input.split(',')]


        # create and manipulate slab1 model
        slab1 = cut(D1, a=a, b=b, c=c)
        rotate(slab1, slab1.cell[0], (0,1,0), slab1.cell[1], (1,0,0))
        vacuum_thickness = float(input("Enter the vacuum thickness in Angstroms:"))  # in Angstroms
        vacuum_direction = int(input("Enter the vacuum direction 0 for a, 1 for b and 2 for c): "))
        slab1.center(vacuum=vacuum_thickness, axis=vacuum_direction)
        os.makedirs(os.path.join(current_path, 'OR1'))
        slab1.write(os.path.join(current_path, 'OR1', 'OR1.vasp'), sort=True, vasp5=True)
        print("Slab created successfully!")

        #view(slab1)



    #Creating planes
    if angle_selection == 3:
        
        filename1 = input("Enter the input bulk 1 file name (with extension): ")
        filename2 = input("Enter the input bulk 2 file name (with extension): ")

        D1 = ase.io.read(filename1)
        D2 = ase.io.read(filename2)



        # define lattice directions for slab1
        a_input = input(Fore.GREEN + Style.BRIGHT + "Enter three comma-separated values for domain 1 lattice direction a: " + Style.RESET_ALL)
        a = [float(value.strip()) for value in a_input.split(',')]

        b_input = input(Fore.GREEN + Style.BRIGHT + "Enter three comma-separated values for domain 1 lattice direction b: " + Style.RESET_ALL)
        b = [float(value.strip()) for value in b_input.split(',')]

        c_input = input(Fore.GREEN + Style.BRIGHT + "Enter three comma-separated values for domain 1 lattice direction c: " + Style.RESET_ALL)
        c = [float(value.strip()) for value in c_input.split(',')]

        # define lattice directions for slab2
        a_input = input(Fore.BLUE + Style.BRIGHT + "Enter three comma-separated values for domain 2 lattice direction a: " + Style.RESET_ALL)
        a1 = [float(value.strip()) for value in a_input.split(',')]

        b_input = input(Fore.BLUE + Style.BRIGHT + "Enter three comma-separated values for domain 2 lattice direction b: " + Style.RESET_ALL)
        b1 = [float(value.strip()) for value in b_input.split(',')]

        c_input = input(Fore.BLUE + Style.BRIGHT + "Enter three comma-separated values for domain 2 lattice direction c: " + Style.RESET_ALL)
        c1 = [float(value.strip()) for value in c_input.split(',')]

        # create and manipulate slab1 model
        slab1 = cut(D1, a=a, b=b, c=c)
        rotate(slab1, slab1.cell[0], (0,1,0), slab1.cell[1], (1,0,0))
        vacuum_thickness = float(input("Enter the vacuum thickness for slab1 in Angstroms:"))  # in Angstroms
        vacuum_direction = int(input("Enter the vacuum direction 0 for a, 1 for b and 2 for c): "))
        slab1.center(vacuum=vacuum_thickness, axis=vacuum_direction)
        os.makedirs(os.path.join(current_path, 'OR1-2'))
        slab1.write(os.path.join(current_path, 'OR1-2', 'slab1.vasp'), sort=True, vasp5=True)
        print("Slab1 created successfully!")

        #view(slab1)

        # create and manipulate slab2 model
        slab2 = cut(D2, a=a1, b=b1, c=c1)
        rotate(slab2, slab2.cell[0], (0,1,0), slab2.cell[1], (1,0,0))
        vacuum_thickness = float(input("Enter the vacuum thickness for slab2 in Angstroms:"))  # in Angstroms
        vacuum_direction = int(input("Enter the vacuum direction 0 for a, 1 for b and 2 for c): "))
        slab2.center(vacuum=vacuum_thickness, axis=vacuum_direction)
        slab2.write(os.path.join(current_path, 'OR1-2', 'slab2.vasp'), sort=True, vasp5=True)
        print("Slab2 created successfully!")

        #view(slab2)

    

