#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from ase.io import read, write
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.structure import Structure
import spglib
from ase.spacegroup import get_spacegroup
from ase.io import read

# Get the current directory
dir_path = os.getcwd()

# Read in atomic coordinates from POSCAR
with open(os.path.join(dir_path, 'POSCAR'), 'r') as f:
    lines = f.readlines()

pos_atoms = []
for line in lines[8:]:
    coords = line.split()[:3]
    if len(coords) == 3:
        pos_atoms.append(list(map(float, coords)))
    else:
        break  # stop reading if line does not contain atomic coordinates

pos_atoms = np.array(pos_atoms)

# Read in atomic coordinates from CONTCAR
with open(os.path.join(dir_path, 'CONTCAR'), 'r') as f:
    lines = f.readlines()

cont_atoms = []
for line in lines[8:]:
    coords = line.split()[:3]
    if len(coords) == 3:
        cont_atoms.append(list(map(float, coords)))
    else:
        break  # stop reading if line does not contain atomic coordinates

cont_atoms = np.array(cont_atoms)

# Read the POSCAR file
poscar = read('POSCAR')
contcar = read('CONTCAR')

# Get the space group
sp_CONTCAR = get_spacegroup(contcar,symprec=1e-5)
print(sp_CONTCAR)

# Get the lattice constants
a, b, c, alpha, beta, gamma = contcar.get_cell_lengths_and_angles()

# Print the lattice constants
print('a =', a)
print('b =', b)
print('c =', c)
print('alpha =', alpha)
print('beta =', beta)
print('gamma =', gamma)


# Get the space group
sp_POSCAR = get_spacegroup(poscar,symprec=1e-5)
print(sp_POSCAR)

# Get the lattice constants
a, b, c, alpha, beta, gamma = poscar.get_cell_lengths_and_angles()

# Print the lattice constants
print('a =', a)
print('b =', b)
print('c =', c)
print('alpha =', alpha)
print('beta =', beta)
print('gamma =', gamma)

# Extract the lattice vector and coordinate information
lattice = np.array([float(lines[i].split()[j]) for j in range(3) for i in range(2, 5)]).reshape(3,3)
species = lines[5].split()
num_atoms = [int(i) for i in lines[6].split()]
coord_start = 8
coord_end = coord_start + sum(num_atoms)

# Calculate differences in fractional coordinates for each axis, taking into account periodic boundary conditions
dx = cont_atoms[:, 0] - pos_atoms[:, 0]
dy = cont_atoms[:, 1] - pos_atoms[:, 1]
dz = cont_atoms[:, 2] - pos_atoms[:, 2]

dx = np.where(dx > 0.9, (dx - 1) * a, dx * a)
dy = np.where(dy > 0.9, (dy - 1) * b, dy * b)
dz = np.where(dz > 0.9, (dz - 1) * c, dz * c)

dx -= np.round(dx)
dy -= np.round(dy)
dz -= np.round(dz)

# Write differences to file
with open('difference.txt', 'w') as f:
    for i in range(len(dx)):
        f.write(f"{dx[i]:.6f} {dy[i]:.6f} {dz[i]:.6f}\n")

# Plot contour plots of differences along each axis
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4))
im1 = ax1.tricontourf(pos_atoms[:, 1], pos_atoms[:, 2], dx, cmap='rainbow', levels=500)
im2 = ax2.tricontourf(pos_atoms[:, 0], pos_atoms[:, 2], dy, cmap='rainbow', levels=500)
im3 = ax3.tricontourf(pos_atoms[:, 0], pos_atoms[:, 1], dz, cmap='rainbow', levels=500)

# Set plot labels and titles
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_title('dx')

ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_title('dy')

ax3.set_xlabel('x')
ax3.set_ylabel('y')
ax3.set_title('dz')

# Add a colorbar
divider = make_axes_locatable(ax3)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(im3, cax=cax)
cbar.set_label('Difference')

# Get the nominal charge values from the user for each species of atoms
charge_values = {}
for element in set(species):
    charge_values[element] = int(input(f"Enter the nominal charge for {element}: "))

# Apply the nominal charge values to each species of atoms
new_lines = lines[:coord_start]
for i in range(len(species)):
    element = species[i]
    charge = charge_values[element]
    num = num_atoms[i]
    for j in range(num):
        line = lines[coord_start + j]
        x, y, z = [float(k) for k in line.split()[:3]]
        new_lines.append(f"{x} {y} {z} {charge}\n")
    coord_start += num

# Write the modified CHG file
with open('CHARGE', 'w') as f:
    f.writelines(new_lines)

plt.tight_layout()
plt.show()


with open('CHARGE') as f:
    lines = f.readlines()

# extract column 4 starting from line 8
charge = np.array([list(map(float, line.split()))[3] for line in lines[8:]])

# Read dipole_moment.txt file and calculate polarization
diff_data = np.loadtxt('difference.txt')

Conv = 1602.176; # Value of Conv
dipole_moment = np.multiply(diff_data[:, :3], charge[:, np.newaxis])
polarization = np.sum(dipole_moment, axis=0) / (a * b * c)*Conv  # divide by cell volume to get polarization
polarization = np.abs(polarization)  # take absolute value to get positive polarization

# Save polarization array to file and print output
np.savetxt('dipole_moments.txt', dipole_moment)
np.savetxt('polarization.txt', polarization)
print("Final polarization along x-axis: {:.2f} µC/cm²".format(polarization[0]))
print("Final polarization along y-axis: {:.2f} µC/cm²".format(polarization[1]))
print("Final polarization along z-axis: {:.2f} µC/cm²".format(polarization[2]))
