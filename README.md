# DWBuilder: Domain Wall Builder
Ferroelectric/ferroelastic domain_wall_builder creates domain walls and interface structure for the perovskites and ceramics. The code reads the structure in the vasp format, and then creates domain walls with user-defined sizes and angles.

# Overview
There are three main scripts [dwbuilder.py](scripts/dwbuilder.py), [dbuilder.py](scripts/dbuilder.py) and [hibuilder.py](scripts/hibuilder.py) which build domain wall, single domain structures and experimentally observed atomic interface strutures with known orientation relationships (ORs), respectively. For the domain wall and interfaces, you only need to run [dwbuilder.py](scripts/dwbuilder.py). These scripts can be executed from the command line. I also have included examples of ferroelectric domain walls.

Dependencies
------------

The package requires:

-  **numpy** - http://www.numpy.org/
-  **metplotlib** - https://matplotlib.org/
-  **ASE** - https://wiki.fysik.dtu.dk/ase/

# Installation guide
For installation simply clone or download the code in your terminal and in the main directory of the package type:
```
> pip install .
```
This will copy the modules to your active python site-packages, thereby making them importable in any python script and will put the scrpits in the python bin, thereby making them executable in the shell.

# Usage
To execute the code, one should run the dwbuilder command from the command line. Upon running, the user will be prompted to input the name of the input structure file, as well as the domain wall angle and size (expressed in number of unit cells). This script builds the domain wall structure based on the space group of the input structure. For structures with the R3m or P4mm space groups, the polarization directions are assumed to align with the [001] direction. For structures with the R3c space group, a pseudo-cubic structure is first constructed by translating the input rhombohedral axes along ```a=[101], b=[-111], and c=[0-11]```. In the translated axis , the polarization direction is assumed to align with the ```[-110]``` direction. To construct 109 and 71 domain walls, each crystal system is translated according to specific orientation relationships in order to meet the polarization angle requirement.

## Example of using [dwbuilder.py](scripts/dwbuilder.py)
_ex:_

```
Enter the input file name (with extension): BiFeO3.vasp
Space group number: 161
International symbol: R3c
Lattice type: rhombohedral
Select the domain wall angle:
1 - 109 degrees
2 - 71 degrees
3 - 180 degrees
4 - All angles
4
Enter the domain wall size (in number of unit cells): 6
Domain 1 Orientation (109 DW): a = [0, 1, 0), b = [1, 0, 0], c = (0, 0, 1)
Domain 2 Orientation (109 DW): a = [0, -1, 0), b = [-1, 0, 0], c = (0, 0, 1)
Warning: This code might create domain wall artifacts, like missing or duplicate atoms, at the domain wall. Therefore, manual adjustment may be needed.
strain along a (%): 0.0
strain along b (%): 0.0
strain along c (%): 0.0
Domain 1 Orientation (71 DW): a = [1, -1, 0), b = [0, 0, 1], c = (1, 1, 0)
Domain 2 Orientation (71 DW): a = [-1, 1, 0), b = [0, 0, -1], c = (1, 1, 0)
Warning: This code might create domain wall artifacts, like missing or duplicate atoms, at the domain wall. Therefore, manual adjustment may be needed.
strain along a (%): 0.0
strain along b (%): 0.0
strain along c (%): 0.0
Domain 1 Orientation (180 DW): a = [1, 1, 0), b = [0, 0, -1], c = (-1, 1, 0)
Domain2 Orientation (180 DW): a = [-1, -1, 0), b = [0, 0, -1], c = (-1, 1, 0)
Warning: This code might create domain wall artifacts, like missing or duplicate atoms, at the domain wall. Therefore, manual adjustment may be needed.
strain along a (%): 0.0
strain along b (%): 0.0
strain along c (%): 0.0
Enter size of supercell (in number of unit cells):  
```

Once you have developed the desired domain wall structure, you can visualize in vesta or ase to further refine domain wall artifacts.  In the above examples, the domain walls are parallel along (001)//(001) planes. 

## Example of using [hibuilder.py](scripts/hibuilder.py)
In this example, I have reproduced an orientation relationship reported in the following article [CMS](https://www.sciencedirect.com/science/article/pii/S0927025621000446). 

![hibuilder](https://user-images.githubusercontent.com/52278972/235938564-5fb21315-a334-4754-ac1c-395b5daffc8f.png)

[hibuilder.py](scripts/hibuilder.py) generates an initial configuration that needs additional refinement along a, b and c-directions. To determine the optimal interface distance along the a, b, and c axes, you can either perform static DFT or forcefield simulations. Additionally, you may need to explore various interfacial configurations to identify the thermodynamically stable structure.

It is also important to know that all the domain wall and interface structures built in this package follow the orthogonality  condition, which means lattice vectors a,b and c are mutually perpendicular (i.e., orthognal) to each other and fulfill these conditions: 

```
a x b = c
a x c = b
b x c = a
```
So far the [dwbuilder.py](scripts/dwbuilder.py) develops the domain walls for the following space groups along with the assumed polarization directions:

```
R3c (Polarization direction (1-11))
R3m (Polarization direction (001))
P4mm (Polarization direction (001))
p6_3cm (Polarization direction (-1-10))
Pnma (Polarization direction (110))
Pmc2_1 (user-defined)
Amm2 (user-defined)
```


# Questions/Contributions
If you have any questions, find an issue you can get in contact with [me](mailto:zeeshan.khalid039@gmail.com).


# License
[MIT](./LICENSE).  




