![GitHub release (latest by date)](https://img.shields.io/github/v/release/mzkhalid039/DWBuilder)
![Python Version](https://img.shields.io/badge/python-3.11-blue)
[![status](https://joss.theoj.org/papers/eeae372c75342085351f92af4ca165db/status.svg)](https://joss.theoj.org/papers/eeae372c75342085351f92af4ca165db)


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
-  **Pymatgen** - https://pymatgen.org/

# Installation guide
For installation simply clone or download the code in your terminal and in the main directory of the package type:
```
> pip install .
```
This will copy the modules to your active python site-packages, thereby making them importable in any python script and will put the scrpits in the python bin, thereby making them executable in the shell.

After installation, you can start the software by running the following command in your terminal:
```
> dwbuilder
```
# Method for the development of FDW
## R3m
Ferroelectric materials with a rhombohedral space group of R3m exhibit three types of ferroelectric domain walls (FDWs): R71, R109, and R180. Due to neutrality and mechanical compatibility, these FDWs develop along specific planes. For the R180 FDW, the ```{1-10}``` plane is used, while for R71 and R109, the FDWs are determined by the sum of the two polarization vectors.

- **R71 FDW**: Lies parallel to the diagonal ```(1-10)``` plane in the primitive cell.
- **R180 FDW**: Developed along the ```{1-10}``` plane.
- **R109 FDW**: The R109 FDW lies along the ```(100)``` plane, and FDWs are simply developed by stacking the primitive cell along this direction.

To analyze the R71 and R180 FDWs, the pseudocubic rhombohedral unit cell is transformed into a 10-atom cell by rotating it 45° around the z-axis relative to the parent cubic unit cell. The transformation matrix applied is given as:

        [ 1  0  1 ]                
    T = [ 1  0 -1 ]       
        [ 0  1  0 ]        

In the script the above transformation matrix is applied in lattice vectors form as given below:

```
a=[1, 1, 0], b=[0, 0, 1], c=[1, -1, 0]

```

**R71**

The R71 domain walls are developed by stacking along the normal direction to the ```[1-10]``` plane, which aligns with the a-axis in the rotated coordinate system. Additionally, the mirror symmtry was developed by tranforming ```[1-10]``` to ```[-110]``` and from ```[001]``` to ```[00-1]```.

For the R71 FDWs, the orientation relationships for both domains (D1 and D2) are:

         [ 1  0  1 ]         [ 1  0 -1 ]          
    D1 = [ 1  0 -1 ]    D2 = [ 1  0  1 ]    
         [ 0  1  0 ]         [ 0 -1  0 ]  

**R180**

Similar to the R71 domain wall, the R180 ```{1-10}``` wall is also parallel to the diagonal plane in the primitive cell. For the R180 FDW:
- A supercell is made by stacking the rotated 10-atom unit cell.
- The final FDWs are developed by stacking along the ```[1-10]``` plane.

```
         [ 1  0  1 ]         [ -1  0   1 ]          
    D1 = [ 1  0 -1 ]    D2 = [ -1  0  -1 ]    
         [ 0  1  0 ]         [  0 -1   0 ]  

```

The transformation ensures that the polarization vectors in the two domains are oriented 180 degrees apart, achieving the necessary mirror symmetry for the R180 domain wall. 

**R109**

 The R109 FDW lies along the ```(100)``` plane. For this FDW:
 
- The domain walls are developed by stacking the primitive cell along the ```(100)``` direction.

This straightforward stacking ensures that the polarization vectors align as required for the R109 domain wall configuration.

```
         [ 1  0  0 ]         [  1  0  0 ]          
    D1 = [ 0  1  0 ]    D2 = [  0 -1  0 ]    
         [ 0  1  1 ]         [  0  0 -1 ]  

```

## Reference: Taherinejad M, Vanderbilt D, Marton P, Stepkova V, Hlinka J. Bloch-type domain walls in rhombohedral BaTiO3. Physical Review B. 2012 Oct 22;86(15):155138.

## R3c
Similar to the R3m, the R3c structure also possesses three FDWs: R71, R109, and R180.  The R3c domain walls are developed by converting the reference R3c structure into a pseudocubic cell. This is done by cutting along the ```(1-11)``` plane, which is also the polarization axis.

The transformation vectors used to convert the R3c structure into a pseudocubic cell are given by the transformation matrix Y:

```
        [ 1 -1  1 ]                
    T = [ 1  1 -1 ]       
        [-1  1  1 ]        

```

The R71 , R109 , and R180 domain walls are found to adopt ```(001)```, ```(110)```, and ```(-110)``` pseudocubic domain wall planes, respectively. The transformation matrix to develop respective domains and domain walls are given as:

```

             [ 1  0  1 ]              [-1  0  1 ]          
    R71_D1 = [-1  0  1 ]    R71_D2 =  [ 1  0  1 ]    
             [ 0  1  0 ]              [ 0 -1  0 ]  

              [ 0  1  0 ]              [ 0 -1  0 ]          
    R109_D1 = [ 1  0  0 ]    R109_D2 = [-1  0  0 ]    
              [ 0  0  1 ]              [ 0  0  1 ]  
                      
              [ 1  0 -1 ]              [-1  0 -1 ]          
    R180_D1 = [ 1  0  1 ]    R180_D2 = [-1  0  1 ]    
              [ 0 -1  0 ]              [ 0  1  0 ]  

```
The R109 domain walls are formed by stacking the D1 and D2 domains along the ```[001]```, while R71 and R180 are stacked along the ```[110]``` and ```[-110]``` directions, respectively. 

## Reference: Wang Y, Nelson C, Melville A, Winchester B, Shang S, Liu ZK, Schlom DG, Pan X, Chen LQ. BiFeO 3 domain wall energies and structures: a combined experimental and density functional theory+ U study. Physical review letters. 2013 Jun 24;110(26):267601.

## P4mm
In the P4mm space group, with the polarization direction along the `[001]` direction, there are two types of domain walls allowed and developed in `dwbuilder`:

1. **T180 Domain Walls (T180)**:
   - **Orientation**: These domain walls are on the `(100)` planes.
   - **Polarization**: The polarization vectors in neighboring domains are parallel but oriented in opposite directions.

2. **T90 Domain Walls**:
   - ***Orientation**: These domain walls are on the `(101)` planes.
   - ***Polarization**: The polarization vectors on either side of the domain wall are almost perpendicular to each other.

The transformation matrix used for both domain walls are given as:

```

              [ 1  0  0 ]               [ -1  0   0 ]          
    T180_D1 = [ 0  1  0 ]    T180_D2 =  [  0  1   0 ]    
              [ 0  0  1 ]               [  0  0  -1 ]  

             [ 0  -1  1 ]             [ 0  -1   1 ]          
    T90_D1 = [ 1   0  0 ]    T90_D2 = [ -1  0   0 ]    
             [ 0   1  1 ]             [ 0  -1   1 ]  

```

The T180 and T90 FDWs are develed by stacking along the ```[010]``` and ```[101]``` direction.  

## Reference: Meyer B, Vanderbilt D. Ab initio study of ferroelectric domain walls in PbTiO 3. Physical Review B. 2002 Mar 1;65(10):104111.

## P6_3cm

The hexagonal manganites RMnO3 ( R= Sc, Y, Dy-Lu) belongs to the P6_3cm hexagonal space group. This space group is non-centrosymmetric, thus lacking inversion symmetry, which is a necessary condition for ferroelectricity. The polar axis for this group is along the `[001]` direction.

### Origin of Polarization

The polarization within hexagonal manganites originates from the structural tilt trimerization of the MnO5 polyhedra. Trimerization is characterized by two angles:
1. **Tilt Angle (θ)**:
   - This is the angle between the `[001]` direction and the Mn-apical O bond.
   - It gives the magnitude of the tilting and can be regarded as the order parameter.
2. **Azimuthal Angle (φ)**:
   - This angle describes the orientation of the tilting.
   - During the phase transition, (φ) adopts one of six values separated by 60°, corresponding to trimerization around three possible origins with either in- or out-orientation of the tilting.

Based on the stacking site, there exist two kinds of domain walls within hexagonal manganites:

1. **Neutral Domain Walls (NDW)**:
   - **Formation**: NDWs are developed by stacking the positive and negative domains perpendicular to the `[001]` axis.
   - **Polarization**: This stacking creates net zero polarization at the domain walls due to the cancellation of the positive and negative polarization of both domains.

2. **Charged Domain Walls (CDW)**:
   - **Formation**: CDWs are formed by stacking two different domains along the `[001]` direction.
   - **Polarization**: This results in a net polarization at the domain walls, leading to charged domain walls.

### Developing Domain Walls with `dwbuilder`

In `dwbuilder`, we have developed two scripts, `ndw.py` and `cdw.py`, to create these domain walls while maintaining hexagonal symmetry. To run these scripts, users need to provide the input domain structures for both domains.

- **`ndw.py`**: This script is used to develop neutral domain walls (NDW).
- **`cdw.py`**: This script is used to develop charged domain walls (CDW).

By providing the necessary input structures, these scripts will generate the respective domain walls within the hexagonal manganites. 

**Note**: You do not explicitly need to run these scripts directly. `dwbuilder` automatically detects the space group and, based on the space group, prompts users to provide details of the input structure and parameters.

## Reference: Kumagai Y, Spaldin NA. Structural domain walls in polar hexagonal manganites. Nature communications. 2013 Feb 26;4(1):1540.

# Usage
To execute the code, one should run the dwbuilder command from the command line. Upon running, the user will be prompted to input the name of the input structure file, as well as the domain wall angle and size (expressed in number of unit cells). This script builds the domain wall structure based on the space group of the input structure. For structures with the R3m or P4mm space groups, the polarization directions are assumed to align with the ```[001]``` direction. For structures with the R3c space group, a pseudo-cubic structure is first constructed by translating the input rhombohedral axes along ```a=[101], b=[-111], and c=[0-11]```. In the translated axis , the polarization direction is assumed to align with the ```[-110]``` direction.  To construct 109 and 71 domain walls, each crystal system is translated according to specific orientation relationships in order to meet the polarization angle requirement.

However, if you are unsure about the polarization dirction, you can run ```polarization.py``` to determine the polarization direction. This script uses the point charge model and you need to have both optimized ```CONTCAR``` and high symmtry structure ```POSCAR``` files to determine the polarization.

## Example of using [dwbuilder.py](scripts/dwbuilder.py)
_ex:_

```
> dwbuilder

    ***************************************************
    *                                                 *
    *                Welcome to DWBuilder             *
    *                                                 *
    ***************************************************
    *        Hey, you must know what you are doing.   *
    *  Otherwise, you might get wrong results.        *
    ***************************************************
    *     Core Developer: M.Z.Khalid                  *
    *     Main Contributors: S.M.Selbach              *
    *     Email: zeeshan.khalid039@gmail.com          *
    ***************************************************
    
Select a script to run:
1: dwbuilder.py                5: polarization.py
2: dbuilder.py                 6: supercell.py
3: hibuilder.py                7: vasp2cif.py
4: slab.py                     8: bondanalysis.py
Enter the number of the script to run: 1
Enter the input file name (with extension): BiFeO3_R3c.vasp
Space group number: 161
International symbol: R3c
Lattice type: rhombohedral
Select the domain wall angle:
R180 - R180
R71 - R71
R109 - R109
ALL - ALL
ALL
Enter the domain wall size (in number of unit cells): 3
Enter the size of the supercell along the a direction: 1
Enter the size of the supercell along the b direction: 1
Enter the size of the supercell along the c direction: 1
Enter the cutoff distance in angstroms to remove close atoms: 0
Lattice strain for R109 DW:
strain along a (%): 0.00
strain along b (%): 0.00
strain along c (%): 0.00
Lattice strain for R71 DW:
strain along a (%): 0.00
strain along b (%): 0.00
strain along c (%): 0.00
Lattice strain for R180 DW:
strain along a (%): 0.00
strain along b (%): 0.00
strain along c (%): 0.00
Domain wall structures and supercells created successfully! Log file written to LOGFILE.txt
```

Once you have developed the desired domain wall structure, you can visualize in vesta or ase to further refine domain wall artifacts.  In the above examples, the domain walls are parallel along ```(001)//(001)``` planes.

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


## Example of using [hibuilder.py](scripts/hibuilder.py)
In this example, I have reproduced an orientation relationship reported in the following article [CMS](https://www.sciencedirect.com/science/article/pii/S0927025621000446). 

```
> dwbuilder

    ***************************************************
    *                                                 *
    *                Welcome to DWBuilder             *
    *                                                 *
    ***************************************************
    *        Hey, you must know what you are doing.   *
    *  Otherwise, you might get wrong results.        *
    ***************************************************
    *     Core Developer: M.Z.Khalid                  *
    *     Main Contributors: S.M.Selbach              *
    *     Email: zeeshan.khalid039@gmail.com          *
    ***************************************************
    
Select a script to run:
1: dwbuilder.py                5: polarization.py
2: dbuilder.py                 6: supercell.py
3: hibuilder.py                7: vasp2cif.py
4: slab.py                     8: bondanalysis.py
Enter the number of the script to run: 3
Enter the bulk phase 1file name (with extension): Fe_unitcell.vasp
Enter the bulk phase 2 file name (with extension): Fe2Al5.vasp
Bulk phase 1 Space group number: 229
Bulk phase 1 International symbol: lm-3m
Bulk phase 1 Lattice type: cubic
Bulk phase 2 Space group number: 63
Bulk phase 2 International symbol: Cmcm
Bulk phase 2 Lattice type: orthrhombic
Enter three comma-separated values for bulk 1 lattice direction a: 1,0,1
Enter three comma-separated values for bulk 1 lattice direction b: 1.5,1.5,-1.5
Enter three comma-separated values for bulk 1 lattice direction c: -2,4,2
Enter three comma-separated values for bulk 2 lattice direction a: 0,0,-1
Enter three comma-separated values for bulk 2 lattice direction b: 1,0,0
Enter three comma-separated values for bulk 2 lattice direction c: 0,-2,0
Enter the stacking direction 0 for a, 1 for b and 2 for c): 2
strain along a (%): 1.1378
strain along b (%): -0.5492
strain along c (%): 8.2378
Angular strain along a (radians): 0.0
Angular strain along b (radians): 0.0
Angular strain along c (radians): 0.0
Warning: This code could generate interface/domain wall artifacts (i.e., Oxygen atoms, duplicate atoms etc,) at the interface, thus it requires manual adjustment.
```
[hibuilder.py](scripts/hibuilder.py) generates an initial configuration that needs additional refinement along a, b and c-directions. To determine the optimal interface distance along the a, b, and c axes, you can either perform static DFT or forcefield simulations. Additionally, you may need to explore various interfacial configurations to identify the thermodynamically stable structure.

The above code generates the following interface structure (the periodic image) visualized by [VESTA](https://jp-minerals.org/vesta/en/):

![image](https://github.com/mzkhalid039/DWBuilder/assets/52278972/19081982-1949-4e16-a734-ca27a6d733bb)

It is also important to know that all the domain wall and interface structures built in this package follow the orthogonality condition, which means lattice vectors a,b and c are mutually perpendicular (i.e., orthognal) to each other and fulfill these conditions: 

```
a x b = c
a x c = b
b x c = a
```

- _**A note on the microscopic degrees of freedom:**_

To find the global minimum of an interface structure, it is essential to calculate the optimal distance between the different bulk phases at the interface. You also need to identify the termination site that corresponds to the structure with the lowest energy, as well as consider the homogeneous interface structures of the periodic images.

To achieve these goals, I often use first-principles methods such as density functional theory (DFT), specifically with the Vienna Ab-initio Simulation Package [VASP](https://www.vasp.at/) at 0 Kelvin. This computational approach helps determine the initial structures for further first-principles calculations. 

## Relevant Papers using ```DWBuilder```

1. Khalid MZ, Williamson BAD, Eggested K, Selbach SM. **Hydrogen Doping at Proper/Improper Ferroelectric Domain Walls.** To be submitted. 2024.

2. [Khalid MZ, Friis J, Ninive PH, Marthinsen K, Ringdalen IG, Strandlie A. **First-principles study of tensile and shear strength of an Fe2Al5//Fe interface.** Computational Materials Science. 2021 May 1;192:110319](https://www.sciencedirect.com/science/article/pii/S0927025621000446).

3. [Khalid MZ, Friis J, Ninive PH, Marthinsen K, Strandlie A. **First-principles study of tensile and shear strength of Fe-Al and α-AlFeSi intermetallic compound interfaces.** Computational Materials Science. 2021 Feb 1;187:110058](https://www.sciencedirect.com/science/article/pii/S0927025620305498).

4. [Khalid MZ, Friis J, Ninive PH, Marthinsen K, Strandlie A. **Ab-initio study of atomic structure and mechanical behaviour of Al/Fe intermetallic interfaces.** Computational Materials Science. 2020 Mar 1;174:109481](https://www.sciencedirect.com/science/article/pii/S0927025619307803).

5. [Khalid MZ, Friis J, Ninive PH, Marthinsen K, Strandlie A. **A first-principles study of the Al (001)/Fe (0-11) interface.** In Materials Science Forum 2019 Jan 27 (Vol. 941, pp. 2349-2355). Trans Tech Publications Ltd.](https://www.scientific.net/MSF.941.2349)


# Questions/Contributions
If you have any questions, find an issue you can get in contact with [me](mailto:zeeshan.khalid039@gmail.com).


# License
[MIT](./LICENSE).  




