---
title: 'DWBuilder: A code to generate ferroelectric/ferroelastic domain walls and multi-material atomic interface structures'
tags:
  - Python
  - domain walls
  - crystallography
  - multi-material inrerfaces
  - atomistic simulations (VASP, LAMMPS)
authors:
  - name: M.Z.Khalid
    orcid: https://orcid.org/0000-0002-7866-3870
    affiliation: 1
  - name: S.M.Selbach
    affiliation: 1
affiliations:
 - name: Department of Materials Science and Engineering, Norwegian University of Science and Technology, Trondheim, Norway
   index: 1
date: 22 April 2024
bibliography: paper.bib
---

# Summary

Ferroelectric materials have an order parameter called polarization that can be toggled using an external electric field. Regions within the material that exhibit consistent polarization are known as domain, and the boundaries between different domains are referred to as domain walls (DWs).  These domain walls, which are only a few nanometers wide, possess unique properties with potential technological applications. DWs show promise for nanoscale electronic circuit elements and enable innovative design concepts because they can be created, erased, and manipulated using electric fields [@meier2022ferroelectric; @catalan2012domain; @meier2015functional; @bednyakov2018physics] . DWs can also replicate the functions of key electronic components such as diodes [@whyte2015diode], transistors [@mundy2017functional], and random access memory (RAM) [@sharma2017nonvolatile] .

Due to the small sizes and exciting properties of DWs, there has been significant interest in studying how to control and manipulate them using atomistic simulations [@schultheiss2020intrinsic; @smaabraaten2020domain; @smaabraaten2018charged]. However, developing atomic DW structures is challenging and requires knowledge and understanding of the order parameter and DW types in ferroelectric materials. DWs can be ferroelectric, antiferroelectric, and/or ferroelastic, and they can vary depending on the allowed symmetry of the ferroelectric material. For instance, ferroelectric BiFeO$_3$ exists at room temperature as a rhombohedrally distorted perovskite with space group R3c and an internal angle of 89.23°. The spontaneous polarization is oriented along the [111]$_P$ axis[@ederer2005effect], [@wang2003epitaxial]. The ferroelectric phase transition in BiFeO$_3$ causes lattice distortion along the <111>$_P$ polarization direction, resulting in four types of DWs with polarization direction changes of 71°, 109°, or 180°[@wang2003epitaxial]. Similarly, other domain wall types have been identified in ferroelectric materials such as BaTiO$_3$[@taherinejad2012bloch], PbTiO$_3$[@meyer2002ab], YMnO$_3$[@smaabraaten2018charged], and ferroelastic DWs in CaTiO$_3$[@barone2014improper].

``DWBuilder`` code is designed as a command-line tool to create DWs and interface structures from specific input unit cell geometries, as described in detail in the README file of the repository. The code comprises two main components: (i) a domain wall builder for similar materials and (ii) a heterogeneous interface builder for multi-material atomic interfaces. Figure 1 explains the structure and workflow of the ``DWBuilder`` package. 

The first part, handled by the scripts ``dwbuilder.py`` and ``dbuilder.py``,  produces domain walls by first analyzing the input unit cell geometry and determining the space group of the structure. The space group is identified using the open-source Python library Pymatgen. If the space group matches the specified type, the script offers a range of possible domain wall types and ultimately creates the DW structures. If the space group of the input structure does not match, the script allows you to choose a desired space group type or manually define the domains by specifying lattice vectors. To generate different domains separately, you can use ``dbuilder.py`` to develop distinct domains, which can be useful for bulk and surface calculations of individual domains.
  
 The second part of the code involves creating a heterogeneous interface structure of multi-material compounds, which is handled by the script ``hibuilder.py``. This script requires two input structures, named ``bulk1`` and ``bulk2``. To develop compatible interfaces, you must define the orientation relationship (OR) between the two bulk phases. This definition is necessary to address any lattice and angular mismatches that arise from differences in space groups and/or atomic structures of the two phases. 

Currently, the script cannot predict the ORs that would result in a low lattice mismatch between the two interfaces. However, theoretical studies and methods such as edge-to-edge [@zhang2005edge] and face-to-face [@khalid2020atomistic] matching techniques can help predict low lattice misfit for interface construction. This script assumes that the user is already familiar with the appropriate ORs to construct the interface structure. For instance, the ORs of interfaces reported in referenced papers  [@khalid2020ab; @khalid2021imc; @khalid2021first] can be replicated using this script. Additionally, the script can generate atomic interfaces if you know the OR from experiments, and it can predict the atomic interface structure and lattice mismatch between the two bulk phases.

 ![Structure of the `DWBuilder` package.\label{fig:scheme}](dwbuilder.png){width="100%"}

# Statement of need:
``DWBuilder`` is an interactive toolbox for developing atomic domain walls and interface structures of homogeneous and heterogeneous material compounds, making it suitable for high-throughput calculations. Its target audience includes students and scientists in materials science and physics at any level of expertise. ``DWBuilder`` utilizes the NumPy library extensively, which speeds up execution, particularly when working with large structures. Users are guided through the process of identifying and creating the desired domain walls in a step-by-step manner. The code is designed to be user-friendly and educational, with a focus on plane orientation and electric polarization switching.


# Acknowledgements

Funding from the Norwegian research counsil  (Grant agreement No. 90544501) is gratefully acknowledged.

# References



