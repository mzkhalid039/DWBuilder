#!/usr/bin/env python
"""
Script for bonding analysis in a VASP format structure.

This script reads one or two VASP format structure files, determines the bonding between
different elements based on a cutoff distance, and plots the bonding network
and a histogram of bond lengths.

Dependencies:
- ase
- matplotlib
- networkx
- numpy

Usage:
- Ensure the required dependencies are installed and the script is executed in a Python environment.
- Run the script and provide the input file name when prompted.
"""

import ase.io
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from ase.neighborlist import neighbor_list
from matplotlib.lines import Line2D
from collections import defaultdict

def read_vasp_structure(file_name):
    """Read the VASP structure file."""
    return ase.io.read(file_name)

def determine_bonding(structure, cutoff):
    """Determine bonds based on a cutoff distance, considering periodic boundary conditions."""
    indices_i, indices_j, offsets = neighbor_list('ijS', structure, cutoff)
    bond_dict = defaultdict(list)
    bond_lengths_dict = defaultdict(list)
    for i, j, S in zip(indices_i, indices_j, offsets):
        distance = structure.get_distance(i, j, mic=True)
        if distance <= cutoff:
            elem1 = structure[i].symbol
            elem2 = structure[j].symbol
            bond_type = tuple(sorted([elem1, elem2]))
            bond_dict[bond_type].append((i, j, elem1, elem2))
            bond_lengths_dict[bond_type].append(distance)
    return bond_dict, bond_lengths_dict

def plot_bonding_network(structure, bond_dict, shortest_bond=None, longest_bond=None, title_suffix=""):
    """Plot the bonding network."""
    G = nx.Graph()
    positions = {}

    for atom in structure:
        G.add_node(atom.index, element=atom.symbol)
        positions[atom.index] = atom.position[:2]

    for bond_type, bonds in bond_dict.items():
        color = np.random.rand(3,)  # Random color for each bond type
        for bond in bonds:
            i, j, elem1, elem2 = bond
            if shortest_bond and ((i, j) == shortest_bond or (j, i) == shortest_bond):
                G.add_edge(i, j, color='green', width=3.0)
            elif longest_bond and ((i, j) == longest_bond or (j, i) == longest_bond):
                G.add_edge(i, j, color='red', width=3.0)
            else:
                G.add_edge(i, j, color=color)

    fig, ax = plt.subplots(figsize=(10, 10))
    elements = list(set(atom.symbol for atom in structure))

    colors = plt.cm.jet(np.linspace(0, 1, len(elements)))
    color_map = {element: color for element, color in zip(elements, colors)}

    node_colors = [color_map[G.nodes[node]['element']] for node in G.nodes]
    edge_colors = [G[u][v]['color'] for u, v in G.edges]

    nx.draw(G, pos=positions, node_color=node_colors, with_labels=True, ax=ax, node_size=500, font_size=10, edge_color=edge_colors)

    legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, label=element)
                       for element, color in color_map.items()]
    
    bond_legend_elements = [Line2D([0], [0], color=color, lw=2, label='-'.join(bond_type)) for bond_type, color in zip(bond_dict.keys(), edge_colors)]
    if shortest_bond:
        bond_legend_elements.append(Line2D([0], [0], color='green', lw=2, label='Shortest Bond'))
    if longest_bond:
        bond_legend_elements.append(Line2D([0], [0], color='red', lw=2, label='Longest Bond'))

    ax.legend(handles=legend_elements + bond_legend_elements, loc='upper right')

    plt.title(f"Bonding Network {title_suffix} (Node numbers represent atom indices)")
    plt.show()

def plot_bond_length_histogram(bond_lengths_dict, title_suffix=""):
    """Plot a histogram of bond lengths."""
    plt.figure(figsize=(10, 6))
    for bond_type, bond_lengths in bond_lengths_dict.items():
        plt.hist(bond_lengths, bins=30, alpha=0.5, label='-'.join(bond_type))
    plt.xlabel('Bond Length (Å)')
    plt.ylabel('Frequency')
    plt.title(f'Histogram of Bond Lengths {title_suffix}')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_combined_bond_length_histogram(bond_lengths_dict1, bond_lengths_dict2, title_suffix1="", title_suffix2=""):
    """Plot combined histogram of bond lengths."""
    plt.figure(figsize=(10, 6))
    for bond_type, bond_lengths in bond_lengths_dict1.items():
        plt.hist(bond_lengths, bins=30, alpha=0.5, label=f'{title_suffix1} - {"-".join(bond_type)}', color='blue')
    for bond_type, bond_lengths in bond_lengths_dict2.items():
        plt.hist(bond_lengths, bins=30, alpha=0.5, label=f'{title_suffix2} - {"-".join(bond_type)}', color='orange')
    plt.xlabel('Bond Length (Å)')
    plt.ylabel('Frequency')
    plt.title('Combined Histogram of Bond Lengths')
    plt.legend()
    plt.grid(True)
    plt.show()

def analyze_bond_statistics(bond_lengths_dict):
    """Analyze and display bond length statistics."""
    all_bond_lengths = [length for lengths in bond_lengths_dict.values() for length in lengths]
    mean_length = np.mean(all_bond_lengths)
    median_length = np.median(all_bond_lengths)
    std_length = np.std(all_bond_lengths)
    shortest_bond_length = np.min(all_bond_lengths)
    longest_bond_length = np.max(all_bond_lengths)

    print(f"Mean bond length: {mean_length:.2f} Å")
    print(f"Median bond length: {median_length:.2f} Å")
    print(f"Standard deviation of bond lengths: {std_length:.2f} Å")
    print(f"Shortest bond length: {shortest_bond_length:.2f} Å")
    print(f"Longest bond length: {longest_bond_length:.2f} Å")

    # Identify shortest and longest bonds
    shortest_bond = None
    longest_bond = None
    for bond_type, lengths in bond_lengths_dict.items():
        for length in lengths:
            if length == shortest_bond_length:
                shortest_bond = bond_type
            if length == longest_bond_length:
                longest_bond = bond_type

    return shortest_bond, longest_bond

def compare_bond_lengths(bond_lengths_dict1, bond_lengths_dict2):
    """Compare bond lengths between two structures and highlight differences."""
    common_bond_types = set(bond_lengths_dict1.keys()).intersection(bond_lengths_dict2.keys())
    bond_length_changes = {}

    print("Percentage change in bond lengths:")
    for bond_type in common_bond_types:
        lengths1 = np.array(bond_lengths_dict1[bond_type])
        lengths2 = np.array(bond_lengths_dict2[bond_type])
        percentage_change = ((lengths2 - lengths1) / lengths1) * 100
        bond_length_changes[bond_type] = percentage_change
        print(f"{'-'.join(bond_type)}: Mean change = {np.mean(percentage_change):.2f}%, Median change = {np.median(percentage_change):.2f}%")

    unique_bond_types1 = set(bond_lengths_dict1.keys()).difference(bond_lengths_dict2.keys())
    unique_bond_types2 = set(bond_lengths_dict2.keys()).difference(bond_lengths_dict1.keys())

    if unique_bond_types1:
        print("\nUnique bonds in structure 1:")
        for bond_type in unique_bond_types1:
            print(f"{'-'.join(bond_type)}")

    if unique_bond_types2:
        print("\nUnique bonds in structure 2:")
        for bond_type in unique_bond_types2:
            print(f"{'-'.join(bond_type)}")

    return bond_length_changes

def compare_structures(structure1, structure2, cutoff):
    """Compare the bonding of two structures."""
    bond_dict1, bond_lengths_dict1 = determine_bonding(structure1, cutoff)
    bond_dict2, bond_lengths_dict2 = determine_bonding(structure2, cutoff)

    total_bonds1 = sum(len(bonds) for bonds in bond_dict1.values())
    total_bonds2 = sum(len(bonds) for bonds in bond_dict2.values())
    print(f"Found {total_bonds1} bonds in structure 1.")
    print(f"Found {total_bonds2} bonds in structure 2.")

    shortest_bond1, longest_bond1 = analyze_bond_statistics(bond_lengths_dict1)
    shortest_bond2, longest_bond2 = analyze_bond_statistics(bond_lengths_dict2)

    plot_bonding_network(structure1, bond_dict1, shortest_bond1, longest_bond1, title_suffix="(Structure 1)")
    plot_bonding_network(structure2, bond_dict2, shortest_bond2, longest_bond2, title_suffix="(Structure 2)")
    plot_combined_bond_length_histogram(bond_lengths_dict1, bond_lengths_dict2, title_suffix1="Structure 1", title_suffix2="Structure 2")

    bond_length_changes = compare_bond_lengths(bond_lengths_dict1, bond_lengths_dict2)

def main():
    """Main function to run the bonding analysis."""
    comparison = input("Do you want to compare two structures? (yes/no): ").strip().lower()
    
    if comparison == "yes":
        file_name1 = input("Enter the first VASP structure file name (with extension): ")
        file_name2 = input("Enter the second VASP structure file name (with extension): ")
        structure1 = read_vasp_structure(file_name1)
        structure2 = read_vasp_structure(file_name2)
        cutoff = float(input("Enter the cutoff distance for bonding analysis (in angstroms): "))
        compare_structures(structure1, structure2, cutoff)
    else:
        file_name = input("Enter the VASP structure file name (with extension): ")
        structure = read_vasp_structure(file_name)
        cutoff = float(input("Enter the cutoff distance for bonding analysis (in angstroms): "))
        bond_dict, bond_lengths_dict = determine_bonding(structure, cutoff)
        total_bonds = sum(len(bonds) for bonds in bond_dict.values())
        print(f"Found {total_bonds} bonds.")
        
        shortest_bond, longest_bond = analyze_bond_statistics(bond_lengths_dict)
        plot_bonding_network(structure, bond_dict, shortest_bond, longest_bond)
        plot_bond_length_histogram(bond_lengths_dict)

if __name__ == "__main__":
    main()
