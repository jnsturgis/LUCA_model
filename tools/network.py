"""
This file describes a metabolic network class and ancillary functions.

The objective of this work is to build whole organism (or meta-organism)
metabolic models and study how these can evolve and adapt to different
conditions.

In this file and the associated class we will try to adhere to the following
scheme.

Instance methods of the class modify the data.
Class methods operate on the class or are alternatice constructors
Standalone functions operate with multiple instances or are very generic.
"""

import sys

import networkx as nx

import reaction as rxn
import compound as cpd

class Network:
    """
    A class to describe metabolic networks
    """
    def __init__(self, *argv):
        match len(argv):
            case 2:
                self.compounds = argv[0]
                self.reactions = argv[1]
            case 0:
                self.compounds = {}
                self.reactions = {}
            case _:
                raise ValueError("Need 0 or 2 dictionaries to make a metabolic network.")

    @classmethod
    def from_csv(cls, csv_name):
        """
        Use a csv file to setup a metabolic network.
        """
        compounds = {}
        reactions = {}

        with (sys.stdin if csv_name == "-" else open(csv_name, "r", encoding='utf-8')) as f:
            for line in f:
                line = line.rstrip()
                parts = line.split(',',1)
                match line[0]:
                    case '#':           # Comment
                        continue
                    case 'R':           # Reaction line
                        reactions[parts[0]] = rxn.Reaction.from_text(line)
                    case 'M':           # Metabolite line
                        compounds[parts[0]] = cpd.Compound.from_text(line)
                    case _:             # Oops error
                        raise ValueError(f'Invalid input line \n{line}')
        return cls(compounds, reactions)

    def add_compound( self, compound_id, compound_data ):
        """
        Add a new compound to a metabolic network.
        """
        self.compounds[compound_id] = compound_data

    def add_reaction(self, reaction_id, reaction_data ):
        """
        Add a new reaction to a metabolic network.
        """
        self.reactions[reaction_id] = reaction_data

def find_modifiers(network):
    """
    Return the set of all modifiers that feature in the reactions of the network.
    """
    modifiers = []
    for _, reaction in network.reactions.items():
        modifiers.append(reaction.modifiers)
    return set(modifiers)

def make_graph( network):
    """
    Make a network graph from the metabolic network information.
    """
    graph = nx.Graph()
    for r_id, reaction in network.reactions.items():
        for c_id in reaction.substrates():
            graph.add_edge(c_id, r_id)
        for c_id in reaction.products():
            graph.add_edge(r_id, c_id)
    return graph

def add_atoms( dict1, dict2 ):
    """
    Sum by atom type the values in dict1 and dict2
    """
    keys = set(dict1) | set(dict2)  # Union of keys
    result = {k: dict1.get(k, 0) + dict2.get(k, 0) for k in keys}
    return result

def check_rxn_balance(reaction, compounds):
    """
    Check the balance of a reaction for atoms and charges
    """
    charge = 0.0
    atoms = {}
    for substrate in reaction.subst:
        c_id = substrate[0]
        num  = substrate[1]
        compound = compounds[c_id]
        charge += num * float(compound.charge)
        new_atoms = compound.atoms()
        for k, v in new_atoms.items():
            new_atoms[k] = v * num
        atoms = add_atoms(atoms, new_atoms)
    for product in reaction.prod:
        c_id = product[0]
        num  = product[1] * -1
        compound = compounds[c_id]
        charge += num * float(compound.charge)
        new_atoms = compound.atoms()
        for k, v in new_atoms.items():
            new_atoms[k] = v * num
        atoms = add_atoms(atoms, new_atoms)
    if not (charge == 0.0 and all(value == 0 for value in atoms.values())):
        s1 = f'Reaction {reaction.r_id} unbalanced: '
        if charge != 0.0:
            s1 += 'charge {charge}'
        for k, v in atoms:
            if v != 0.0:
                s1 += f",{k} {v}"
        print( s1 )
        return False
    return True

def check_connectivity(network, exclude):
    """
    Check that the metabolic network is connected and indicate divisions if there are any

    network is a Network structure
    exclude is a list of ids to remove from the graph before checking for sub_graphs.

    Returns the networkx undirected graph of the network based on reactions without
    any unused compounds that might feature as modifiers.
    """
    graph = make_graph(network)
    nodes = list(graph.nodes)
    unrefed = find_modifiers(network).union(nodes)
    for c_id, _ in network.compounds.items():
        if c_id not in unrefed:
            print(f'Compound {c_id} is not used by reactions.')
    for item in exclude:
        if graph.has_node(item):
            graph.remove_node(item)
    sub_graphs = list(nx.connected_components(graph))
    print(f'Network has {len(sub_graphs)} connected components.')
    if len(sub_graphs) > 1:
        for item in sub_graphs:
            print(item)
    return graph

def check_balance(network):
    """
    Check that the reactions are balanced for charge and ChemicalFormula and warn of problems
    """
    balanced = True
    for _, reaction in network.reactions.items():
        balanced = balanced and check_rxn_balance(reaction, network.compounds)
    if balanced:
        print( 'All reactions are chemically balanced.')

def analyse_pathway(network, exclude, *args):
    """
    Perform various analyses on the network! A bit generic.

    1. Find loops
    2. Notify if no loop metabolite is a reaction modulator
    """
    if len(args) == 1:
        graph = args[0]
    else:
        graph = make_graph( network)
    for item in exclude:
        if graph.has_node(item):
            graph.remove_node(item)
    cycles = list(nx.cycle_basis(graph))
    print(f'Network has {len(cycles)} base cycles.')
    modifiers = find_modifiers(network)
    for ring in cycles:
        mods = modifiers.intersection(ring)
        if mods:
            print (f"Modifiers {mods} appear in cycle: {ring}")
        else:
            print (f"No modifiers in cycle: {ring}")
    return graph

def main():
    """
    The main routine for the network file
    """
    print("This file is used to define a network class")

if __name__ == '__main__':
    main()
