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

# pylint: disable=W0511
# TODO Needs unit testing
import sys

import numpy as np

import networkx as nx
from networkx.algorithms import bipartite

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
        old_line = ""

        sep = ";"

        with (sys.stdin if csv_name == "-" else open(csv_name, "r", encoding='utf-8')) as f:
            for line in f:
                line = line.rstrip()
                if len(old_line) > 0:
                    line = old_line + line
                if line[-1] != '\\':        # Continuation character
                    parts = line.split(sep,1)
                    match line[0]:
                        case '#':           # Comment
                            continue
                        case 'R':           # Reaction line
                            if parts[0] in reactions:
                                print(f'Duplicate reaction {parts[0]}.')
                            reactions[parts[0]] = rxn.Reaction.from_text(line, sep)
                        case 'M':           # Metabolite line
                            if parts[0] in compounds:
                                print(f'Duplicate compound {parts[0]}.')
                            compounds[parts[0]] = cpd.Compound.from_text(line, sep)
                        case _:             # Oops error
                            raise ValueError(f'Invalid input line \n{line}')
                    old_line = ""
                else:
                    old_line = line[:-1]
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
        modifiers.extend(reaction.modifiers())
    return set(modifiers)

def make_graph( network ):
    """
    Make a network graph from the metabolic network information using networkx.

    In the resulting graph, produced from the reactions, the nodes correspond
    to the reactions and compounds featuring as substrates and products of the
    reactions. The graph is undirected and unweighted. The node labels
    correspond to the reaction and compound id's is used to set the bipartite
    attrbute (metabolites = 0, other things = 1) allowing consideration as a
    bipartite graph.

    The initial letter codes for id's are:
        * Mx_ Metabolite in compartment x (x=i or e)
        * R_  Reaction normal type.
        * X_  Exchange reaction for input from the environment.
        * T_  Transport reaction between compartments.
        * B_  Biomass reaction.

    """
    graph = nx.Graph()
    for r_id, reaction in network.reactions.items():
        for c_id in reaction.substrates():
            graph.add_edge(c_id, r_id)
        for c_id in reaction.products():
            graph.add_edge(r_id, c_id)
    # set attribute bipartite for nodes based on first letter of node_id...
    nx.set_node_attributes(graph,
        {node: 0 if node.startswith("M") else 1 for node in graph.nodes}, "bipartite")
    assert nx.is_bipartite(graph)
    return graph

def make_dgraph( network, **kwargs):
    """
    Make a directed network graph from the metabolic network information using networkx.

    In the resulting graph, produced from the reactions, the nodes correspond
    to the reactions and compounds featuring as substrates and products of the
    reactions.
    The graph is directed and weighted based on the standard free energy change
    associated with the reaction in the forward and reverse directions using
    the metropolis scheme.
    TODO different weighting options (flux, dG°'m etc)
    The node labels correspond to the reaction and compound id's is used to set
    the bipartite attrbute (reations = 1, other things = 0) allowing
    consideration as a bipartite graph.
    """
    graph = nx.DiGraph()

    options = {
        'weighting': 'dG0'
    }
    if kwargs:
        options.update({k: v for k, v in kwargs.items() if k in options})

    for r_id, reaction in network.reactions.items():

        match options['weighting']:
    # TODO add weighting algorithm and check others
            case 'dG0':
                forward_weight = 1.0        # This should be based of DG°'m'
                reverse_weight = 0.0
            case 'flux':
                flux = reaction.dict.get('flux', 0.0)
                if flux > 0:
                    forward_weight = flux   # This should be based of DG°'m'
                    reverse_weight = 0.0
                else:
                    forward_weight = 0.0
                    reverse_weight = -flux
            case '_':
                raise ValueError("Unrecognized weighting scheme.")

        for c_id in reaction.substrates():
            graph.add_weighted_edges_from([c_id, r_id, forward_weight])
            graph.add_weighted_edges_from([r_id, c_id, reverse_weight])
        for c_id in reaction.products():
            graph.add_weighted_edges_from([r_id, c_id, forward_weight])
            graph.add_weighted_edges_from([c_id, r_id, reverse_weight])
    # set attribute bipartite for nodes based on first letter of node_id...
    nx.set_node_attributes(graph,
        {node: 0 if node.startswith("M") else 1 for node in graph.nodes}, "bipartite")
    assert nx.is_bipartite(graph)
    return graph

def reaction_adjacency_graph( network ):
    """
    Calculate the reaction adjacency matrix for a metabolic network.

    This is the projection of the metabolic network bipartite graph onto the
    reaction node, with edges weighted by the number of connections.
    """
    # TODO check weighting algorithm
    bi_graph = make_graph( network )
    result = bipartite.weighted_projected_graph(
        bi_graph, nodes = {n for n, d in bi_graph.nodes(data=True) if d["bipartite"] == 1})
    return result

def normalized_flow_graph( network ):
    """
    Calculate the normalized_flow_graph for a metabolic network based on weighted graph

    This is the projection of the weighted directions metabolic network
    bipartite graph onto the reaction node, with edges weighted by the edge
    weightings.
    """
    # TODO check algorithm
    bi_graph = make_dgraph( network )
    result = bipartite.weighted_projected_graph(
        bi_graph, nodes = {n for n, d in bi_graph.nodes(data=True) if d["bipartite"] == 1})
    return result

def mass_flow_graph( network ):
    """
    Calculate the reaction adjacency matrix for a metabolic network.

    This is the projection of the metabolic network bipartite graph onto the
    reaction node, with edges weighted by the number of connections.
    """
    # TODO check algorithm
    bi_graph = make_dgraph( network, weighting='flux')
    result = bipartite.weighted_projected_graph(
        bi_graph, nodes = {n for n, d in bi_graph.nodes(data=True) if d["bipartite"] == 1})
    return result

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
            s1 += f'charge {charge}'
        for k, v in atoms.items():
            if v != 0.0:
                s1 += f",{k} {v}"
        print( s1 )
        return False
    return True

def check_connectivity(network):
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
    for c_id, cmpd in network.compounds.items():
        if c_id not in unrefed:
            print(f"Compound {c_id} '{cmpd.name}'is not used in any reactions.")
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
    print(f'Modifiers are :{modifiers}')
    for ring in cycles:
        mods = modifiers.intersection(ring)
        if mods:
            print (f"Modifiers {mods} appear in cycle: {ring}")
        else:
            print (f"No modifiers in cycle: {ring}")
    return graph

def stochiometric_matrix( network ):
    """
    Calculate the stochiometric matrix S for a metabolic network

    S(i,j) is the stochiometry of production of compound(i) by reaction (j),
    thus reagents get negative stochiometries.
    Compound indexes can be found from ??? and reaction indexes from ???
    """
    reaction_keys = list(network.reactions.keys())
    compound_keys = list(network.compounds.keys())
    s_matrix = np.zeros((len(compound_keys),len(reaction_keys)))
    # pylint: disable=C0200
    for j in range(0,len(reaction_keys)):
        reaction = network.reactions[reaction_keys[j]]
        for substrate in reaction.subst:
            val  = - substrate[1]
            i    = compound_keys.index(substrate[0])
            s_matrix[i][j] = val
        for substrate in reaction.prod:
            val  = substrate[1]
            i    = compound_keys.index(substrate[0])
            s_matrix[i][j] = val
    # pylint: enable=C0200
    return s_matrix, reaction_keys, compound_keys

def unit_test():
    """
    Test the functions and structures in this file.
    """
    print("Running the tests...")
    print("Testing Reactions...")
    rxn.unit_test()
    print("Testing Compounds...")
    cpd.unit_test()
    print("Finished testing.")

def main():
    """
    The main routine for the network file
    """
    if len(sys.argv) > 1:
        if sys.argv[1] == 'test':
            unit_test()
            sys.exit()

    print("This file is used to define a network class")
    print("Run with argument 'test' to do tests.")

if __name__ == '__main__':
    main()
