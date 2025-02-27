"""
This file defines some routines for handling pathway tables, and checks the table given as stdin.

"""

import network

def main():
    """
    Verify the table provided on input as a pathway table.

    Read the input as a pathway file and perform the following checks:
    1. All reagents in reactions appear as compounds in the file.
    2. All reactions in the file are equilibrated (same numbers of atoms in and out)
    3. Analyse the graph for products, substrates and recycled molecules
        Warn if no elements a cycle feature as modifiers

    TODO for cycles need to remove overcommon substrates as they make lots of little
    cycles.
    """
    source = "-"
    # For connectivity and cycle detection remove common metabolites
    exclude = ['Mi_C00001', 'Mi_C00002', 'Mi_C00003', 'Mi_C00004', 'Mi_C00005',
        'Mi_C00080', 'Mi_C00011', 'Mi_C00006', 'Mi_C00008', 'Mi_C00009', 'Mi_C00013',
        'Mi_C00020']
#    handle_argv()
    graph = network.Network.from_csv(source)
    network.check_connectivity(graph)
    network.check_balance(graph)
    network.analyse_pathway(graph, exclude)
    network.stochiometric_matrix( graph )
    network.reaction_adjacency_graph( graph )

if __name__ == '__main__':
    main()
