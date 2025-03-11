#!/usr/bin/env python3
"""
This file defines some routines for handling pathway tables, and checks the table given as stdin.
"""

import argparse
import network

def main():
    """
    Verify the table provided on input as a pathway table.

    Read the input as a pathway file and perform the following checks:
    1. All reagents in reactions appear as compounds in the file.
    2. All reactions in the file are equilibrated (same numbers of atoms in
       and out)
    3. Analyse the graph for products, substrates and recycled molecules
       Warn if no elements a cycle feature as modifiers

    TODO for cycles need to remove overcommon substrates as they make lots of little
    cycles.
    """

    parser = argparse.ArgumentParser(
        prog = 'csv_table_verify.py',
        description = 'Analyse a csv file describing a metabolic network.'
    )
    parser.add_argument('-b','--balance',
        action='store_true',
        help='Check the equations are balanced')
    parser.add_argument('-c','--connectivity',
        action='store_true',
        help='Analyse graph connectivity')
    parser.add_argument('-v','--verbose',
        action='store_true',
        help='Provide extra output about programme execution')
    parser.add_argument('-r','--rings',
        help= "String containing metabolite id's to ignore in cycle analysis.")
    parser.add_argument('filenames', nargs='*', default=['-'],
        help='List of files to analyse where "-" is stdin (,defaults to stdin).')
    args = parser.parse_args()

    for filename in args.filenames:
        if len(args.filenames) > 1:
            print(f"Analysing {filename}.")
        graph = network.Network.from_csv(filename)
        if args.connectivity:
            if args.verbose:
                print('Analysing connectivity of graph.')
            network.check_connectivity(graph)
        if args.balance:
            if args.verbose:
                print('Checking reaction balance.')
            network.check_balance(graph)
        if args.rings:
            if args.verbose:
                print('Analysing cycles in the graph.')
            network.analyse_pathway(graph, args.rings.split())

if __name__ == '__main__':
    main()
