#!/usr/bin/env python3
"""
This file generates an sbml model based on csv file.

Version 3. Using network classes becomes very simple.
"""

# pylint: disable='fixme' 'invalid-name' 'bare-except'

import argparse
import network

def main():
    """
    Main routine using network object
    """
    # Handle command line if pertinent
    parser = argparse.ArgumentParser(
        description = 'Convert a csv file describing a metabolic network to sbml.'
    )
    parser.add_argument('-o','--outfile', default='-',
        help= "Destination file, '-' means stdout.")
    parser.add_argument('infile', nargs='?', default=['-'],
        help='File to translate where "-" is stdin (,defaults to stdin).')
    args = parser.parse_args()

    network.Network.from_csv(args.infile).write_sbml(args.outfile)

if __name__ == '__main__':
    main()
