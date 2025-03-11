"""
Check that the reactions in sbml file are all balanced.

Usage:
        check_balance <sbml_file>

This program takes one required argument, the name of an sbml file that will be
read and checked.

The program will give a warning for non integer formulae - this is to be
expected for some "aggregation" metabolites. However unbalanced equations
should just be the objective biomass reaction, and the environment exchange
reactions, these are not checked as they are in model.boundary.

The program relies on the `cobra` python library.
"""

# Allow todo comments and dont insist on capitals if initialized with a default.
# pylint: disable='fixme' 'bare-except'

# TODO: expand to more general sanity check on sbml file.

import sys
import cobra as cb

def main():
    """
    The main function does it all.
    """
    if len(sys.argv) != 2:
        print(f'Usage: {sys.argv[0]} <sbml_file>' )
        sys.exit()

    # Read the sbml file.

    try:
        model = cb.io.read_sbml_model(sys.argv[1])
    except:
        print( f"Failed to read sbml model from '{sys.argv[1]}'.")
        sys.exit()


    print( f"Read sbml model from '{sys.argv[1]}'." )
    print(  "Checking reactions are balanced.")
    for reaction in set(model.reactions) - set(model.boundary):
        if len(reaction.check_mass_balance()) > 0 :
            print(reaction.id + str(reaction.check_mass_balance()) )
    print( "Done.")

if __name__ == '__main__':
    main()
