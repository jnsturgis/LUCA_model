"""
Check that the reactions in sbml file are all balanced.

This program takes one argument, an sbml file that is going to be read
and checked.
This will give a warning for non integer formulae - this is expected for the
some "aggregation" metabolites. However unbalanced equations should just be
the objective, and the exchange reactions, these are not checked as they are
in model.boundary
"""

# Allow todo comments and dont insist on capitals if initialized with a default.
# pylint: disable='fixme' 'bare-except'

# TODO: expand to more general sanity check on sbml file.

import sys
import cobra as cb

if len(sys.argv) != 2:
    print(f'Usage: {sys.argv[0]} sbml_file' )
    sys.exit()

# Read the sbml file.

try:
    model = cb.io.read_sbml_model("LUCA.sbml")
except:
    print( f"Failed to read sbml model from '{sys.argv[1]}'.")
    sys.exit()


print( f"Read sbml model from '{sys.argv[1]}'." )
print(  "Checking reactions are balanced.")
for reaction in set(model.reactions) - set(model.boundary):
    if len(reaction.check_mass_balance()) > 0 :
        print(reaction.id + str(reaction.check_mass_balance()) )
print( "Done.")
