"""
Modify an sbml file by keeping or removing selected reactions.

This program takes 2 or 3 command line options. An optional filename to use as
the source document, if none is given stdin will be used. A flag -r or -k to
remove or keep certain reactions, if no flag is given -k is assumed. Finally
a regular expression that will be tested against the reaction names.

The program
1. Read in the sbml model.
2. Trim the list of reactions according to the regular expression.
3. Trim the list of species to only keep those included in the remaining
   reactions.
4. Output the resulting sbml model to stdout.
"""

# Allow todo comments and dont insist on capitals if initialized with a default.
# pylint: disable='fixme' 'invalid-name' 'bare-except'

import sys
import re
from bs4 import BeautifulSoup

def usage():
    """
    Print a brief usage line, as the program was incorrectly called.
    """
    print(f'Usage: {sys.argv[0]} [sbml_file] [-r|k] regex' )
    sys.exit()

def inreactions_sp(reactions, species_id) :
    """
    Is the species_id refered to in the reaction list?

    References can either be as a speciesReference in reactant or product lists
    or as a modifierSpeciesReference in a listOfModifiers.
    """

    for my_species in reactions.find_all('speciesReference'):
        if my_species['species'] == species_id:
            return True

    for my_species in reactions.find_all('modifierSpeciesReference'):
        if my_species['species'] == species_id:
            return True

    return False

def inreactions_id(reactions, reaction_id ):
    """
    Is the reaction_id included in the reaction list?
    """
    for my_reaction in reactions.find_all('reaction'):
        if my_reaction['id']==reaction_id:
            return True

    return False

# Command line handling

if len(sys.argv) > 4 or len(sys.argv) <= 1:
    usage()

# Defaults
filename = None
fileflag = False
keepflag = True

i = len(sys.argv)-1
regex = re.compile(sys.argv[i])
i = i-1
while i > 0:
    if sys.argv[i][0] == '-':           # Handle flags
        keep = sys.argv[i]
        if keep[1] == 'r':
            keepflag = False
        elif keep[1] != 'k':
            usage()
    else:
        if fileflag:
            usage()
        filename = sys.argv[i]          # Handle filename
        fileflag = True
    i = i-1

# Read the sbml file.

if fileflag:
    try:
        with open( filename, encoding="UTF8") as fp:
            model = BeautifulSoup( fp, 'xml' )
    except:
        print( f"Failed to read sbml model from '{sys.argv[1]}'.")
        sys.exit()
else:
    try:
        model = BeautifulSoup(sys.stdin, 'xml' )
    except:
        print( "Failed to read sbml from stdin.")
        sys.exit()

# Edit the Reactions

for reaction in model.find_all('reaction'):
    if regex.match(reaction['id']) or regex.match(reaction['name']):
        if not keepflag:
            reaction.decompose()
    else:
        if keepflag:
            reaction.decompose()

# Edit the species

for species in model.find_all('species'):
    if not inreactions_sp(model.find('listOfReactions'), species['id']):
        species.decompose()

# Edit the fbc:FluxObjectives

for reaction in model.find_all('fbc:fluxObjective'):
    if not inreactions_id (model.find('listOfReactions'),
                            reaction['fbc:reaction']):
        reaction.decompose()

# if there are no fluxObjective's in the list remove the list!
if not model.find('fbc:fluxObjective'):
    fbc_list = model.find('fbc:listOfObjectives')
    fbc_list.decompose()

# Output the result

print(model.prettify(formatter="minimal"))
