"""
This program takes two arguments. First an sbml file that is going to be read and edited,
and second a tsv file containing a table linking reaction names to cofactor codes. The
program should for each line in the table:
1. If all the species do not already exist add new new species to the SBML
2. If necessary create a list of modulators for the indicated reaction.
3. Set the modulators of the reaction to correspond to the indicated species.
"""

# Allow todo comments and dont insist on capitals if initialized with a default.
# pylint: disable='fixme' 'invalid-name' 'bare-except'

import sys
from bs4 import BeautifulSoup

if len(sys.argv) != 3:
    print(f'Usage: {sys.argv[0]} sbml_file cofactor_table' )
    sys.exit()

# Read the sbml file.

try:
    with open( sys.argv[1], encoding="UTF8") as fp:
        model = BeautifulSoup( fp, 'xml' )
except:
    print( f"Failed to read sbml model from '{sys.argv[1]}'.")
    sys.exit()

# Read the cofactor table.

my_table = []
i = 0
try:
    with open( sys.argv[2], encoding="UTF8") as fp:
        for line in fp:
            elements = line.split('\t')
            if i > 0 and len(elements[1].strip()) > 0:
                my_table.append(( elements[0], elements[1].strip() ))
            i += 1
except:
    print( f"Failed to read cofactor table from '{sys.argv[2]}'." )
    sys.exit()

# Now process the cofactors.
# 1. Check the species are in the sbml.
species_list = []
for element in my_table:
    for species in element[1].split(','):
        if species.strip() not in species_list :
            species_list.append(species.strip())

for species in model.find_all('species'):
    # pylint: disable='cell-var-from-loop'
    specid = species['id']
    if specid[4:] in species_list:
        species_list[:] = list(filter(lambda a: a!=specid[4:], species_list ))

if len(species_list) > 0:
    print( 'The following species need to be added to the sbml model.' )
    for species in species_list:
        print( species )
    sys.exit()

# 2. Modify reactions as necessary the list of modulators for each reaction

my_dict = dict(my_table)

for reaction in model.find_all('reaction'):
    reacid = reaction['id']
    if reacid in my_dict:
        modifier   = []
        modlist = reaction.find('listOfModifiers')
        if modlist is None :
            prodlist = reaction.find('listOfProducts')
            assert prodlist is not None
            modlist = model.new_tag('listOfModifiers')
            prodlist.insert_after(modlist)
        for cofactor in my_dict[reacid].split(','):
            cof_tag = model.new_tag('modifierSpeciesReference')
            cof_tag['species'] = f'M_{cofactor.strip()}' 
            modlist.append( cof_tag )
            pass

print(model.prettify(formatter="minimal"))
