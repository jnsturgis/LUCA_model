"""
This program is desigend to perform a flux balance analysis and analyse the results
using them to examine the thermodynamics and ensure that the model BIOMASS function
is complete.

Step 1. Read the sbml file and do the fba.
Step 2. Checks and information
    2.1 Are all cofactors used in the reactions with flux in the BIOMASS?
    2.2 Are all running reaction cycles filled?
    2.3 What unfavourable reactions have high flux?
    2.4 Report on any unused reactions
"""

# Allow todo comments and dont insist on capitals if initialized with a default.
# pylint: disable='fixme' 'invalid-name' 'bare-except'

import sys
import cobra
from cobra.io import read_sbml_model
from bs4 import BeautifulSoup

if len(sys.argv) != 2:
    print(f'Usage: {sys.argv[0]} sbml_file' )
    sys.exit()

# Read the sbml file both into a model and for cobra.

source = sys.argv[1]
try:
    with open( source, encoding="UTF8") as fp:
        model = BeautifulSoup( fp, 'xml' )
except:
    print( f"Failed to read sbml model from '{sys.argv[1]}'.")
    sys.exit()

cobra_model = read_sbml_model( source )

# Do the fba

solution = cobra_model.optimize()
print( f'Solved fba with objective_value:{solution.objective_value}')
cofactors = []
energy_flux = []
if solution.status == 'optimal':
    for item in solution.fluxes.items():
        reaction = model.find('reaction', id=f"R_{item[0]}")
        # Add cofactors to list of cofactors
        for species in reaction.find_all('modifierSpeciesReference'):
            cofactors.append(species['species'])
        # Add flux to annotations (replacing old one if present)
        old_annotation = None
        e_flux = 0.0
        for annotation in reaction.find_all('rdf:li'):
            if annotation['rdf:resource'].split('/')[3] == 'dG0' :
                e_flux = float(annotation['rdf:resource'].split('/')[4]) * float( item[1] )
            if annotation['rdf:resource'].split('/')[3] == 'fba-flux' :
                old_annotation = annotation
        new_annotation = model.new_tag('rdf:li')
        new_annotation['rdf:resource'] = f'https://identifiers.org/fba-flux/{item[1]}'
        annotation.insert_after( new_annotation )
        if old_annotation is not None:
            old_annotation.decompose()
        # Add dG0.flux to energy_flux in dataframe
        energy_flux.append(e_flux)
    
    # Deal with cofactors    
    print( "These cofactors are required for this solution." )
    print( list(set(cofactors)) )
    # Search for unmakable used compounds

print(model.prettify(formatter="minimal"))
