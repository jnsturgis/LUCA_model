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
import math

from cobra.io import read_sbml_model
from bs4 import BeautifulSoup

def sign(x : float ) -> float:
    """Implement a sign function, returns +1 or -1 depending on sign of x"""
    return math.copysign(1, x)

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

error = 0

# Do the fba

solution = cobra_model.optimize()
print( f'Solved fba with objective_value:{solution.objective_value}')
cofactors = []
compounds = []

energy_flux = []
if solution.status == 'optimal':
    for item in solution.fluxes.items():
        reaction = model.find('reaction', id=f"R_{item[0]}")
        if reaction == None:
            print( f'Unable to find reaction R_{item[0]}.')
            sys.exit()
        # Add required cofactors to list of cofactors
        for species in reaction.find_all('modifierSpeciesReference'):
            cofactors.append(species['species'])
        # Collect all used and produced reagents
        for species in reaction.find_all('speciesReference'):
            compounds.append(species['species'])
        # Add flux to annotations (replacing old one if present)
        old_annotation = None
        r_energy = 0.0
        for annotation in reaction.find_all('rdf:li'):
            if annotation['rdf:resource'].split('/')[3] == 'dG0' :
                r_energy = float(annotation['rdf:resource'].split('/')[4])
            if annotation['rdf:resource'].split('/')[3] == 'fba-flux' :
                old_annotation = annotation
        new_annotation = model.new_tag('rdf:li')
        new_annotation['rdf:resource'] = f'https://identifiers.org/fba-flux/{item[1]}'
        # pylint: disable='undefined-loop-variable'
        if annotation is None:
            pass                                 # Should insert first annotation here
        else:
            annotation.insert_after( new_annotation )
        if old_annotation is not None:
            old_annotation.decompose()
        # Add dG0.flux to energy_flux in dataframe
        e_flux    = r_energy * float(item[1])
        r_energy *= sign(float(item[1]))
        energy_flux.append((f'R_{item[0]}',e_flux, r_energy ))

    # Finished looping over reactions and fluxes. Analyse the resulting lists and model.
    # Deal with required cofactors
    cofactor_list = list(set(cofactors))
    compound_list = list(set(compounds))

    reaction = model.find('reaction', id='R_BIOMASS' )
    # pylint: disable='undefined-loop-variable'
    assert reaction is not None
    for species in reaction.find_all('speciesReference'):
        try:
            cofactor_list.remove( species['species'] )
        except ValueError:
            pass
    if len(cofactor_list) > 0:
        print( f"These {len(cofactor_list)} cofactors are required for this solution." )
        print( cofactor_list )
        error += 1

    # Energetic problems reactions of more than 4kbT positive and a reasonable flux
    energy_flux.sort(key=lambda x: x[2], reverse=True)
    for i, e_flux in enumerate(energy_flux):
        if e_flux[2] <= 10:
            end = i
            break
    energy_flux = energy_flux[:i]

    energy_flux.sort(key=lambda x: x[1], reverse=True)
    for i, e_flux in enumerate(energy_flux):
        if e_flux[1] <= 10:
            end = i
            break
    energy_flux = energy_flux[:i]
    print( energy_flux )

    # Search for unmakable used compounds
    # A subset of species that for all reactions with flux have coefficients, after
    # ajusting coefficients for flux direction, that sum to 0.
    # How to search efficiently
    if len(compound_list) > 0:
        print( f"These {len(compound_list)} compounds feature in the flux.")
#        print( compound_list )
#        error += 1

# Oops solver failed
else:
    error += 1
    print( "Failed to find optimal solution.")

# If consistent then save everything
if error == 0:
    print(model.prettify(formatter="minimal"))
