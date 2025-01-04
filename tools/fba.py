"""
This program performs a flux balance analysis on a sbml model file.

The program should read in an sbml file and then analyse the model in order to
perform a flux balance analysis. This analysis involves several steps to ensure
the model and analysis are coherent and self-consistent.

1. Are all cofactors used in the reactions with flux included in the BIOMASS
   reaction at appropriate stochiometries? If not modify the BIOMASS reaction
   appropriately.
2. Are all running reaction cycles filled? Check that at least one intermediate
   in all closed reaction cycle appears at an appropriate stochiometry in the
   BIOMASS reaction? If not make the necessary modifications.
3. Does the fba run producing a growth rate? If not analyse why not by
   isolating parts of biomass and backtracking through metabolism.
4. Return to step 1 until no more modifications need to be made.
5. Add annotations on reaction fluxes from the fba to the model.
6. Output the resulting model (-o filename).
7. Identify and report on unused reactions and species (-r)
8. Identify and report on thermodynamic hurdles in the model (-r)

The different modifications made to the input sbml file should be logged.
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

def usage():
    """
    Print a usage message
    """
    print( f'Usage: {sys.argv[0]} [-rh] sbml_file [-o destination]')
    print(  '-r include report on resulting fluxes.')
    print(  '-h print this usage information.')
    print(  '-o destination output resulting sbml file to destination.')
    sys.exit()

def main():
"""
Main routine to perform the various steps.
"""
# Command line handling
    report = False
    destination = None
    source = None

    for i in range(1,len(sys.argv)):
        item = sys.argv[i]
        if item[0] = '-':
            for j in range(1,len(item)):
                if item[j] == 'h':
                    usage()
                elif item[j] == 'r':
                    report = True
                elif item[j] == 'o':
                    if i+1 == len(sys.argv):
                        usage()
                    destination = sys.argv[i+1]
                    i = i+1
                else:
                    usage()
        else:
            if source:
                usage()
            source = item

    if not source:
        usage()

# 0. Read the sbml file both into a model and for cobra.

    try:
        with open( source, encoding="UTF8") as fp:
            model = BeautifulSoup( fp, 'xml' )
    except:
        print( f"Failed to read sbml model from 'source'.")
        sys.exit()

    cobra_model = read_sbml_model( source )

    while not finished:
# 1. Are all cofactors used in the reactions with flux included in the BIOMASS
#    reaction at appropriate stochiometries? If not modify the BIOMASS reaction
#    appropriately.
# 2. Are all running reaction cycles filled? Check that at least one species
#    in all closed reaction cycle appears at an appropriate stochiometry in the
#    BIOMASS reaction? If not make the necessary modifications.
# 3. Does the fba run producing a growth rate? If not analyse why not by
#    isolating parts of biomass and backtracking through metabolism.
        solution = cobra_model.optimize()
        print( f'Solved fba with objective_value:{solution.objective_value}')
# 4. Return to step 1 until no more modifications need to be made.

#   end of while not finished loop

# 5. Add annotations on reaction fluxes from the fba to the model.

# 6. Output the resulting model (-o filename).
    if destination:
        with open(destination, 'w') as f:
            f.write(model.prettify(formatter="minimal"))

# 7. Identify and report on unused reactions and species (-r)
# 8. Identify and report on thermodynamic hurdles in the model (-r)

def holding():
    cofactors = []
    compounds = []
    energy_flux = []
    if solution.status == 'optimal':
        for item in solution.fluxes.items():
            reaction = model.find('reaction', id=f"R_{item[0]}")
            if reaction is None:
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

    # Search for unmakeable used compounds
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


if __name__ == '__main__':
    main()
