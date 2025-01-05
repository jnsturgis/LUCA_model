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

# TODO: setup routines that operate on the BeautifulSoup structure.
# Less efficient but more portable than current version.

import sys
import math

from cobra.io import read_sbml_model
from bs4 import BeautifulSoup

# Constants

R_BIOMASS = 'R_BIOMASS'

# Utility functions

def sign(x : float ) -> float:
    """Implement a sign function, returns +1 or -1 depending on sign of x"""
    return math.copysign(1, x)

# Program control and user interface

def usage():
    """
    Print a usage message
    """
    print( f'Usage: {sys.argv[0]} [-rh] sbml_file [-o destination]')
    print(  '-r include report on resulting fluxes.')
    print(  '-h print this usage information.')
    print(  '-o destination output resulting sbml file to destination.')
    sys.exit()

def command_line( args ):
    """
    Handle the command line options returning the report flag and filenames.
    """
    my_report = False
    infile = None
    outfile = None
    skip = False

    for i in range(1,len(args)):
        if skip:
            skip = False
            continue
        item = args[i]
        if item[0] == '-':
            for j in range(1,len(item)):
                if item[j] == 'h':
                    usage()
                elif item[j] == 'r':
                    my_report = True
                elif item[j] == 'o':
                    if i+1 == len(args):
                        usage()
                    outfile = args[i+1]
                    skip = True
                else:
                    usage()
        else:
            if infile:
                usage()
            infile = item
    if not infile:
        usage()
    return my_report, infile, outfile

# Routines that manipulate the sbml model

def update_annotations( model, solution ):
    """
    Update annotations for each reaction to include the fba-flux.
    """
    for item in solution.fluxes.items():
        reaction = model.find('reaction', id=f"R_{item[0]}")
        if reaction is None:
            print( f'Unable to find reaction R_{item[0]}.')
            sys.exit()
        old_annotation = None
        for annotation in reaction.find_all('rdf:li'):
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

def exchange_species( model ):
    """
    Return a list of species that are exchanged with the environment.

    Go through the list of reactions collecting species that appear in a
    reaction with no substrates or no products.
    """
    exchanged = []
    for reaction in model.find_all('reaction'):
        if not reaction['id'] == R_BIOMASS:
            if ((not reaction.find('listOfReactants')) or
                (not reaction.find('listOfProducts'))):
                for species in reaction.find_all('speciesReference'):
                    exchanged.append(species['species'])
    return list(set(exchanged))

def collect_cofactors( model, solution ):
    """
    Collect the cofactors (modifiers) used in reactions with flux.
    """
    cofactors = []
#   compounds = []
    for item in solution.fluxes.items():
        if float(item[1]) != 0.0 :
            reaction = model.find('reaction', id=f"R_{item[0]}")
            if reaction is None:
                print( f'Unable to find reaction R_{item[0]}.')
                sys.exit()
            for species in reaction.find_all('modifierSpeciesReference'):
                cofactors.append(species['species'])
    return list(set(cofactors))

# Reports on the model

def thermodynamic_report(model, solution):
    """
    Print a report on the thermodynamic properties of the fba solution.
    """
    energy_flux = []
    for item in solution.fluxes.items():
        if float(item[1]) != 0.0 :
            reaction = model.find('reaction', id=f"R_{item[0]}")
            if reaction is None:
                print( f'Unable to find reaction R_{item[0]}.')
                sys.exit()

            r_energy = 0.0
            for annotation in reaction.find_all('rdf:li'):
                if annotation['rdf:resource'].split('/')[3] == 'dG0' :
                    r_energy = float(annotation['rdf:resource'].split('/')[4])
                    break
            e_flux    = r_energy * float(item[1])
            r_energy *= sign(float(item[1]))
            energy_flux.append((f'R_{item[0]}',e_flux, r_energy ))
    i=0
    energy_flux.sort(key=lambda x: x[2], reverse=True)
    for i, e_flux in enumerate(energy_flux):
        if e_flux[2] <= 10:
            break
    #pylint: disable='undefined-loop-variable'
    energy_flux = energy_flux[:i]
    energy_flux.sort(key=lambda x: x[1], reverse=True)

    print( f"The solution uses {len(energy_flux)} reactions with dG0 > 10")
    print( "Reaction       flux*dG0       dG0")
    for item in energy_flux:
        print( f'{item[0]:<11}  {item[1]:>10.3f} {item[2]:>10.1f}')

def unused_report( model, solution ):
    """"
    Make a report on the unused reactions and species.
    """
    # TODO: Collect list of unused reactions and unfeatured species.
    print("The following reactions are not used in the fba solution:")
    print("The following species do not feature in the fba solution:")

def main():
    """
    Main routine to perform the various steps.
    """
    # Command line handling
    report, source, destination = command_line(sys.argv)

# 0. Read the sbml file both into a model and for cobra.

    try:
        with open( source, encoding="UTF8") as fp:
            model = BeautifulSoup( fp, 'xml' )
    except:
        print( f"Failed to read sbml model from '{source}'.")
        sys.exit()

    cobra_model = read_sbml_model( source )

    finished = False
    solution = cobra_model.optimize()
    print( f'Initial fba with objective_value:{solution.objective_value}')

    while not finished:
# 1. Are all cofactors used in the reactions with flux included in the BIOMASS
#    reaction at appropriate stoichiometry? If not modify the BIOMASS reaction
#    appropriately.
        if solution.status == 'optimal':
            cofactor_list = collect_cofactors( model, solution )
            print( "Used cofactors")
            print( cofactor_list )
            reaction = model.find('reaction', id = R_BIOMASS )
            # pylint: disable='undefined-loop-variable'
            assert reaction is not None
            for species in reaction.find_all('speciesReference'):
                try:
                    cofactor_list.remove( species['species'] )
                except ValueError:
                    pass

            # If necessary modify the BIOMASS reaction adding the lines to the
            # end of listOfReactants
            # <speciesReference constant="true" species="M_i_Fe2S2" stoichiometry="1"/>
            # setup another attempt.

            # DONE: remove from list if provided by environment
            cofactor_list_too = cofactor_list.copy()
            for item in cofactor_list_too:
                if item in exchange_species(model):
                    cofactor_list.remove(item)

            # TODO: recover stoichiometry from species initial concentration

            if len(cofactor_list) > 0:
                print( f"These {len(cofactor_list)} cofactors are required for this solution." )
                print( cofactor_list )

            if len(cofactor_list) > 0:
                for item in cofactor_list:
                    new_species = model.new_tag('speciesReference')
                    new_species['constant'] = "true"
                    new_species['species']  = f'{item}'
                    new_species['stoichiometry']  = "1"
                    species.insert_after(new_species)

                print( "These have been added to the BIOMASS reaction." )
                print( cofactor_list )
                finished = False

# 2. Are all running reaction cycles filled? Check that at least one species
#    in all closed reaction cycle appears at an appropriate stochiometry in the
#    BIOMASS reaction? If not make the necessary modifications.
        if solution.status == 'optimal':
            # TODO: Search for cycles
            # TODO: Ensure cycles are filled in BIOMASS reaction
            pass
# 3. Does the fba run producing a growth rate? If not analyse why not by
#    isolating parts of biomass and backtracking through metabolism.
        if solution.status != 'optimal':
            pass

# 4. Return to step 1 until no more modifications need to be made.

        if not finished:
            # TODO: make a new cobra model from the BeautifulSoup
            cobra_model = read_sbml_model( model.prettify(formatter="minimal"))
            solution = cobra_model.optimize()
            print( f'Solved fba with objective_value:{solution.objective_value}')
            finished = True

#   end of while not finished loop

# 5. Add annotations on reaction fluxes from the fba to the model.
    if solution.status == 'optimal':
        update_annotations( model, solution )

# 6. Output the resulting model (-o filename).
    if destination:
        with open(destination, 'w', encoding="UTF8") as f:
            f.write(model.prettify(formatter="minimal"))

# 7. Identify and report on unused reactions and species (-r)
# 8. Identify and report on thermodynamic hurdles in the model (-r)
    if report:
        unused_report( model, solution)
        thermodynamic_report( model, solution )

# Call the main routine
if __name__ == '__main__':
    main()
