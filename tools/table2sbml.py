"""
This file generates an sbml model based on Tables of compounds and reactions.

The table formats:
Reaction table is a tab separated list of:
* ID for reaction if possible from KEGG or another database.
* Name Reaction name
* Reagents list [Compound ID's - usually KEGG or CHEBI]
* Products list [Compound ID's - usually KEGG or CHEBI]
* Modifier list [Compound ID's - usually KEGG or CHEBI]
* geneRule rule using gene product names, and, or and parentheses for parsing.
* EC# for an enzyme doint the reaction (for annotations).
* FreeEnergy A list [dG0, uncertainty] at pH 7.0 (for annotations).
* Dbases A list of dbase references [ dbase, id, ...] (for annotations).

Compounds table is a tab separated list of:
* ID for compound if possible from KEGG or another database.
* Name Compound name
* Compartment (defaults to 'i')
* initialConcentration (defaults to '1.0')
* Charge (defaults to 0)
* ChemicalFormula in (C/H/alphabetical or alphabetical (if no C))
* Dbases A list of dbase references [ dbase, id, ...] (for annotations).

For chemical formulae define:
* Dna (DNA)
* Rna (RNA)
* Pro (Protein)
* Acp (Acyl carrier protein)
"""

# pylint: disable='fixme' 'invalid-name' 'bare-except'
# TODO: define table formats for efficiency and possible completeness.

import sys
import pandas as pd
import numpy as np
from bs4 import BeautifulSoup

COMPARTMENT_COLUMN = 2
GENE_COLUMN = 5

def usage():
    """
    Print usage info to stdout.
    """
    print(f'Usage: {sys.argv[0]} [-o output_file] compounds reactions')
    sys.exit()

def command_line( args ):
    """
    Handle the command line filenames and flags.
    """
    compounds = None
    reactions = None
    output = "-"

    skip = False
    for i in range(1, len(args)):
        if not skip:
            if args[i] == '-o':
                skip = True
                if i+1 >= len(args):
                    usage()
                output = args[i+1]
            else:
                if not compounds:
                    compounds = args[i]
                elif not reactions:
                    reactions = args[i]
                else:
                    usage()

    if not compounds:
        compounds = "-"
    if not reactions:
        reactions = "-"

    return compounds, reactions, output

def read_table( source ):
    """
    Read a tab separated file into a pandas DataFrame
    """
    try:
        with open(source, 'r', encoding="UTF8") if source != '-' else sys.stdin as fp:
            data_table = pd.read_csv ( fp, sep = '\t')
    except:
        print( f"Failed to read table from '{source}'.")
        sys.exit()
    return data_table

def gene_elements( rule_list ):
    """
    Return the list of genes referenced in a list of gene rules.
    """
    elements = []
    for gene_rule in rule_list:
        if not isinstance( gene_rule, float ):
            elements.extend(gene_rule.split())
    return list(set(elements) - set(['and','or','AND','OR']))

def add_listOfUnitDefinitions( model, units_definitions ):
    pass

def add_listOfCompartments( model, compartments ):
    pass

def add_listOfSpecies( model, species ):
    pass

def add_listOfReactions( model, reactions ):
    pass

def add_listOfObjectives( model, objectives ):
    pass

def add_listOfGenes( model, genes ):
    pass

def main():
    """
    Main routine for table2sbml.

    Read in 2 data tables and organize an sbml file from them.
    """
    compounds_file, reactions_file, output_file = command_line( sys.argv )

    compounds = read_table( compounds_file )
    print(compounds)
    try:
        compartments = list(set(compounds['Compartment'].tolist()))
    except KeyError:
        compartments = None
    print(compartments)
    reactions = read_table( reactions_file )
    print(reactions)
    try:
        genes = gene_elements(reactions['geneRule'].tolist())
    except KeyError:
        genes = None
    print(genes)
    objectives = None
    units_definitions = ["mmol_per_gDW_per_hr"]

    model = BeautifulSoup()
    add_listOfUnitDefinitions( model, units_definitions )
    add_listOfCompartments( model, compartments )
    add_listOfSpecies( model, compounds )
    add_listOfReactions( model, reactions )
    add_listOfObjectives( model, objectives )
    add_listOfGenes( model, genes )

    with open(output_file, 'w', encoding="UTF8") if output_file != '-' else sys.stdout as fp:
        fp.write(model.prettify(formatter="minimal"))

if __name__ == '__main__':
    main()
