"""
This file generates an sbml model based on Tables of compounds and reactions.

The table formats:
Reaction table is a tab separated list of:
* Reaction ID if possible from KEGG or another database.
* Reaction name
* Reagents list [Compound ID's - usually KEGG or CHEBI]
* Products list [Compound ID's - usually KEGG or CHEBI]
* Modifier list [Compound ID's - usually KEGG or CHEBI]
* Gene product rule using gene product names, and, or and parentheses for parsing.
* The EC# for an enzyme doint the reaction (for annotations).
* A list [dG0, uncertainty] at pH 7.0 (for annotations).
* A list of dbase references [ dbase, id, ...] (for annotations).

Compounds table is a tab separated list of:
* Compound ID if possible from KEGG or another database.
* Compound name
* Compartment (defaults to 'i')
* initialConcentration (defaults to '1.0')
* Charge (defaults to 0)
* ChemicalFormula in (C/H/alphabetical or alphabetical (if no C))
* A list of dbase references [ dbase, id, ...] (for annotations).

For chemical formulae define:
* Dna (DNA)
* Rna (RNA)
* Pro (Protein)
* Acp (Acyl carrier protein)
"""

# pylint: disable='fixme' 'invalid-name' 'bare-except'
# TODO: define table formats for efficiency and possible completeness.

import sys
from bs4 import BeautifulSoup

def main():
    compounds_file, reactions_file, output_file = command_line( sys.argv )

    compounds = read_compounds( compounds_file )
    compartments = find_compartments( compounds )
    reactions = read_reactions( reactions_file )
    genes = find_genes( reactions )
    objectives = None
    units_definitions = ["mmol_per_gDW_per_hr"]

    model = BeautifulSoup()
    add_listOfUnitDefinitions( model, units_definitions )
    add_listOfCompartments( model, compartments )
    add_listOfSpecies( model, species )
    add_listOfReactions( model, reactions )
    add_listOfObjectives( model, objectives )
    add_listOfGenes( model, genes )

    with open(output_file) as fp:
        write( fp, model )

if __name__ == '__main__':
    main()
