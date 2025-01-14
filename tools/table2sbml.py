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
# TODO: make robust against missing columns and data
# TODO: avoid double entry

import sys
import pandas as pd
from bs4 import BeautifulSoup

SBML_INITIAL = """<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" xmlns:fbc="http://www.sbml.org/sbml/level3/version1/fbc/version2" sboTerm="SBO:0000624" level="3" version="1" fbc:required="false">
  <model metaid="meta_model_name" id="model_name" fbc:strict="true">
  </model>
</sbml>
"""

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

def add_annotation( soup, item ):
    """
    Add to 'item' the xml bits to hold annotation key:value pairs
    """
    new_tag = soup.new_tag('annotation')
    item.append(new_tag)
    new_tag.append(soup.new_tag('rdf:RDF', attrs={
        "xmlns:bqbiol":"http://biomodels.net/biology-qualifiers/",
        "xmlns:bqmodel":"http://biomodels.net/model-qualifiers/",
        "xmlns:dcterms":"http://purl.org/dc/terms/",
        "xmlns:rdf":"http://www.w3.org/1999/02/22-rdf-syntax-ns#",
        "xmlns:vCard":"http://www.w3.org/2001/vcard-rdf/3.0#",
        "xmlns:vCard4":"http://www.w3.org/2006/vcard/ns#"}))
    new_tag = new_tag.find('rdf:RDF')
    new_tag.append(soup.new_tag('rdf:Description', attrs={'rdf:about':f'#{item["metaid"]}'}))
    new_tag = new_tag.find('rdf:Description')
    new_tag.append(soup.new_tag('bqbiol:is'))
    new_tag = new_tag.find('bqbiol:is')
    new_tag.append(soup.new_tag('rdf:Bag'))

def add_annotations( soup, item, item_list ):
    """
    add key value pairs from the 'item_list' to the annotation of 'item'
    """
    new_tag = item.find('rdf:Bag')
    for label, value in zip(item_list.split()[::2], item_list.split()[1::2]):
        new_tag.append(soup.new_tag('rdf:li',attrs={
            'rdf:resource': f'https://identifiers.org/{label}/{value}'
        }))

def add_listOfUnitDefinitions( model, units_definitions ):
    """
    Add to the model a list of unit definitions as necessary.
    """
    # TODO: This needs parsing of units expression (read definition)
    # if necessary add a listOfUnitDefinitions
    if units_definitions:
        for term in units_definitions:
            # if necessary add the unitDefinition to the list
            pass

def add_listOfCompartments( soup, compartments ):
    """
    Add to the model a list of compartments as necessary.
    """
    model = soup.find("model")
    mylist = model.find("listOfCompartments")
    if not mylist:
        new_tag = soup.new_tag("listOfCompartments",)
        model.append( new_tag )
        mylist = new_tag
    if compartments:
        for compartment in compartments:
            old = mylist.find("compartment", attrs={"id":compartment})
            if not old:
                new_tag = soup.new_tag( "compartment", attrs={"id":compartment}, constant="true")
                mylist.append(new_tag)

def add_listOfSpecies( soup, species ):
    """
    Add to the model a list of species as necessary from the DataFrame.
    """
    model = soup.find('model')
    mylist = model.find('listOfSpecies')
    if not mylist:
        new_tag = soup.new_tag('listOfSpecies')
        model.append(new_tag)
        mylist = new_tag
    for row in species.itertuples():
        new_species = soup.new_tag('species',
            attrs={
                "boundaryCondition":"false",
                "compartment":row.Compartment,
                "constant":"false",
                "fbc:charge":row.Charge,
                "fbc:chemicalFormula":row.ChemicalFormula,
                "hasOnlySubstanceUnits":"false",
                "id":row.Id,
                "initialConcentration":row.InitialConcentration,
                "metaid":f'meta_{row.Id}',
                "name":row.Name
            } )
        mylist.append(new_species)
        if len(row.Dbases) > 0:
            add_annotation( soup, new_species)
            add_annotations( soup, new_species, row.Dbases )

def add_listOfReactions( soup, reactions ):
    """
    Add to the model a list of reactions as necessary from the DataFrame.
    """
    model = soup.find('model')
    mylist = model.find('listOfReactions')
    if not mylist:
        new_tag = soup.new_tag('listOfReactions')
        model.append(new_tag)
        mylist = new_tag
    for row in reactions.itertuples():
        new_reaction = soup.new_tag('reaction',
            attrs={
                "fast":"false",                          # Required
                "reversible":"true",                     # Required
                "id":row.Id,                             # Required
                "metaid":f'meta_{row.Id}',               # Optional
                "name":row.Name,                         # Optional
                "fbc:lowerFluxBound":"cobra_default_lb", # Required fbc:strict
                "fbc:upperFluxBound":"cobra_default_ub"  # Required fbc:strict
            } )
        mylist.append(new_reaction)

        if len(row.Reagents):
            new_tag=soup.new_tag('listOfReactants')
            for species in row.Reagents.split():
                new_tag.append(soup.new_tag('speciesReference', attrs={
                    "constant":"true",
                    "species":species,
                    "stoichiometry":"1"
                }))
            new_reaction.append(new_tag)

        if len(row.Products):
            new_tag=soup.new_tag('listOfProducts')
            for species in row.Products.split():
                new_tag.append(soup.new_tag('speciesReference', attrs={
                    "constant":"true",
                    "species":species,
                    "stoichiometry":"1"
                }))
            new_reaction.append(new_tag)

        if len(row.Modifiers):
            new_tag=soup.new_tag('listOfModifiers')
            for species in row.Modifiers.split():
                new_tag.append(soup.new_tag('modifierSpeciesReference', attrs={
                    "species":species,
                }))
            new_reaction.append(new_tag)

        if len(row.GeneRules) > 0:
            # TODO: This needs parsing of gene rules expression (read definition)
            pass

        annotations = row.Dbases
        if row.FreeEnergy:
            temp = row.FreeEnergy.split()
            annotations = f'dG0 {temp[0]} dG0_uncertainty {temp[1]} {annotations}'
        if row.EC:
            annotations=f'ec-code {row.EC} {annotations}'
        if len(annotations) > 0:
            add_annotation( soup, new_reaction)
            add_annotations( soup, new_reaction, annotations )

def add_listOfParameters( soup, parameters ):
    """
    Add to the model a list of parameters.
    """
    model = soup.find('model')
    mylist = model.find('listOfParameters')
    if not mylist:
        new_tag = soup.new_tag('listOfParameters')
        model.append(new_tag)
        mylist = new_tag
    for label, value in zip(parameters[::2], parameters[1::2]):
        new_tag.append(soup.new_tag('parameter',attrs={
            'constant': 'true',
            'id': label,
            'sboTerm': 'SBO:0000626',
            'value': value
        }))

def add_listOfObjectives( soup, objectives ):
    """
    Add to the model a list of objectives as necessary.
    """
    # TODO: Do we really need this, and how is it parsed in the list? (read documentation)
    if objectives:
        pass

def add_listOfGenes( soup, genes ):
    """
    Add to the model a list of genes as necessary.
    """
    if genes:
        for gene in genes:
            pass

def main():
    """
    Main routine for table2sbml.

    Read in 2 data tables and organize an sbml file from them.
    """
    compounds_file, reactions_file, output_file = command_line( sys.argv )

    compounds = read_table( compounds_file )
    try:
        compartments = list(set(compounds['Compartment'].tolist()))
    except KeyError:
        compartments = None
    reactions = read_table( reactions_file )
    try:
        genes = gene_elements(reactions['geneRule'].tolist())
    except KeyError:
        genes = None

    parameters = ["cobra_default_lb", "-1000", "cobra_default_ub", "1000"]

    soup = BeautifulSoup( SBML_INITIAL ,"xml")
    add_listOfUnitDefinitions( soup, ["mmol_per_gDW_per_hr"] )
    add_listOfCompartments( soup, compartments )
    add_listOfSpecies( soup, compounds )
    add_listOfParameters(soup, parameters)
    add_listOfReactions( soup, reactions )
    add_listOfObjectives( soup, None )
    add_listOfGenes( soup, genes )

    with open(output_file, 'w', encoding="UTF8") if output_file != '-' else sys.stdout as fp:
        fp.write(soup.prettify(formatter="minimal"))

if __name__ == '__main__':
    main()
