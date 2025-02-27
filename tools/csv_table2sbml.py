"""
This file generates an sbml model based on tables of compounds and reactions.

Version 3. Using network, reaction and compound classes
"""

# pylint: disable='fixme' 'invalid-name' 'bare-except'

import sys

from bs4 import BeautifulSoup

import network

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
    source = "-"
    output = "-"

    skip = False
    for i in range(1, len(args)):
        if not skip:
            if args[i] == '-o':
                skip = True
                if i+1 >= len(args):
                    usage()
                output = args[i+1]
            elif args[i] == '-h':
                usage()
            else:
                if not source:
                    source = args[i]
                else:
                    usage()

    return source, output

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
            # Dummy
            model += term

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

def add_listOfSpecies( soup, metab_model ):
    """
    Add to the model a list of species as necessary from the DataFrame.
    """
    model = soup.find('model')
    mylist = model.find('listOfSpecies')
    if not mylist:
        new_tag = soup.new_tag('listOfSpecies')
        model.append(new_tag)
        mylist = new_tag
    for key, species in metab_model.compounds.items():
        new_species = soup.new_tag('species',
            attrs={
                "boundaryCondition":"false",
                "compartment":species.compartment,
                "constant":"false",
                "fbc:charge":int(species.charge),
                "fbc:chemicalFormula":species.formula,
                "hasOnlySubstanceUnits":"false",
                "id":key,
                "initialConcentration":species.concentration,
                "metaid":f'meta_{key}',
                "name":species.name
            } )
        mylist.append(new_species)
#        if not pd.isna(row.Dbases) and len(row.Dbases) > 0:
#            add_annotation( soup, new_species)
#            add_annotations( soup, new_species, row.Dbases )

def add_listOfReactions( soup, matab_model ):
    """
    Add to the model a list of reactions as necessary from the DataFrame.
    """
    model = soup.find('model')
    mylist = model.find('listOfReactions')
    if not mylist:
        new_tag = soup.new_tag('listOfReactions')
        model.append(new_tag)
        mylist = new_tag
    for key, reaction in matab_model.reactions.items():
        new_reaction = soup.new_tag('reaction',
            attrs={
                "fast":"false",                          # Required
                "reversible":"true",                     # Required
                "id": key,                               # Required
                "metaid":f'meta_{key}',                  # Optional
                "name":reaction.name,                    # Optional
                "fbc:lowerFluxBound":"cobra_default_lb", # Required fbc:strict
                "fbc:upperFluxBound":"cobra_default_ub"  # Required fbc:strict
            } )
        mylist.append(new_reaction)

        if len(reaction.subst):
            new_tag=soup.new_tag('listOfReactants')
            for item in reaction.subst:
                new_tag.append(soup.new_tag('speciesReference', attrs={
                    "constant":"true",
                    "species":item[0],
                    "stoichiometry":item[1]
                }))
            new_reaction.append(new_tag)

        if len(reaction.prod):
            new_tag=soup.new_tag('listOfProducts')
            for item in reaction.prod:
                new_tag.append(soup.new_tag('speciesReference', attrs={
                    "constant":"true",
                    "species":item[0],
                    "stoichiometry":item[1]
                }))
            new_reaction.append(new_tag)

        if len(reaction.modif):
            new_tag=soup.new_tag('listOfModifiers')
            for species in reaction.modif:
                new_tag.append(soup.new_tag('modifierSpeciesReference', attrs={
                    "species":species,
                }))
            new_reaction.append(new_tag)

#        if not pd.isna(row.GeneRules) and len(row.GeneRules) > 0:
            # TODO: This needs parsing of gene rules expression (read definition)
#            pass

        # TODO this needs better parsing to include all possibilities and add them.
#        annotations = ""
#        if not pd.isna(row.FreeEnergy) and row.FreeEnergy:
#            temp = row.FreeEnergy.split()
#            annotations = f'dG0 {temp[0]} dG0_uncertainty {temp[1]} {annotations}'
#        if not pd.isna(row.EC) and row.EC:
#            annotations=f'ec-code {row.EC} {annotations}'
#        if len(annotations) > 0:
#            add_annotation( soup, new_reaction)
#            add_annotations( soup, new_reaction, annotations )

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
    # Dummy
    if objectives:
        soup += objectives

def add_listOfGenes( soup, genes ):
    """
    Add to the model a list of genes as necessary.

    Dummy
    """
    if genes:
        for gene in genes:
            soup += gene

def main():
    """
    Main routine for table2sbml.

    Read in data and organize an sbml file from them.
    """
    source_file, output_file = command_line( sys.argv )

    data = network.Network.from_csv(source_file)

    parameters = ["cobra_default_lb", "-1000", "cobra_default_ub", "1000"]

    soup = BeautifulSoup( SBML_INITIAL ,"xml")

    add_listOfUnitDefinitions( soup, ["mmol_per_gDW_per_hr"] )
    add_listOfCompartments( soup, network.find_compartments( data ) )
    add_listOfSpecies( soup, data )
    add_listOfParameters(soup, parameters)
    add_listOfReactions( soup, data )
    add_listOfObjectives( soup, None )
    add_listOfGenes( soup, None )

    with open(output_file, 'w', encoding="UTF8") if output_file != '-' else sys.stdout as fp:
        fp.write(soup.prettify(formatter="minimal"))

if __name__ == '__main__':
    main()
