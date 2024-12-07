"""
This program reads in an sbml file and generates a series of tables from it,
these include a table of reactions, of compounds and of genes. In each of
these tables the different desired information for the model is collected.
"""

from bs4 import BeautifulSoup

			# These lines control the program functions
PRINT_TABLE = True 	# print tables after building them
SAVE_SBML = False	# print the xml tree after modifications
EMPTY = ""

def annot2key_value( text: str ) -> tuple:
    """Parse the annotation into the key and value as a tuple"""
    elements = text.split('/')
    return (elements[3],elements[4])

with open("LUCA.sbml", encoding="UTF8") as fp:
    model = BeautifulSoup( fp, 'xml' )

for reaction in model.find_all('reaction'):
    pass

for species in model.find_all('species'):
    specid = species['id']
    name   = species['name']
    try:
        metaid = species['metaid']
    except KeyError:
        metaid = 'Meta_'+specid
        if not PRINT_TABLE:
            print( f'Warning no metaid for {specid}.' )
    try:
        formul = species['fbc:chemicalFormula']
    except KeyError:
        # pylint: disable="invalid-name"
        formul = EMPTY
        if not PRINT_TABLE :
            print( f'Warning no formula for {specid}.' )
    charge = species['fbc:charge']
    compar = species['compartment']
    info   = []
    for annotation in species.find_all('rdf:li'):
        info.append( annot2key_value(annotation['rdf:resource'] ))
    infodict = dict(info)
    if 'kegg' in infodict :
        if specid[:5] == 'M_i_C':
            infodict['kegg-compound'] = specid[4:]
        del infodict['kegg']
    if PRINT_TABLE:
        print( specid + '\t' + name + '\t' + formul + '\t', end='' )
        print( charge + '\t' + compar + '\t' , infodict, sep="" )

for gene in model.find_all('gene'):
    pass

if SAVE_SBML:
    print(model.prettify(formatter="minimal"))
