"""
This program reads in an sbml file and generates a series of tables from it,
these include a table of reactions, of compounds and of genes. In each of
these tables the different desired information for the model is collected.
"""

from bs4 import BeautifulSoup

			# These lines control the program functions
PRINT_TABLE = True 	# print tables after building them
SAVE_SBML   = False	# print the xml tree after modifications
VERBOSE     = False	# warnings about formulae and meta_id's
EMPTY = ""

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
        if VERBOSE:
            print( f'Warning no metaid for {specid}.' )
    try:
        formul = species['fbc:chemicalFormula']
    except KeyError:
        # pylint: disable="invalid-name"
        formul = EMPTY
        if VERBOSE :
            print( f'Warning no formula for {specid}.' )
    charge = species['fbc:charge']
    compar = species['compartment']
    info   = []
    for annotation in species.find_all('rdf:li'):
        el  = annotation['rdf:resource'].split('/')		# Fix kegg-compound annotations
        if el[3] == 'kegg' :
            if specid[:5] == 'M_i_C':
                el[3] = 'kegg-compound'
                el[4] = specid[4:]
                annotation['rdf:resource'] = f'{el[0]}/{el[1]}/{el[2]}/{el[3]}/{el[4]}'
        el  = annotation['rdf:resource'].split('/')
        info.append( (el[3],el[4]) )
    if PRINT_TABLE:
    	print( specid, name, formul, charge, compar, info, sep='\t')

for gene in model.find_all('gene'):
    pass

if SAVE_SBML:
    print(model.prettify(formatter="minimal"))
