"""
This program reads in an sbml file and generates a series of tables from it,
these include a table of reactions, of compounds and of genes. In each of
these tables the different desired information for the model is collected.

Note not everything survives a passage through cnapy.
Notably reversibility (perhaps it is due to fluxLowerBound?)
"""

# pylint: disable='fixme'
from bs4 import BeautifulSoup

			# These lines control the program functions
PRINT_STABLE = False 	# print information in tables after building them
PRINT_RTABLE = False
PRINT_GTABLE = False
SAVE_SBML    = False	# save the xml tree after modifications
VERBOSE      = True	# warnings about formulae and meta_id's

EMPTY        = ""
MINIMUM      = -15.0

with open("LUCA.sbml", encoding="UTF8") as fp:
    model = BeautifulSoup( fp, 'xml' )

for reaction in model.find_all('reaction'):
    reacid = reaction['id']
    name   = reaction['name']
    revers = reaction['reversible']
    substrate  = []
    sublist = reaction.find('listOfReactants')
    if sublist is not None :
        for species in sublist.find_all('speciesReference'):
            substrate.append((species['species'],species['stoichiometry']))
    product    = []
    prodlist = reaction.find('listOfProducts')
    if prodlist is not None :
        for species in prodlist.find_all('speciesReference'):
            product.append((species['species'],species['stoichiometry']))
    modifier   = []
    modlist = reaction.find('listOfModifiers')
    if modlist is not None :
        for species in modlist.find_all('modifierSpeciesReference'):
            product.append(species['species'])
    info       = []
    for annotation in reaction.find_all('rdf:li'):
        el  = annotation['rdf:resource'].split('/')
        if el[3] == 'ec #':                             # fix ec-code annotations
            if el[4][:5] == 'EC #:' :
                el[4] = el[4][5:]
            el[3] = 'ec-code'
            annotation['rdf:resource'] = f'{el[0]}/{el[1]}/{el[2]}/{el[3]}/{el[4]}'
        if el[3] == 'dG0' :                             # Fix reversibility based on dG0
            if float( el[4] ) > MINIMUM :               # TODO fbc:lowerFluxBound
                # pylint: disable="invalid-name"
                revers = 'true'
                reaction['reversible'] = revers
                try:
                    reaction['fbc:lowerFluxBound'] = 'cobra_default_lb'
                except KeyError:
                    pass
        el  = annotation['rdf:resource'].split('/')
        info.append( (el[3],el[4]) )

    infodict = dict(info)
    if 'dG0' not in infodict:
        if VERBOSE:
            print( f'Warning {reacid} has no dG0.')
    if 'kegg-reaction' not in infodict:                 # add kegg-reaction annotation
        if name[:5] == 'KEGG ':
            # pylint: disable='undefined-loop-variable'
            el[3] = 'kegg-reaction'
            el[4] = name[5:]
            infodict[el[3]] = el[4]                    # update dictionary
            info.append((el[3],el[4]))                 # update list
            new_tag = model.new_tag('rdf:li')          # update sbml tree
            new_tag['rdf:resource'] = f"https://identifiers.org/{el[3]}/{el[4]}"
            annotation.insert_after(new_tag)
        else :
            if VERBOSE:
                print(f'Warning {reacid} has no kegg identifier.')

    if PRINT_RTABLE:
        print( reacid, name, revers, substrate, product, modifier, info, sep='\t')

for species in model.find_all('species'):		# TODO ensure species have initialConcentration
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
    try:
        conc = float(species['initialConcentration'])
    except KeyError:
        # pylint: disable="invalid-name"
        species['initialConcentration'] = '1.0'
        conc = 1.0
    info   = []
    for annotation in species.find_all('rdf:li'):
        el  = annotation['rdf:resource'].split('/')	# Fix kegg-compound annotations
        if el[3] == 'kegg' :
            if specid[:5] == 'M_i_C':
                el[3] = 'kegg-compound'
                el[4] = specid[4:]
                annotation['rdf:resource'] = f'{el[0]}/{el[1]}/{el[2]}/{el[3]}/{el[4]}'
        el  = annotation['rdf:resource'].split('/')
        info.append( (el[3],el[4]) )
    if PRINT_STABLE:
        print( specid, name, formul, charge, compar, info, sep='\t')

for gene in model.find_all('gene'):
    pass

if SAVE_SBML:
    print(model.prettify(formatter="minimal"))
