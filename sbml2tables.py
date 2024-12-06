"""
This program reads in an sbml file and generates a series of tables from it,
these include a table of reactions, of compounds and of genes. In each of
these tables the different desired information for the model is collected.
"""

from bs4 import BeautifulSoup

			# These lines control the program functions
PRINT_TABLE = True 	# print tables after building them
SAVE_SBML = False	# print the xml tree after modifications

with open("LUCA.sbml", encoding="UTF8") as fp:
    model = BeautifulSoup( fp, 'xml' )

for reaction in model.find_all('reaction'):
    pass

for compound in model.find_all('compound'):
    pass

for gene in model.find_all('gene'):
    pass

if SAVE_SBML:
    print(model.prettify(formatter="minimal"))
