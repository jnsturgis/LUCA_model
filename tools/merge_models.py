"""
Program to merge different sbml files to form a larger model

The program takes a list of .sbml files and merges them together taking species
and reactions and joining the listOfReactions and the listOfSpecies from the
multiple models. The other parts of the sbml file are taken from the first
sbml file in the list.
If no filenames are given on the command line then stdin will be used.
"""

import sys
from bs4 import BeautifulSoup

def usage():
    """
    Print a brief usage line, as the program was incorrectly called.
    """
    print(f'Usage: {sys.argv[0]} [-h]  [-o destination] [sbml_files...]' )
    sys.exit()

def add_model(old_model, new_part):
    pass

def main():
    pass

if __name__ == '__main__':
    main()
