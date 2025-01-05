"""
Rename a sbml model (modify the model id and meta_id fields.)

The model can either be from a command line argument or stdin if none is
specified on the command line.
"""
import sys
from bs4 import BeautifulSoup

def usage():
    """
    Print a brief usage line, as the program was incorrectly called.
    """
    print(f'Usage: {sys.argv[0]} [-h] -n "name" [-o destination][sbml_file]' )
    sys.exit()

def rename_model(model, name):
    pass

def main():
    pass

if __name__ == '__main__':
    main()
