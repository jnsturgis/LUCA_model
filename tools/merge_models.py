"""
Program to merge different sbml files to form a larger model

The program takes a list of .sbml files and merges them together taking species
and reactions and joining the listOfReactions and the listOfSpecies from the
multiple models. The other parts of the sbml file are taken from the first
sbml file in the list.
If no filenames are given on the command line then stdin will be used.
"""

# Allow todo comments and dont insist on capitals if initialized with a default.
# pylint: disable='fixme' 'invalid-name' 'bare-except'

import sys
from bs4 import BeautifulSoup

def usage():
    """
    Print a brief usage line, as the program was incorrectly called.
    """
    print(f'Usage: {sys.argv[0]} [-h] [-o destination] sbml_files...' )
    print( '-h             print this usage information.')
    print( '-o destination output resulting sbml file to destination,')
    print( '               if none given output is to stdout.')
    print( 'sbml_files..   A list of files to read if "-" is given as')
    print( '               a file name stdin will be used.')
    sys.exit()

def command_line( args ):
    """
    Handle the command line options returning the report flag and filenames.
    """
    infiles = []
    outfile = None
    skip = False

    for i in range(1,len(args)):
        if skip:
            skip = False
            continue
        item = args[i]
        if item[0] == '-':
            for j in range(1,len(item)):
                if item[j] == 'h':
                    usage()
                elif item[j] == 'o':
                    if i+1 == len(args):
                        usage()
                    outfile = args[i+1]
                    skip = True
                else:
                    usage()
        else:
            infiles.append(item)

    if not infiles:
        usage()
    return infiles, outfile

def insert_new_member( my_list, member, list_of ):
    """
    Inset the member into the list is it is not already present (append)
    """
    # pylint: disable='undefined-loop-variable'
    for item in my_list.find_all(list_of):
        if item['id'] == member['id']:
            return
    item.insert_after(member)


def add_model(old_model, new_part):
    """
    Merge the new_part into the old_model by adding new species and reactions.
    """
    first_list = old_model.find('listOfSpecies')
    second_list = new_part.find('listOfSpecies')
    for species in second_list.find_all('species'):
        insert_new_member(first_list, species, 'species')

    first_list = old_model.find('listOfReactions')
    second_list = new_part.find('listOfReactions')
    for reaction in second_list.find_all('reaction'):
        insert_new_member(first_list, reaction, 'reaction')

    return old_model

def main():
    """
    Main routine
    """
    infiles, dest = command_line(sys.argv)

    source = infiles.pop(0)
    try:
        with open(source, 'r', encoding="UTF8") if source != '-' else sys.stdin as fp:
            old_model = BeautifulSoup( fp, 'xml' )
    except:
        print( f"Failed to read sbml model from '{source}'.")
        sys.exit()

    while infiles:
        source = infiles.pop(0)
        try:
            with open(source, 'r', encoding="UTF8") if source != '-' else sys.stdin as fp:
                new_part = BeautifulSoup( fp, 'xml' )
        except:
            print( f"Failed to read sbml model from '{source}'.")
            sys.exit()

        old_model = add_model( old_model, new_part )

    if dest:
        with open(dest, 'w', encoding="UTF8") as fp:
            fp.write(old_model.prettify(formatter="minimal"))
    else:
        print(old_model.prettify(formatter="minimal"))

if __name__ == '__main__':
    main()
