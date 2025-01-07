"""
Rename a sbml model (modify the model id and meta_id fields.)

The model can either be from a command line argument or stdin if none is
specified on the command line.
"""

# Allow todo comments and dont insist on capitals if initialized with a default.
# pylint: disable='fixme' 'invalid-name' 'bare-except'

import sys
from bs4 import BeautifulSoup

def usage():
    """
    Print a brief usage line, as the program was incorrectly called.
    """
    print(f'Usage: {sys.argv[0]} [-h] -n "name" [-o destination][sbml_file]' )
    print( '-h             print this usage information.')
    print( '-n name        the new name for the model.')
    print( '-o destination output resulting sbml file to destination.')
    print( 'sbml_file      the file to read if none is given use stdin.')
    sys.exit()

def command_line( args ):
    """
    Handle the command line options returning the report flag and filenames.
    """
    # pylint: disable='too-many-branches'
    new_name = None
    infile = None
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
                elif item[j] == 'n':
                    if i+1 == len(args):
                        usage()
                    new_name = args[i+1]
                    skip = True
                else:
                    usage()
        else:
            if infile:
                usage()
            infile = item

    if not new_name:
        usage()
    return new_name, infile, outfile

def rename_model(model, name):
    """
    Replace the id with the new name and the metaid as well
    """
    my_model = model.find('model')
    my_model['id'] = name
    my_model['metaid'] = f'meta_{name}'

def main():
    """
    Main routine.
    1. Handle the command line
    2. Read the sbml file.
    3. Change the id.
    4. Write the result.
    """
    new_name, source, dest = command_line(sys.argv)
    try:
        with open(source, 'r', encoding="UTF8") if source else sys.stdin as fp:
            model = BeautifulSoup( fp, 'xml' )
    except:
        print( f"Failed to read sbml model from '{source}'.")
        sys.exit()

    rename_model( model, new_name)

    if dest:
        with open(dest, 'w', encoding="UTF8") as fp:
            fp.write(model.prettify(formatter="minimal"))
    else:
        print(model.prettify(formatter="minimal"))

if __name__ == '__main__':
    main()
