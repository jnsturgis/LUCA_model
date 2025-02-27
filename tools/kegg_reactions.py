"""
Retrieve reaction information from kegg based on a list of reactions

The program takes as input a file containing reaction a list of identifiers and
interrogates the kegg database to recover the name equation and possibly ec
number(s) for the reaction.

The program then recovers, for the different compounds that appear in the
reactions, the name charge and chemical formula from the kegg database.

This information is the written as a csv file suitable for making into a
metabolic model, by the other tools.
"""

import sys
import argparse
import requests

# Constants definitions
URL = 'https://rest.kegg.jp/get/'           # URL of kegg database api.

def first_word(line: str) -> str:
    """
    Get the first white space separated word of a line
    """
    return line.split(None)[0]

def second_word(line:str) -> str:
    """
    Get the second word of a line, possibly dividing on ;
    """
    words = line.split(None,1)[1]
    return words.split(';')[0]

def parse_equation(line: str):
    """
    Parse an equation line
    """
    equation = line.split(None,1)[1]
    parts = equation.split('=')
    substrates = parts[0][:-2]
    substrates = " ".join([word if word != '+' else '' for word in substrates.split()])
    products   = parts[1][2:]
    products   = " ".join([word if word != '+' else '' for word in products.split()])
    return substrates, products

def rename( line, prefix ):
    """
    Add prefix to compounds.
    """
    return " ".join([prefix + word if word.startswith('C') else word for word in line.split()])

def fetch_compounds( url: str, all_compounds: set) -> dict:
    """
    Fetch the compounds information for database
    """
    compounds = {}
    for c_id in all_compounds:
        if c_id[0] != 'C' or len(c_id) != 6:
            print(f'Parsing error - {c_id} not a valid kegg compound id. **IGNORED**')
        else:
            response = requests.get(url+c_id, timeout=5)
            try:
                response.raise_for_status()
            except requests.exceptions.HTTPError:
                print(f'Error retrieving data for "{c_id}" **IGNORED**')
                continue

            lines = response.text.strip().split("\n")
            c_name = ""
            c_formula = ""
            for line in lines:
                if first_word(line) == "NAME":
                    c_name = second_word(line)
                if first_word(line) == "FORMULA":
                    c_formula = second_word(line)
            compounds[c_id] = [c_id, c_name, c_formula]
    return compounds

def fetch_reactions( url: str, reaction_set: set ) -> dict:
    """
    Fetch the reaction information from the database
    """
    reactions = {}

    for r_id in reaction_set:
        if r_id[0] != 'R' or len(r_id) != 6:
            print(f'Input format error - {r_id} not a valid kegg reaction id. **IGNORED**')
        else:
            response = requests.get(url+r_id, timeout=5)
            try:
                response.raise_for_status()
            except requests.exceptions.HTTPError:
                print(f'Error retrieving data for "{r_id}" **IGNORED**')
                continue

            lines = response.text.strip().split("\n")
            r_name = ""
            r_substrates = []
            r_products = []
            for line in lines:
                if first_word(line) == "NAME":
                    r_name = second_word(line)
                if first_word(line) == "EQUATION":
                    r_substrates, r_products = parse_equation(line)
            reactions[r_id] = [r_id,r_name,r_substrates,r_products]
    return reactions

def main():
    """
    Main routine does the work.
    """
    # Handle the command line using argparse.
    parser = argparse.ArgumentParser(prog="kegg_reactions",
        description="Interrogate kegg to build metabolic network csv file")
    parser.add_argument("infile", type=str, help="Filename or use '-' for stdin")
    parser.add_argument("outfile", type=str, help="Filename or use '-' for stdout")
    args = parser.parse_args()

    # Read words from input
    with (sys.stdin if args.infile == "-" else open(args.infile,
        "r", encoding="utf-8")) as f:
        reaction_set = set(f.read().split())

    reactions = fetch_reactions(URL, reaction_set )
    compound_set = []
    for _, rxn in reactions.items():
        compound_set.extend(rxn[2].split())
        compound_set.extend(rxn[3].split())
    compounds = fetch_compounds(URL, set(compound_set))

    # Output the csv file
    with (sys.stdout if args.outfile == "-" else open(args.outfile,
        "w", encoding='utf-8')) as f:
        for _ , compound in sorted(compounds.items()):
            f.write(f'Mi_{compound[0]};{compound[1]};1;0;{compound[2]}\n')
        for _ , rxn in sorted(reactions.items()):
            f.write(f'Ri_{rxn[0]};{rxn[1]};{rename(rxn[2],'Mi_')};{rename(rxn[3],'Mi_')};\n')

if __name__ == '__main__':
    main()
