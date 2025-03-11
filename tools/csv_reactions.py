"""
Retrieve reaction information from kegg based on a list of reactions

The program takes as input a csv file containing reaction identifiers and
interrogates the
kegg database to recover the name equation and possibly ec number(s) for the
reaction.
The program then recovers, for the different compounds that appear in the
reactions, the name charge and chemical formula from the kegg database.
This information is the written as a csv file suitable for making into a
metabolic model, by the other tools.
"""

import sys
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
        response = requests.get(url+c_id, timeout=5)
        response.raise_for_status()
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

def main():
    """
    Main routine does the work.
    """
    reactions = {}
    compounds = {}
    all_compounds = []
    all_enzymes = []
    for in_line in sys.stdin:
        if len(in_line)<=1 or in_line[0] == '#':
            continue
        for r_id in in_line.split():
            if r_id[0] != 'R':
                print(f'Input format error - {r_id} unexpected.')
            response = requests.get(URL+r_id, timeout=5)
            response.raise_for_status()
            lines = response.text.strip().split("\n")
            r_name = ""
            r_substrates = []
            r_products = []
            r_enzymes = []
            for line in lines:
                if first_word(line) == "NAME":
                    r_name = second_word(line)
                if first_word(line) == "EQUATION":
                    r_substrates, r_products = parse_equation(line)
                if first_word(line) == "ENZYME":
                    r_enzymes = line.split()[1:]

            reactions[r_id] = [r_id,r_name,r_substrates,r_products,r_enzymes]
            all_compounds.extend(r_substrates.split())
            all_compounds.extend(r_products.split())
            all_enzymes.extend(r_enzymes)

    all_compounds = set(all_compounds)
    all_compounds = {item for item in all_compounds if item.startswith('C')}
    all_enzymes = set(all_enzymes)

    compounds = fetch_compounds(URL, all_compounds)

    for _ , compound in sorted(compounds.items()):
        print(f'Mi_{compound[0]};{compound[1]};1;0;{compound[2]}')

    for _ , rxn in sorted(reactions.items()):
#            Should we write all enzymes or the first one?
#        print(f'Ri_{rxn[0]};{rxn[1]};{rename(rxn[2],'Mi_')};{rename(
#            rxn[3],'Mi_')};;{' '.join(rxn[4])}')
#            Optional items should only appear if present
        print(f'Ri_{rxn[0]};{rxn[1]};{rename(rxn[2],'Mi_')};{rename(
            rxn[3],'Mi_')};;{rxn[4][0]}')

if __name__ == '__main__':
    main()
