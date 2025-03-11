#!/usr/bin/env python3
"""
Retrieve enzyme information from uniprot based on a csv file containing network
information.

Using the ec-code as a starting point collect proteins in uniprot that
correspond, with a preference for those from model organisms and not too many
perhaps 20 maximum.

Then for the protein entries that have been found check if the enzyme has a
cofactor requirement. If so find the corresponding compound if it is already
in the model using the kegg identifier to find it.
"""

import argparse
import requests
import network

def uniprot_find_proteins( ec_code ):
    """
    Find a list of proteins that correspond to the passed code, return a maximum
    of 20 protein codes preferably from model organisms in order of preference.
    """
    url = ( f'https://rest.uniprot.org/uniprotkb/'
            f'search?query=reviewed:true+AND+EC:{ec_code}'
            f'+AND+taxonomy_id:2+OR+taxonomy_id:2157&format=tsv')
    with requests.get(url, timeout=10) as request:
        request.raise_for_status()
        lines=request.text.strip().split("\n")[1:]
    codes = []
    for line in lines:
        codes.append(line.split('\t')[0])
        if len(codes) >= 10:
            return codes
    if len(codes)== 0:
        return None
    return codes

def uniprot_find_cofactors( pr_code ):
    """
    Check if enzyme has a cofactor requirement, and if so add the cofactors to
    the list of cofactors (initially as CHEBI_id's but would prefer KEGG).
    """
    url = ( f'https://rest.uniprot.org/uniprotkb/'
            f"search?query={pr_code}&format=txt")
    print(url)
    codes = []
    with requests.get(url, timeout=10) as request:
        request.raise_for_status()
        lines=request.text.strip().split("\n")
    for counter, line in enumerate(lines):
        if "COFACTOR" in line:
            while '-!-' not in lines[counter + 1]:
                counter += 1
                print(lines[counter])
                if 'Name' in lines[counter]:
                    codes.append(f'ChEBI:{lines[counter].split(';')[1].split(':')[2]}')
    if len(codes) == 0:
        return None
    return set(codes)

def main():
    """
    The main routine.
    """
    parser = argparse.ArgumentParser(
        prog = 'csv_get_enzyme_info.py',
        description = 'Fetch enzyme info for a metabolic network in a csv file.'
    )
    parser.add_argument('-u','--uniprot',
        action='store_true',
        help='Fill enzyme info from uniprot database.')
    parser.add_argument('-o','--outfile', default='-',
        help= "Destination file, '-' means stdout.")
    parser.add_argument('infile', nargs='?', default=['-'],
        help='File to translate where "-" is stdin (,defaults to stdin).')
    args = parser.parse_args()

    my_network = network.Network.from_csv(args.infile)
    if args.uniprot :
        for ec_code, enzyme in my_network.enzymes.items():
            enzyme.add_proteins(uniprot_find_proteins(ec_code))
            if len(enzyme.proteins) > 0:
                pr_code = enzyme.proteins[0]
                enzyme.add_cofactors(uniprot_find_cofactors(pr_code))

    for _, enzyme in sorted(my_network.enzymes.items()):
        print(enzyme)

if __name__ == '__main__':
    main()
