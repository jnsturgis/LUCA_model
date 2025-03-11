"""
Perform set logic operations on a pair of files.
"""

import sys
import argparse

OPERATIONS = ['UNION','INTERSECTION','DIFFERENCE','SYM_DIFFERENCE']

def get_words(filename):
    """
    Read in a set of words from a file and return the set.

    If the passed filename is '-' then stdin is used.
    """
    with (sys.stdin if filename == "-" else open(filename, "r", encoding="utf-8")) as f:
        return set(f.read().split())  # Read, split into words, and store in a set

def main():
    """
    The main routine handles the logic.
    """
    # Handle the command line using argparse.
    parser = argparse.ArgumentParser(prog="kegg_logic",
        description="Logic operations on sets of words")
    parser.add_argument("operation", type=str,
        help = f'Known operations are {OPERATIONS}')
    parser.add_argument("infile1", type=str, help="Filename or use '-' for stdin")
    parser.add_argument("infile2", type=str, help="Filename or use '-' for stdin")
    parser.add_argument("outfile", type=str, help="Filename or use '-' for stdout")
    args = parser.parse_args()
    if args.operation not in OPERATIONS:
        parser.print_help()
        sys.exit()

    # Read words from both files
    file1_words = get_words(args.infile1)
    file2_words = get_words(args.infile2)
    # Combine the two sets
    match args.operation:
        case "UNION":
            result_words = file1_words.union(file2_words)
        case "INTERSECTION":
            result_words = file1_words.intersection(file2_words)
        case "DIFFERENCE":
            result_words = file1_words.difference(file2_words)
        case "SYM_DIFFERENCE":
            result_words = file1_words.symmetric_difference(file2_words)
        case _:
            # Should be unreachable
            result_words = set()

    with (sys.stdout if args.outfile == "-" else open(args.outfile,
        "w", encoding='utf-8')) as f:
        f.write("\n".join(sorted(result_words)))

# Print results


if __name__ == '__main__':
    main()
