"""
This file describes an enzyme class
"""

class Enzyme:
    """
    A class to describe enzymes in a metabolic network.
    """
    def __init__(self, e_id ):
        self.e_id      = e_id       # The ec number of the enzyme
        self.reactions = []         # List of reaction id's catalyzed
        self.cofactors = []         # List of bound cofactors
        self.proteins  = []         # List of proteins uniprot references
        self.gene_rule = []         # To be determined.

    def __str__(self)-> str:
        """
        print the enzyme info
        """
        return ( f'** ec-code ** :{self.e_id}\n'
                 f'reactions     :{self.reactions}\n'
                 f'cofactors     :{self.cofactors}\n'
                 f'proteins      :{self.proteins}\n'
                 f'gene rule     :{self.gene_rule}\n')

    def add_reaction( self, r_id ):
        """
        Add the id of a reaction that is catalysed by this enzyme
        """
        self.reactions.append(r_id)

    def add_proteins( self, prot_list ):
        """
        Add the list of proteins in prot_list to the proteins, should be None or
        a list of uniprot identifiers.
        """
        if prot_list:
            self.proteins.extend(prot_list)

    def add_cofactors( self, cof_list ):
        """
        Add the list of cofactors in cof_list to the cofactors, should be None
        or a list of chebi/kegg identifiers.
        """
        if cof_list:
            self.cofactors.extend(cof_list)
