"""
This file describes a reaction class
"""

class Reaction:
    """
    A class to describe reactions in a metabolic network.
    """
    def __init__(self, items, **kwargs):
        self.r_id  = items[0]          # An id string
        self.name  = items[1]          # A string
        self.subst = items[2]          # A list of lists [compound id, stochiometry]
        self.prod  = items[3]          # A list of lists [compound id, stochiometry]
        self.modif = items[4].split()  # A list of compound_id (split on whitespace)
        self.dict  = {}
        for k, v in kwargs.items():
            self.dict[k] = v

    @classmethod
    def from_text(cls, csv_line):
        """
        Create a raction based on a csv line containing the information.
        """
        items = csv_line.split(',')
        items[2] = rxn_parse(items[2])
        items[3] = rxn_parse(items[3])
        return cls(items)

    def substrates(self):
        """
        Return a list of reaction substrates
        """
        return [r[0] for r in self.subst]

    def products(self):
        """
        Return a list of reaction products
        """
        return [r[0] for r in self.prod]

    def modifiers(self):
        """
        Return a list of the reaction modifiers
        """
        return self.modif


def rxn_parse( mystring ):
    """
    Parse the reagent/product list to a list of tuples.
    """
    result = []
    words = mystring.split()
    stochiometry = 1.0
    for word in words:
        if word[0]=='M':               # Word is a compound_id
            result.append([word, stochiometry])
            stochiometry = 1
        elif float(word) > 0.0:        # Word is a number
            stochiometry = float(word)
        else:
            raise ValueError("Reagent list:'{mystring}' error.")
    return result
