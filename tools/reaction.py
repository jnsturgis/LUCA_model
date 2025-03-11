"""
This file describes a reaction class
"""

from utils import add_annotation

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
    def from_text(cls, csv_line, sep ):
        """
        Create a reaction based on a csv line containing the information.
        (note the following copied to/from network.py)

        The lines for reactions contain the following information:
	       1. an ID (Rc_XXXXXX where c is compartment and XXXXXX is possibly the
	          kegg id of the equivalent reaction).
	       2. a name for the reaction.
	       3. a set of words representing the reaction substrates possibly
	          preceded by a numerical stoichiometry.
	       4. a set of words representing the reaction products possibly
	          preceded by a numerical stoichiometry.
	       5. a set of words representing any reaction modulators (cofactors and
	          regulators)
	       6. potentially a pair of numbers representing the DG°'m standard free
	          energy change at 1mM standard state in aqueous solution at pH7.0
              (and pMg2+ of 50mM and ionic strength of about 400mM) and the
              precision of this estimate. (in kJ/mol) [nn.nn n.nn]
	       7. potentially a set of EC numbers identifying the enzyme(s)
              responsable for the reaction [x.y.z.zz]
           8. potentially a pair of numbers for the maximum fluxes in forward
              and reverse directions. [nn.nn -nn.nn]
           9. potentially a flux from the last recorded fba [nn.nn]

           N.B. The format of optional fields allows determination of the
           interpretation split() gives 1 (ec,fba) or 2 (dG,flux)
        """

        items = csv_line.split( sep )
        items[2] = _rxn_parse(items[2])
        items[3] = _rxn_parse(items[3])
        kwargs = {}
        if len(items) > 5:
            for item in items[5:]:
                match _info_type(item):
                    case "dG":
                        parts = item.split()
                        kwargs["dG°"] = parts[0]
                        kwargs["dG-uncertainty"] = parts[1]
                    case "EC":
                        kwargs["ec-code"] = item
                    case "FL":
                        parts = item.split()
                        kwargs["max_flux"] = parts[0]
                        kwargs["min_flux"] = parts[1]
                    case "flux":
                        kwargs["flux"] = item
                    case _:
                        pass
        return cls(items[:5], **kwargs)

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

    def __str__(self)-> str:
        """
        print the reaction
        """
        return ( f'** Reaction id:{self.r_id}\n'
                 f'name          :{self.name}\n'
                 f'substrates    :{self.subst}\n'
                 f'products      :{self.prod}\n'
                 f'modifiers     :{self.modif}\n'
                 f'annotations   :{self.dict}\n' )

    def as_sbml(self, soup):
        """
        Return a parsed sbml object for the reaction
        """
        new_reaction = soup.new_tag('reaction',
            attrs={
                "fast":"false",                          # Required
                "reversible":"true",                     # Required
                "id": self.r_id,                         # Required
                "metaid":f'meta_{self.r_id}',            # Optional
                "name":self.name,                        # Optional
                "fbc:lowerFluxBound":self.dict.get("min_flux",
                    "cobra_default_lb"),
                "fbc:upperFluxBound":self.dict.get("max_flux",
                    "cobra_default_ub")
            } )

        for key, value in self.dict.items():
            match key:
                case 'dG°':
                    add_annotation( soup, new_reaction, key, value )
                case 'dG-uncertainty':
                    add_annotation( soup, new_reaction, key, value )
                case 'ec-code':
                    add_annotation( soup, new_reaction, key, value )
                case 'flux':
                    add_annotation( soup, new_reaction, 'fba_flux', value )
                case _:
                    pass

        if len(self.subst):
            new_tag=soup.new_tag('listOfReactants')
            for item in self.subst:
                new_tag.append(soup.new_tag('speciesReference', attrs={
                    "constant":"true",
                    "species":item[0],
                    "stoichiometry":item[1]
                }))
            new_reaction.append(new_tag)

        if len(self.prod):
            new_tag=soup.new_tag('listOfProducts')
            for item in self.prod:
                new_tag.append(soup.new_tag('speciesReference', attrs={
                    "constant":"true",
                    "species":item[0],
                    "stoichiometry":item[1]
                }))
            new_reaction.append(new_tag)

        if len(self.modif):
            new_tag=soup.new_tag('listOfModifiers')
            for species in self.modif:
                new_tag.append(soup.new_tag('modifierSpeciesReference', attrs={
                    "species":species,
                }))
            new_reaction.append(new_tag)

        return new_reaction

def _rxn_parse( mystring ):
    """
    Parse the reagent/product list to a list of lists
    [[reagent, stochiometry],[]...]
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

def _info_type( my_string ):
    """
    Determine the type of information in the csv field from the internal
    structure and return the appropriate code:
    "dG"   free energy and its uncertainty
    "EC"   ec-code
    "FL"   flux limits
    "flux" a flux measurement.
    """
    parts = my_string.split()
    if len(parts) == 0:
        return None
    if len(parts) == 1:
        if parts[0].count('.') > 1:
            return "EC"
        return "flux"
    # else
    if parts[1][0] == '-':
        return "FL"
    return "dG"

def unit_test():
    """
    Test the functions and structures in this file.
    """
