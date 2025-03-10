"""
This file describes a metabolic network class and ancillary functions.

The objective of this work is to build whole organism (or meta-organism)
metabolic models and study how these can evolve and adapt to different
conditions.

In this file and the associated class we will try to adhere to the following
scheme.

Instance methods of the class modify the data.
Class methods operate on the class or are alternatice constructors
Standalone functions operate with multiple instances or are very generic.
"""

# pylint: disable=W0511
# TODO Needs unit testing
import sys
import numpy as np

from bs4 import BeautifulSoup

import networkx as nx
from networkx.algorithms import bipartite

import reaction as rxn
import compound as cpd

SBML_INITIAL = """<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" xmlns:fbc="http://www.sbml.org/sbml/level3/version1/fbc/version2" sboTerm="SBO:0000624" level="3" version="1" fbc:required="false">
<model metaid="meta_model_name" id="model_name" fbc:strict="true">
</model>
</sbml>
"""
class Network:
    """
    A class to describe metabolic networks.

    A metabolic network is made of:
    1. A set of compounds as a dictionary.
    2. A set of reactions in a dictionary.
    3. A set of enzymes in a dictionary.
    """
    def __init__(self, *argv):
        match len(argv):
            case 3:
                self.compounds = argv[0]
                self.reactions = argv[1]
                self.enzymes   = argv[2]
            case 0:
                self.compounds = {}
                self.reactions = {}
                self.enzymes   = {}
            case _:
                raise ValueError("Need 0 or 3 dictionaries to make a metabolic network.")

    @classmethod
    def from_sbml(cls, sbml_name):
        """
        Use an sbml file to setup a metabolic network.
        """

    @classmethod
    def from_csv(cls, csv_name):
        """
        Use a csv file to setup a metabolic network.

        `.csv` condensed description of a model in csv format with ';' as the
        separator with lines for metabolites (start with letter M) and reactions
		(start with letter R). the '#' character can be used for comments, and long
		lines can be broken with a '\' before the EOL character. Lines for
		metabolites contain:
	       1. an ID (Mc_XXXXXX where c is compartment and XXXXXX is possibly the kegg
	          id of the equivalent molecule).
	       2. a name for the molecule.
	       3. an initial concentration of the molecule in the model.
	       4. the charge of the molecule (coherent with the reactions and formula).
	       5. the formula of the molecule with the abbreviations Pr, Rn and Dn
	          corresponding to generic protein (polypeptide), RNA and DNA molecules.
           6. Potentially more information that is parsed as: Nothing yet.
        The lines for reactions contain the following information:
	       1. an ID (Rc_XXXXXX where c is compartment and XXXXXX is possibly the kegg
	          id of the equivalent reaction).
	       2. a name for the reaction.
	       3. a set of words representing the reaction substrates possibly preceded by
	          a numerical stoichiometry.
	       4. a set of words representing the reaction products possibly preceded by
	          a numerical stoichiometry.
	       5. a set of words representing any reaction modulators (cofactors and
	          regulators)
	       6. potentially a pair of numbers representing the DG°'m standard free energy
	          change at 1mM standard state in aqueous solution at pH7.0 (and pMg2+ of
			  50mM and ionic strength of about 400mM)
	       7. potentially a set of EC numbers identifying the enzyme(s) responsable for
	          the reaction.
           8. potentially a pair of numbers for the maximum fluxes in forward and
	          reverse directions.
	       9. potentially a flux from the last recorded fba.
        """

        compounds = {}
        reactions = {}
        enzymes   = {}
        old_line = ""

        sep = ";"

        with (sys.stdin if csv_name == "-" else open(csv_name, "r", encoding='utf-8')) as f:
            for line in f:
                line = line.rstrip()
                if len(old_line) > 0:
                    line = old_line + line
                if line[-1] != '\\':        # Continuation character
                    parts = line.split(sep,1)
                    match line[0]:
                        case '#':           # Comment
                            continue
                        case 'R':           # Reaction line
                            if parts[0] in reactions:
                                print(f'Duplicate reaction {parts[0]}.')
                            reactions[parts[0]] = rxn.Reaction.from_text(line, sep)
                        case 'M':           # Metabolite line
                            if parts[0] in compounds:
                                print(f'Duplicate compound {parts[0]}.')
                            compounds[parts[0]] = cpd.Compound.from_text(line, sep)
                        case 'E':           # Enzyme line
                            continue
                        case _:             # Oops error
                            raise ValueError(f'Invalid input line \n{line}')
                    old_line = ""
                else:
                    old_line = line[:-1]
        return cls(compounds, reactions, enzymes)

    def add_compound( self, compound_id, compound_data ):
        """
        Add a new compound to a metabolic network.
        """
        self.compounds[compound_id] = compound_data

    def add_reaction(self, reaction_id, reaction_data ):
        """
        Add a new reaction to a metabolic network.
        """
        self.reactions[reaction_id] = reaction_data

    def compartments(self) -> list:
        """
        Get a list of compartments in the metabolic network.
        """
        return ['i','e']

    def write_sbml(self, outfile):
        """
        Write the model as an sbml file with the passed name.
        """
        parameters = ["cobra_default_lb", "-1000", "cobra_default_ub", "1000"]
        soup = BeautifulSoup( SBML_INITIAL ,"xml")

        add_list_unitdefinitions( soup, ["mmol_per_gDW_per_hr"] )
        add_list_compartments( soup, self.compartments() )
        add_list_species( soup, self.compounds )
        add_list_parameters(soup, parameters)
        add_list_reactions( soup, self.reactions )
        add_list_objectives( soup, None )
        add_list_genes( soup, None )

        with open(outfile, 'w',
            encoding="UTF8") if outfile != '-' else sys.stdout as fp:
            fp.write(soup.prettify(formatter="minimal"))

def find_modifiers(network):
    """
    Return the set of all modifiers that feature in the reactions of the network.
    """
    modifiers = []
    for _, reaction in network.reactions.items():
        modifiers.extend(reaction.modifiers())
    return set(modifiers)

def make_graph( network ):
    """
    Make a network graph from the metabolic network information using networkx.

    In the resulting graph, produced from the reactions, the nodes correspond
    to the reactions and compounds featuring as substrates and products of the
    reactions. The graph is undirected and unweighted. The node labels
    correspond to the reaction and compound id's is used to set the bipartite
    attrbute (metabolites = 0, other things = 1) allowing consideration as a
    bipartite graph.

    The initial letter codes for id's are:
        * Mx_ Metabolite in compartment x (x=i or e)
        * R_  Reaction normal type.

    """
    graph = nx.Graph()
    for r_id, reaction in network.reactions.items():
        for c_id in reaction.substrates():
            graph.add_edge(c_id, r_id)
        for c_id in reaction.products():
            graph.add_edge(r_id, c_id)
    # set attribute bipartite for nodes based on first letter of node_id...
    nx.set_node_attributes(graph,
        {node: 0 if node.startswith("M") else 1 for node in graph.nodes}, "bipartite")
    assert nx.is_bipartite(graph)
    return graph

def make_dgraph( network, **kwargs):
    """
    Make a directed network graph from the metabolic network information using networkx.

    In the resulting graph, produced from the reactions, the nodes correspond
    to the reactions and compounds featuring as substrates and products of the
    reactions.
    The graph is directed and weighted based on the standard free energy change
    associated with the reaction in the forward and reverse directions using
    the metropolis scheme.
    TODO different weighting options (flux, dG°'m etc)
    The node labels correspond to the reaction and compound id's is used to set
    the bipartite attrbute (reations = 1, other things = 0) allowing
    consideration as a bipartite graph.
    """
    graph = nx.DiGraph()

    options = {
        'weighting': 'dG0'
    }
    if kwargs:
        options.update({k: v for k, v in kwargs.items() if k in options})

    for r_id, reaction in network.reactions.items():

        match options['weighting']:
    # TODO add weighting algorithm and check others
            case 'dG0':
                forward_weight = 1.0        # This should be based of DG°'m'
                reverse_weight = 0.0
            case 'flux':
                flux = reaction.dict.get('flux', 0.0)
                if flux > 0:
                    forward_weight = flux   # This should be based of DG°'m'
                    reverse_weight = 0.0
                else:
                    forward_weight = 0.0
                    reverse_weight = -flux
            case '_':
                raise ValueError("Unrecognized weighting scheme.")

        for c_id in reaction.substrates():
            graph.add_weighted_edges_from([c_id, r_id, forward_weight])
            graph.add_weighted_edges_from([r_id, c_id, reverse_weight])
        for c_id in reaction.products():
            graph.add_weighted_edges_from([r_id, c_id, forward_weight])
            graph.add_weighted_edges_from([c_id, r_id, reverse_weight])
    # set attribute bipartite for nodes based on first letter of node_id...
    nx.set_node_attributes(graph,
        {node: 0 if node.startswith("M") else 1 for node in graph.nodes}, "bipartite")
    assert nx.is_bipartite(graph)
    return graph

def reaction_adjacency_graph( network ):
    """
    Calculate the reaction adjacency matrix for a metabolic network.

    This is the projection of the metabolic network bipartite graph onto the
    reaction node, with edges weighted by the number of connections.
    """
    # TODO check weighting algorithm
    bi_graph = make_graph( network )
    result = bipartite.weighted_projected_graph(
        bi_graph, nodes = {n for n, d in bi_graph.nodes(data=True) if d["bipartite"] == 1})
    return result

def normalized_flow_graph( network ):
    """
    Calculate the normalized_flow_graph for a metabolic network based on weighted graph

    This is the projection of the weighted directions metabolic network
    bipartite graph onto the reaction node, with edges weighted by the edge
    weightings.
    """
    # TODO check algorithm
    bi_graph = make_dgraph( network )
    result = bipartite.weighted_projected_graph(
        bi_graph, nodes = {n for n, d in bi_graph.nodes(data=True) if d["bipartite"] == 1})
    return result

def mass_flow_graph( network ):
    """
    Calculate the reaction adjacency matrix for a metabolic network.

    This is the projection of the metabolic network bipartite graph onto the
    reaction node, with edges weighted by the number of connections.
    """
    # TODO check algorithm
    bi_graph = make_dgraph( network, weighting='flux')
    result = bipartite.weighted_projected_graph(
        bi_graph, nodes = {n for n, d in bi_graph.nodes(data=True) if d["bipartite"] == 1})
    return result

def add_atoms( dict1, dict2 ):
    """
    Sum by atom type the values in dict1 and dict2
    """
    keys = set(dict1) | set(dict2)  # Union of keys
    result = {k: dict1.get(k, 0) + dict2.get(k, 0) for k in keys}
    return result

def check_rxn_balance(reaction, compounds):
    """
    Check the balance of a reaction for atoms and charges
    """
    charge = 0.0
    atoms = {}
    for substrate in reaction.subst:
        c_id = substrate[0]
        num  = substrate[1]
        compound = compounds[c_id]
        charge += num * float(compound.charge)
        new_atoms = compound.atoms()
        for k, v in new_atoms.items():
            new_atoms[k] = v * num
        atoms = add_atoms(atoms, new_atoms)
    for product in reaction.prod:
        c_id = product[0]
        num  = product[1] * -1
        compound = compounds[c_id]
        charge += num * float(compound.charge)
        new_atoms = compound.atoms()
        for k, v in new_atoms.items():
            new_atoms[k] = v * num
        atoms = add_atoms(atoms, new_atoms)
    if not (charge == 0.0 and all(value == 0 for value in atoms.values())):
        s1 = f'Reaction {reaction.r_id} unbalanced: '
        if charge != 0.0:
            s1 += f'charge {charge}'
        for k, v in atoms.items():
            if v != 0.0:
                s1 += f",{k} {v}"
        print( s1 )
        return False
    return True

def check_connectivity(network):
    """
    Check that the metabolic network is connected and indicate divisions if there are any

    network is a Network structure
    exclude is a list of ids to remove from the graph before checking for sub_graphs.

    Returns the networkx undirected graph of the network based on reactions without
    any unused compounds that might feature as modifiers.
    """
    graph = make_graph(network)
    nodes = list(graph.nodes)
    unrefed = find_modifiers(network).union(nodes)
    for c_id, cmpd in network.compounds.items():
        if c_id not in unrefed:
            print(f"Compound {c_id} '{cmpd.name}'is not used in any reactions.")
    sub_graphs = list(nx.connected_components(graph))
    print(f'Network has {len(sub_graphs)} connected components.')
    if len(sub_graphs) > 1:
        for item in sub_graphs:
            print(item)
    return graph

def check_balance(network):
    """
    Check that the reactions are balanced for charge and ChemicalFormula and warn of problems
    """
    balanced = True
    for _, reaction in network.reactions.items():
        test = check_rxn_balance(reaction, network.compounds)
        balanced = balanced and test
    if balanced:
        print( 'All reactions are chemically balanced.')

def analyse_pathway(network, exclude, *args):
    """
    Perform various analyses on the network! A bit generic.

    1. Find loops
    2. Notify if no loop metabolite is a reaction modulator
    """
    if len(args) == 1:
        graph = args[0]
    else:
        graph = make_graph( network)
    for item in exclude:
        if graph.has_node(item):
            graph.remove_node(item)
    cycles = list(nx.cycle_basis(graph))
    print(f'Network has {len(cycles)} base cycles.')
    modifiers = find_modifiers(network)
    print(f'Modifiers are :{modifiers}')
    for ring in cycles:
        mods = modifiers.intersection(ring)
        if mods:
            print (f"Modifiers {mods} appear in cycle: {ring}")
        else:
            print (f"No modifiers in cycle: {ring}")
    return graph

def stochiometric_matrix( network ):
    """
    Calculate the stochiometric matrix S for a metabolic network

    S(i,j) is the stochiometry of production of compound(i) by reaction (j),
    thus reagents get negative stochiometries.
    Compound indexes can be found from ??? and reaction indexes from ???
    """
    reaction_keys = list(network.reactions.keys())
    compound_keys = list(network.compounds.keys())
    s_matrix = np.zeros((len(compound_keys),len(reaction_keys)))
    # pylint: disable=C0200
    for j in range(0,len(reaction_keys)):
        reaction = network.reactions[reaction_keys[j]]
        for substrate in reaction.subst:
            val  = - substrate[1]
            i    = compound_keys.index(substrate[0])
            s_matrix[i][j] = val
        for substrate in reaction.prod:
            val  = substrate[1]
            i    = compound_keys.index(substrate[0])
            s_matrix[i][j] = val
    # pylint: enable=C0200
    return s_matrix, reaction_keys, compound_keys

def find_compartments( network ):
    """
    Return a set of compartment codes for the network, extracted from compound id's
    """
    compartment_set = set([])
    for compound in list(network.compounds.keys()):
        compartment_set.add(compound[1])
    return compartment_set

def add_annotation( soup, item ):
    """
    Add to 'item' the xml bits to hold annotation key:value pairs
    """
    new_tag = soup.new_tag('annotation')
    item.append(new_tag)
    new_tag.append(soup.new_tag('rdf:RDF', attrs={
        "xmlns:bqbiol":"http://biomodels.net/biology-qualifiers/",
        "xmlns:bqmodel":"http://biomodels.net/model-qualifiers/",
        "xmlns:dcterms":"http://purl.org/dc/terms/",
        "xmlns:rdf":"http://www.w3.org/1999/02/22-rdf-syntax-ns#",
        "xmlns:vCard":"http://www.w3.org/2001/vcard-rdf/3.0#",
        "xmlns:vCard4":"http://www.w3.org/2006/vcard/ns#"}))
    new_tag = new_tag.find('rdf:RDF')
    new_tag.append(soup.new_tag('rdf:Description', attrs={'rdf:about':f'#{item["metaid"]}'}))
    new_tag = new_tag.find('rdf:Description')
    new_tag.append(soup.new_tag('bqbiol:is'))
    new_tag = new_tag.find('bqbiol:is')
    new_tag.append(soup.new_tag('rdf:Bag'))

def add_annotations( soup, item, item_list ):
    """
    add key value pairs from the 'item_list' to the annotation of 'item'
    """
    new_tag = item.find('rdf:Bag')
    for label, value in zip(item_list.split()[::2], item_list.split()[1::2]):
        new_tag.append(soup.new_tag('rdf:li',attrs={
            'rdf:resource': f'https://identifiers.org/{label}/{value}'
        }))

def add_list_unitdefinitions( soup, units_definitions ):
    """
    Add to the model a list of unit definitions as necessary.

    Passed parameters:
        soup                a BeautifulSoup xml parse tree.
        units_definitions   a list of strings.

    """
    # TODO: This needs parsing of units expression (read definition)
    # if necessary add a listOfUnitDefinitions
    my_model = soup.find("model")
    if units_definitions:
        my_list = soup.new_tag("listOfUnitDefinitions")
        my_model.append(my_list)
        for term in units_definitions:
            # if necessary add the unitDefinition to the list
            new_term = soup.new_tag("unitDefinition", attrs={"id":term})
            my_list.append(new_term)
            new_tag = soup.new_tag("listOfUnits")
            new_term.append(new_tag)
            words = term.split('_')
            # Need to think how to do this properly and not just assume the term.
            if not words[0] == 'mmol':
                pass
            new_unit = soup.new_tag("unit", attrs={
                "kind":"mole" ,
                "exponent":"1",
                "scale":"-3",
                "multiplier":"1"})
            new_tag.append(new_unit)
            new_unit = soup.new_tag("unit", attrs={
                "kind":"gram" ,
                "exponent":"-1",
                "scale":"0",
                "multiplier":"1"})
            new_tag.append(new_unit)
            new_unit = soup.new_tag("unit", attrs={
                "kind":"second" ,
                "exponent":"-1",
                "scale":"0",
                "multiplier":"3600"})
            new_tag.append(new_unit)

def add_list_compartments( soup, compartments ):
    """
    Add to the model a list of compartments as necessary.
    """
    model = soup.find("model")
    mylist = model.find("listOfCompartments")
    if not mylist:
        new_tag = soup.new_tag("listOfCompartments",)
        model.append( new_tag )
        mylist = new_tag
    if compartments:
        for compartment in compartments:
            old = mylist.find("compartment", attrs={"id":compartment})
            if not old:
                new_tag = soup.new_tag( "compartment", attrs={"id":compartment}, constant="true")
                mylist.append(new_tag)

def add_list_species( soup, metab_model ):
    """
    Add to the model a list of species as necessary from the DataFrame.
    """
    model = soup.find('model')
    mylist = model.find('listOfSpecies')
    if not mylist:
        new_tag = soup.new_tag('listOfSpecies')
        model.append(new_tag)
        mylist = new_tag
    for key, species in metab_model.items():
        new_species = soup.new_tag('species',
            attrs={
                "boundaryCondition":"false",
                "compartment":species.compartment,
                "constant":"false",
                "fbc:charge":int(species.charge),
                "fbc:chemicalFormula":species.formula,
                "hasOnlySubstanceUnits":"false",
                "id":key,
                "initialConcentration":species.concentration,
                "metaid":f'meta_{key}',
                "name":species.name
            } )
        mylist.append(new_species)
#        if not pd.isna(row.Dbases) and len(row.Dbases) > 0:
#            add_annotation( soup, new_species)
#            add_annotations( soup, new_species, row.Dbases )

def add_list_reactions( soup, matab_model ):
    """
    Add to the model a list of reactions as necessary from the DataFrame.
    """
    model = soup.find('model')
    mylist = model.find('listOfReactions')
    if not mylist:
        new_tag = soup.new_tag('listOfReactions')
        model.append(new_tag)
        mylist = new_tag
    for key, reaction in matab_model.items():
        new_reaction = soup.new_tag('reaction',
            attrs={
                "fast":"false",                          # Required
                "reversible":"true",                     # Required
                "id": key,                               # Required
                "metaid":f'meta_{key}',                  # Optional
                "name":reaction.name,                    # Optional
                "fbc:lowerFluxBound":"cobra_default_lb", # Required fbc:strict
                "fbc:upperFluxBound":"cobra_default_ub"  # Required fbc:strict
            } )
        mylist.append(new_reaction)

        if len(reaction.subst):
            new_tag=soup.new_tag('listOfReactants')
            for item in reaction.subst:
                new_tag.append(soup.new_tag('speciesReference', attrs={
                    "constant":"true",
                    "species":item[0],
                    "stoichiometry":item[1]
                }))
            new_reaction.append(new_tag)

        if len(reaction.prod):
            new_tag=soup.new_tag('listOfProducts')
            for item in reaction.prod:
                new_tag.append(soup.new_tag('speciesReference', attrs={
                    "constant":"true",
                    "species":item[0],
                    "stoichiometry":item[1]
                }))
            new_reaction.append(new_tag)

        if len(reaction.modif):
            new_tag=soup.new_tag('listOfModifiers')
            for species in reaction.modif:
                new_tag.append(soup.new_tag('modifierSpeciesReference', attrs={
                    "species":species,
                }))
            new_reaction.append(new_tag)

#        if not pd.isna(row.GeneRules) and len(row.GeneRules) > 0:
            # TODO: This needs parsing of gene rules expression (read definition)
#            pass

        # TODO this needs better parsing to include all possibilities and add them.
#        annotations = ""
#        if not pd.isna(row.FreeEnergy) and row.FreeEnergy:
#            temp = row.FreeEnergy.split()
#            annotations = f'dG0 {temp[0]} dG0_uncertainty {temp[1]} {annotations}'
#        if not pd.isna(row.EC) and row.EC:
#            annotations=f'ec-code {row.EC} {annotations}'
#        if len(annotations) > 0:
#            add_annotation( soup, new_reaction)
#            add_annotations( soup, new_reaction, annotations )

def add_list_parameters( soup, parameters ):
    """
    Add to the model a list of parameters.
    """
    model = soup.find('model')
    mylist = model.find('listOfParameters')
    if not mylist:
        new_tag = soup.new_tag('listOfParameters')
        model.append(new_tag)
        mylist = new_tag
    for label, value in zip(parameters[::2], parameters[1::2]):
        new_tag.append(soup.new_tag('parameter',attrs={
            'constant': 'true',
            'id': label,
            'sboTerm': 'SBO:0000626',
            'value': value
        }))

def add_list_objectives( soup, objectives ):
    """
    Add to the model a list of objectives as necessary.
    """
    # TODO: Do we really need this, and how is it parsed in the list? (read documentation)
    # Dummy
    if objectives:
        soup += objectives

def add_list_genes( soup, genes ):
    """
    Add to the model a list of genes as necessary.

    Dummy
    """
    if genes:
        for gene in genes:
            soup += gene

def unit_test():
    """
    Test the functions and structures in this file.
    """
    print("Running the tests...")
    print("Testing Reactions...")
    rxn.unit_test()
    print("Testing Compounds...")
    cpd.unit_test()
    print("Finished testing.")

def main():
    """
    The main routine for the network file
    """
    if len(sys.argv) > 1:
        if sys.argv[1] == 'test':
            unit_test()
            sys.exit()

    print("This file is used to define a network class")
    print("Run with argument 'test' to do tests.")

if __name__ == '__main__':
    main()
