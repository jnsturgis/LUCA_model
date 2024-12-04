"""
Check that the reactions in sbml file are all balanced.
This will give a warning for non integer formulae - this is expected for the protein
metabolite. The unbalanced equations should just be the objectives, and the exchange 
reactions, these are not checked as they are in model.boundary
"""

import cobra as cb

model = cb.io.read_sbml_model("LUCA.sbml")

for reaction in set(model.reactions) - set(model.boundary):
    if len(reaction.check_mass_balance()) > 0 :
        print(reaction.id + str(reaction.check_mass_balance()) )
