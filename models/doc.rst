
SBML Models
===========

LUCA_v0.0.sbml
---------------
original version from article

LUCA.sbml
----------
Version December 2024 built from LUCA_v0.0.sbml by the manual addition of
various reactions and metabolites, essentially using CNapy, and the annotation
of reactions with DG° values as much as possible.

A more automated approach
-------------------------
Building this model and seeing its imperfections and the difficulties of
producing a complete model manually lead to the approach involving the different
models below that is more automated and easier (I hope) to maintain, document
and explain.

LUCA_full.sbml
---------------
This model is the union of LUCA_anabolism.sbml, LUCA_catabolism.sbml and
LUCA_environment.sbml described below using the merge_model.py program.

LUCA_anabolism.sbml
--------------------
The anabolism model contains the reactions for the synthesis of biomass From
the available materials, harvested from the environment. This model has been
constructed from various different elements, beyond those in LUCA.sbml and
LUCA_v0.0.sbml. This model is the union of multiple models included below and
contains xx reactions and yy species.

LUCA_environment.sbml
----------------------
The environment model contains the various exchanges between the presumed
environment of LUCA bringing raw materials and taking away waste products not
included in the biomass.

LUCA_catabolism.sbml
---------------------
The catabolism model contains the reactions necessary to provide the ATP
required my the anabolism model. Below there are several different anabolic
models based on different assumptions. The default model, that is here, is a
copy of the catabolic scheme driven by a respiratory chain and proton motive
force driven ATP synthase. Various alternative catabolic can be used to build
the full model.

----

FAB_pathway.sbml
-----------------
The fatty acid biosynthesis pathway, constructed from the tables `data/FAB_compound.csv`
and `data/FAB_reaction.csv` using `tools/table2sbml.py`. This pathway synthesises
hexadecanoyl-ACP from Acetyl-CoA using NADPH as the reducing source.

WL_pathway.sbml
----------------
The Wood-Ljungdahl pathway.
