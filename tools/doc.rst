

Programme documentation
=======================

Some logic in the naming conventions:

Programmes that begin with `kegg_` operate on a text file containing words that
are might be interpreted as kegg identifiers.
Programmes that operate on `csv` files containing compound and reaction data start
with `csv_`
Programmes that operate on `sbml` files containing a valid model, either an fbc
model or not, start with `sbml_`

`kegg_reactions.py` uses the set of identifiers, identified as reaction in kegg,
to interrogate the kegg database and produces a csv file containing compound and
reaction data.

`kegg_logic.py` does logical operations on two sets of identifiers.

check_balance
-------------

.. automodule:: check_balance
  :members:

fba
---

.. automodule:: fba
  :members:

merge_models
------------

.. automodule:: merge_models
  :members:

rename_model
------------

.. automodule:: rename_model
  :members:

sbml2tables
-----------

.. automodule:: sbml2tables
  :members:

table2sbml
----------

.. automodule:: table2sbml
  :members:
