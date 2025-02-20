
Data files
==========

The different models have been produced from several data tables and pieces of
information. These tables and the data is found here and explained.

Parametrization
----------------

* Concentrations.csv` File with steady state concentrations of various
   metabolites and ions, based on E. coli from the bionumbers website.

This information formed the basis of the initial concentrations for Metabolites,
and the development of the biomass function.

Data tables for pathways
------------------------

* ``FAB_pathway.csv`` contain the data for producing the fatty acid bioynthesis pathway.
* ``WL_pathway.csv`` contain the data for producing the Wood-Ljungdahl carbon
assimilation pathway and energy production pathway, as per Buckel and Thauer 2013
following the Acetobacterium woodii energetics and reactions about 0.25 ATP/Acetate.
* ``WC_pathway.csv`` contain the data for producing the Wolfe cycle carbon
assimilation pathway. As shown in Buckel and Thauer about 0.5 ATP/Methane.
