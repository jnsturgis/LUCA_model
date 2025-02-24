
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
* ``PS_pathway.csv`` the data for protein synthesis reactions from amino acids,
  making and average composition protein. Mr = 11248.39 so 1mM is 1.1% wt/vol.
* ``NAP_pathway.csv`` Nucleic acid polymerase reactions producing DNA (1204.76 g/mol)
  and RNA (1268.76 g/mol 1mM = 0.1269% wt/vol).
* ``Wimmers.csv`` contains the data for the reactions in the article by Wimmers
  et al. The Autotrophic Core: An Ancient Network of 404 Reactions Converts H2,
  CO2, and NH3 into Amino Acids, Bases, and Cofactors.
  Microorganisms 2021, 9, 458. https://doi.org/ 10.3390/microorganisms9020458
