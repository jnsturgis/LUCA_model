## Branch Information ##

* main      - the main branch for merging other bits into.
* fba       - tools for building an appropriate objective function based on
              active parts of metabolism and required co-factors and initial
              concentrations. Also, and related, analysis of fba failures (what
              required species are not available and why) **fba.py**
* analysis  - Analysis of used and unused parts of the metabolism and
              thermodynamic choke points in the flux model.
* anabolism - develop separate anabolism and environment files and tools form
              merging them into a single sbml file.
* external  - develop tools for automatically recovering pertinent information
              from external databases - reactions, species, adjust pH for
              thermodynamic parameters.
* genetics  - improve catabolic model by semi-automatically adding genes so
              deletion studies are possible.
* reactions - add extra reactions for other parts that might/will be needed
              unwrap oleate biosynthesis and add quinones etc. (as far as
              possible, use scripts)

## Structure the Directory ##

* presentation - directory for latex document of article and presentation.
* data         - external data and sources with catalogue
* tools        - python tools for manipulating and analysing sbml models.
* models       - different models and parts of them with catalogue
