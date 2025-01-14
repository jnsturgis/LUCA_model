This directory contains a metabolic model of LUCA.
This model aims to be as complete as possible and provide a self consistent metabolism and objective functions for analysis.
Various choices will be necessary and explained here in the development of the model that is based on the article by Wimmer et al. 2021 (Energy at Origins: Favorable Thermodynamics of Biosynthetic Reactions in the Last Universal Common Ancestor (LUCA) Frontiers in Microbiology 12:793664 doi: 10.3389/fmicb.2021.793664)

The objective with this metabolic model is to be able to then propose an evolutionary path from earlier models while maintaining a coherent metabolic and catalytic scheme.

Step 1 Convert Wimmers model to an smbl model this is LUCA_clean.sbml that was created by the python code in ../Tables2SBML.ipynb convert to pure python code but fetching external info seems to be imperfect.

Step 2 Convert to a model suitable for FBA and a metabolic model:
	1. Need an objective function so add DNA, RNA, Protein, Membrane and cofactors to define the composition of LUCA. This needs some decisions!
	2. Need to define a catabolism to drive the anabolism and merge with the anabolic reations. This needs some decisions.
	3. Need to define exchanges with the environment to provide raw materials.

1.1 Following Jain et al. 2014 propose a membrane composed of a mixture of polyprenyl ether lipids on glycerol-1-phosphate and fatty acyl ester lipids on glycerol-3-phosphate forming a heterochiral membrane. These then have head groups added by the CDP route including serine, ethanolamine and glycerol.

1.2 Following from bionumbers (BNID 111490) and derived from Niedhardt et al. the composition of the cell (recalculated to remove LPS, peptidoglycan and glycogen) is about:
Protein 60% dry wt
RNA 22%
DNA 3%
Lipids 10%
Metabolites and cofactors 4%
Ions 1%

1.3 Within these we will opt for DNA and RNA 50%GC and 50%AT (an arbitrary choice)

1.4 Proteome composition based on data in Brune et al 2018.
W	1	Y	2	F	6	M	2
L	11	I	6	V	7.5	A	9.5
P	4	G	7.5	C	1	Q	4.5
N	4.5	T	5	S	5	E	6.5
D	5.5	K	4	H	2	R	5.5
For amino acid production arginine depends on acetate production and no recycling in initial pathways

1.5 Add membrane lipid biosynthesis pathway both acyl and isoprenyl - need to use ACP not phosphopantatheine and synth of oleyl-ACP and GGPP todo

1.6 Define metabolites and cofactors for growth
1.7 Correct define formulae and charges of different compounds and verify reactions are equilibrated.
1.8 Complete as possible external references dG° and KEGG and ec#

AcetylCoA C00024
AcetoacetylCoA C00332
Coenzyme A
HMGCoA C00356
Mevalonate C00418
Mevalonate-5-P C01107
Mevalonate-5-PP C01143
Isopentenyl-PP C00129
Dimethylallyl-PP C00235
GPP C00341
FPP C00448
GGPP C00353

Modifications made
==================
0. Add biomass function based on Bionumbers (BIOMASS)
1. Polymerisation and aggregation reactions
1.1 DNA polymerase
1.2 RNA polymerase
1.3 Protein synthesis from amino-acyl tRNA's
1.4 Lipid aggregation to make membrane
1.5 (TODO) Metabolite aggregation for cytoplasmic pool of cofactors and metabolites
2. Add reactions to produce membrane lipids, a mixture of archael and bacterial style lipids (Article)
2.1 Isoprenoid synthesis pathway (via mevalonate)
2.2 Fatty acid synthesis pathway
3. Some other bits:
3.1 CMP phosphorylation to CTP using ATP
3.2 Pyrophosphate hydrolysis
3.3 Hydrogenase to reduce NAD+
3.4 Acetyl CoA Synthase
3.5 ApoACP Synthase (RMAN8 2.7.8.7) and Adenosine 3'5' bis phosphate recycling
3.6 Iminoglycine synthesis R10245
3.7 FeS metabolism synthesis of iron sulfur clusters (IscS/IscU) Barras Article
3.8 Recycling/recovery of 5deoxyAdenosine Beaudoin article
4. Added exchange reactions (initially with internal compartment)
4.1 H2O, CO2, H+, NH3, H2S, PO4
4.2 Unlimited energy supply (EX_Energy) regenerates ATP for free!!!
4.3 Protein as Pr to simplify.
5. Remove unused reagents
5.1 Miscellaneous e-

Next steps.
1. Automatic recovery and tidying of annotations for reactions and metabolites
1.1 Data Sources
1.1.1 kegg https://www.genome.jp/entry/R00742 -> EC#, RHEA#
1.1.2 kegg https://www.genome.jp/entry/C00001 -> CHEBI, Mass
1.1.3 rhea https://www.rhea-db.org/rhea/11311 -> Directionality, Reagents charges, CHEBI, SMILES etc
1.2 Information desired
1.2.1 Reaction: EC, RHEA, KEGG, directionality, DG0', DG0' uncertainty
1.2.2 Metabolite: KEGG, CHEBI, Formula (fully protonated acids, deprotonated bases), charge pH7,
2. Addition of modifier species for reactions to indicate cofactors
2.1 Cofactors associated with enzyme. http://www.ebi.ac.uk/thornton-srv/databases/CoFactor/queries.php?ec=1.1.2.10&submit=Go
3. Setup of metabolites based on cofactors and cytoplasmic composition for Biomass function


TODO Estimate dG0 values uncertainty 1000 if not from equilibrator DONE
TODO lookup cofactor requirements (FeS included) for enzymes and add as modulator (check cnapy resistance) DONE (ish)
TODO Create metabolites CoA, NADH, ATP are in reactions from bionumbers data ? to BIOMASS reaction
dG0 for tRNA synthases +8 +/- 10 estimated from Transfer RNA Identity Contributes to Transition State Stabilization During Aminoacyl-tRNA Synthesis (1999) Nucleic Acids Research 27(18):3631-7

9/12/24
Next steps

1. Insert cofactors into reactions as modulators from a table (add_cofactors.py) DONE as far as possible try to check the litterature for the remaining reactions.
2. fba calculation protocol... (using cobra)
2.0 Aggregation reacions adjust stochiometry so Mr = 1500 => 1mM = 1.5mg/ml = 1%dry wt.
2.1 fba for protein, DNA, RNA and lipids (redo aggregations etc to make things easier to understand)
2.2 add to BIOMASS required cofactors for reactions with flux get contribution from annotations
2.3 add to BIOMASS elements from cycle with flux get contribution from annotations
2.4 redo fba and repeat from 2.2 until no changes.
3. work out how to add genes to the mix so can do a deletion analysis of the model...
4. work out thermodynamics from dG° and flux... where are bottlenecks and can we compensate with concentrations.
5. merge models to test different anabolic choices.

Note cnapy mangles:
1. Cofactors as modifiers
2. Initial concentrations
3. Notes

13/12/24
C00120 and C00194 cofactors not synthesised...

20/12/24
What we need to do and in what order

1. Separate Anabolism, Catabolism and Environment
1.1 Catabolism contains all species
1.2 Anabolism can add reactions with existing species to produce ATP, add genes, modify dG° and inital concentrations
1.2.1 Consider free energy from environment (meta_tri_phosphate)
1.2.2 Consider CO2 reduction pathway
1.2.3 Consider ATP synthase and respiratory chain with H2 donor and sulfite acceptor
1.3 Environment adds exchange reactions only
1.3.1 A meta_tri_phosphate environment
1.3.2 A CO2 acetate production environment
1.3.3 A more standard environement
1.4 A program to merge the components and check that the result is "complete"
1.4.1 Merge aspect - including biomass update based on modifiers (check outside fba how enzyme conc is supplied)
1.4.2 Check biomass production
1.4.3 Search for cycles in the reaction graph and add component to BIOMASS
1.4.4 Analysis of model failure and proposition for fixes (what exchanges make work - work back from BIOMASS reaction)

2. Modify the smbl files.
2.1 Semi-automatically add genes for reactions
2.1.1 Unwrap oleate biosynthesis
2.1.2 Based on EC number type labelling using bacterial and archael annotations
2.2 Modify initial concentrations for Bionumbers
2.3 Collect CHEBI and RHEA references automatically

3. Analysis of working model
3.1 Analysis of unused and unnecessary reactions and unreferenced substances
3.2 Analysis of thermodynamic choke points (possibly with initial concentration modification)

14/1/25: reactions.
table2sbml works and produces valid sbml.
Several points to resolve with the program
1. GeneRules and list of genes are not currently output.
2. Objectives are not output where to get them from?
3. Parametres are generic
4. Units are not output and generic
5. Not at all robust to errors in the tables.
