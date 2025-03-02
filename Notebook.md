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

NOTE 5/1/25: Probably do not need FeS cluster exchange as synthesis is in model.

Note 5/1/25: These lines extract the LUCA.sbml catabolism (energy supply)
environment and anabolism segments respectively.
```
python3 tools/select_parts.py  models/LUCA.sbml -k ".*Energy.*"
python3 tools/select_parts.py  models/LUCA.sbml -k ".*exchange.*"
python3 tools/select_parts.py  models/LUCA.sbml -r ".*EX_.*"
```
14/1/25: reactions.
table2sbml works and produces valid sbml.
Several points to resolve with the program
1. GeneRules and list of genes are not currently output.
2. Objectives are not output where to get them from?
3. Parametres are generic
4. Units are not output and generic
5. Not at all robust to errors in the tables.

15/1/25
FAB_reaction.csv needs Dbases and FreeEnergy - should be filled by 'external'
FAB_compound.csv needs Dbases, charge and formula - should be filled by 'external'

16/1/25
Extract kegg-compound and kegg-reaction annotations from Id's automatically if there is
not already the annotation and the Id has the right form `Rx_{kegg-id}` or `Mx_{kegg-id}`.

19/1/25
Added Wood Ljungdahl pathway - and adjusted names in tables to conform to `{kegg-id}`
syntax introduced.

21/2/25
Need to find a way to use databases to automatically fill databases but still
use the previous knowledge

previous knowledge is:
reactions in original set
new reactions to add from kegg
modifications to make (protein, dna, rna synthesis)
transport reactions

24/2/25
```
network.py        - make a (networkx) graph for doing graphy things from .csv
kegg-reactions.py - takes list of reactions id's and recovers info from kegg.
```
Changed network.py to use ';' rather than ',' because of compound names with ','
Requires modification to csv input files (tr "," ";" < old > new)

1. built Wimmers.csv from Supplementary_Table_2
   `cut -f1 -d';' < Supplementary_Table_2.csv | sort | head --lines=-3 | tail -n +2 > toto`
    `python tools/kegg_reactions.py < data/toto > data/Wimmers.csv`
2. Manually add RMAN1 - 3 and methanofuran.
```
# Manual Modifications
Ri_RMAN1;Wimmers manually added reaction;Mi_C21107;Mi_C20559;
Ri_RMAN2;Wimmers manually added reaction;Mi_C20562;Mi_C05927;
Ri_RMAN3;Wimmers manually added reaction;Mi_C21070;Mi_C00862;
Mi_C00862;Methanofuran;1;0;C34H44N4O15
```
3. Add positive charges to (C00080)H+, (C00003)NAD+, (C00006)NADP+
4. Status with `python tools/table_verify.py < data/Wimmers.csv > toto` shows
   graph is in 4 pieces!!! 1 big one and:
	 * 'Ri_R00104' - this should not be
	 * 'Ri_R00127' - this should not be
	 * 'Ri_R00333', 'Mi_C00454', 'Mi_C00201'  - NDP and NTP only appearances but with AMP/ADP
	 * One reaction is charge inbalanced `Ri_R00189 unbalanced: charge -1.0`
   * 382 cycles found in graph.

TODO

Code the reaction names as Qxxxxx if changed even a little from the KEGG standard.
Setup Biomass (B), catabolism (K), Transport (T) and exchange (X) reactions to test.

Indicate modified metabolites as D rather than C  (changed equation) (this should
be done automatically later in the workflow)

Write network to SBML routines in network.py and possibly add construction for
cobra and fba.

Get DG'm from Equilibrator or elsewhere and fill in blanks unknown is not
allowed and then limit reversibility.

Identify missing reactions and decide on options:
* Add from kegg
* Add not from kegg

AUtomatic identification of blocking points.


25/2/25
Created
* Energy00.csv a simple ATP condensation reaction.
* Environment00.csv a series of exchange reactions.
* Biomass_full.csv an initial biomass reaction to get started.

Objective here is to get the combination Energy00, Environment00 and Biomass_full
to work together as a network and fba and its analysis to then be able to add
different aspects and keep working by backtracking.
* Membranes and lipid metabolisms
* Wimmers pathways
* FeS metabolism etc.

26/2/25
Thoughts on rolling back from LUCA.
* DNA is an afterthought
* Some amino acids are afterthoughts
* In rolling back consider protein evolution through duplication and specialization
  where did the enzymes/functionalities come from.
* A family tree of the enzymes for reactions.
* ATP synthase as an energy driven exporter or pump before takeover as ATP source.
* How did pyrophosphatase work as an ATP synthesiser in metabolism?

27/2/25
Thoughts on tools and documenting model building decisions.
1. Script `.sh` to build the model and produce figures for the modules and the
		whole model.
	* Central metabolism (EM, KC, RPPP) (no oxygen) (CM)
	* Proteins - amino acid biosynthesis (AAB) and protein polymerization (PS)
	* Nucleic acids - Purine (PUB), Pyrimidine (PYB) nucleotide biosynthesis and
	  polymerization reactions (NAP)
	* Membrane synthesis - Fatty acids (FAB), Polyprene (GGB), and membrane lipid
	  (ML) biosynthetic pathways.
	* Cofactor biosynthesis pathways and cycles probably several. (C01..C99)
	* Transport reactions (TR) for different transporters.
	* Anabolism, several options (A01..A99)
	* Exchange reactions for different scenarios and tests (X01..X99)
	* Biomass reactions for different scenarios and tests (B01..B99)
2. Model building step 1 start with list of kegg reaction id's produce `.csv`
		with `kegg_model.py`
3. Model building step 2 edit `.csv` file(s) with sed scripts.
4. Model building step 3 convert `.csv` files to `.sbml` models for different
    modules and check that they work and collect a minimal set of `cofactors`
		that appear in loops of the model.
5. Verify that the module `.sbml` files can produce flux in an appropriate
		environmental context.
6. Merge modules to form complete model as `.sbml` file.
7. Analysis:
	* fba
	* graph analysis
	* thermodynamic analysis

Thoughts on file formats and flow of information.
* `.rxn` List of kegg reaction id's as white space separated words.
* `.sed` sed control file for operating on `.csv` files.
* `.csv` condensed description of a model in csv format with ';' as the
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
* `.sbml` an sbml representation of the same data as in the csv file.

PUR.rxn from kegg modules M00048 M00049 M00050 all in Wimmers tables,
and kegg module M00053 with reactions R02014 and R02020 not in Wimmers tables.

PYR.rxn from kegg modules M00051 M00052, all in Wimmers tables except R00570 (ndk),
and kegg module M00938 and R02023, missing from Wimmers tables R02018, R02022,
R02023, R02024(desoxynucleotide synthases), R06613 alternative dUMP->dTMP and
R02331 (desoxy-ndk)

Note - this pathway would be simpler and less energy consuming if UMP could be
converted directly to dUMP or dUDP converted to dTDP. Even dUDP conversion to
dUMP would save energy.

AAB.rxn pathways from kegg modules
  * M00028 Ornithine biosynthesis without R02282 in Wimmers
	* M00844 Arginine biosynthesis in Wimmers
	* M00340 Histidine biosynthesis in Wimmers
	* M00020 Serine biosynthesis in Wimmers
	* M00018 Threonine biosynthesis in Wimmers
	* R00945 R00751 Glycine synthesis routes in Wimmers
	* R01001 R01290 Cysteine synthesis in Wimmers
  * R00243 R00248 R00253 R00258 R00396 R00400 R00355 R00483 R00578 A,D,N,E,Q Biosynthesis in Wimmers
	* M00019 2 oxobutanone biosynthesis in Wimmers
	* M00535 Valine/Isoleucine biosynthesis in Wimmers
	* M00432 Leucine biosynthesis in Wimmers
	* M00022 Shikimate biosynthesis pathway in Wimmers
	* M00025 Tyrosine biosynthesis pathway in Wimmers
	* M00024 Phenylalanine biosynthesis pathway in Wimmers
	* M00023 Tryptophan biosynthesis pathway in Wimmers
  * M00015 Proline biosynthesis pathway in Wimmers
  * M00016 Lysine biosynthesis pathway in Wimmers
	* M00030 Lysine biosynthesis pathway too in Wimmers
	* M00017 Methionine biosynthesis pathway without M04409 in Wimmers

CM.rxn based on kegg maps 10, 20, and 30
  * M00007 PPP Oxidative branch
	* M00006 PPP Non oxidative branch R01067 R08575 absent from Wimmers
	* M00005 PRPP synthesis in Wimmers
	* M00008 Entner-Doudoroff without R02036 R05605 in Wimmers
	* M00001 Glycolysis without R01786 R02189 R09085 R13199 R05805 R00756 R01068 R07159 in Wimmers
	* M00009 TCA cycle without R00621 R03316 R07618 R02570 R00727 R10343 R00361
	* R00709 in Wimmers adding R00351 R01899 R00268 absent from Wimmers

CF.rxn cofactor biosynthesis based on kegg map m01240
  * R10089 PyridoxalP synthesis
	* M00896 Archeal TPP synthesis
	* M00895 Prokaryotic TPP synthesis adding R05636 R07460 in Wimmers except
	  R07463 (added for iminoglycine synthesis)
	* M00115 NAD biosynthesis without R00481 (O2) adding R00104 for NADP as Wimmers
	* M00120 CoA biosynthesis in Wimmers
	* M00119 Pantothenate biosynthesis in Wimmers
  * M00573 BioI Biotin pathway without R10123 (O2) in Wimmers
	* M00572 PimolenoylACP synthesis in Wimmers
  * M00125 FAD biosynthesis
	* M00881 Lipoate biosynthesis NOT in Wimmers
	* M00126 THF biosynthesis plus R11072 in Wimmers
	* THF cycle in Wimmers
  * M00880 Molybdenum cofactor biosynthesis in Wimmers
	* M00121 Heme biosynthesis without R03220 R03222 R12605 R09489 R00310 but
	  with R04178 R11329 R11522. R03197 and R11329 absent from Wimmers.
  * M00846 Siroheme biosynthesis R2864 absent from Wimmers.
	* M00836 Cofactor F430 in Wimmers
	* M00924 Anaerobic Cobalamine biosynthesis in Wimmers
  * M00112 Cobalamin biosynthesis without R06529 as Wimmers but with R12184 R09083
	* M00989 Oxygen independant quinol synthesis absent from Wimmers
	* M00935 Methanofuran biosynthesis in Wimmers incomplete in kegg (RMAN)
	* M00378 Coenzyme F420 biosynthesis in Wimmers
	* M00358 Coenzyme M biosynthesis in Wimmers
  * M00608 Coenzyme B biosynthesis in Wimmers but 3 reactions missing R09720
	  R10391 R10394 in Wimmers and incomplete in kegg

After this effort...
3 reactions in Wimmers.rxn are not included in CM.rxn, CF.rxn, AAB.rxn, PUR.rxn
and PYR.rxn: RMAN1 RMAN2 RMAN3

Reactions not in Wimmers.rxn:
AAB.rxn:  R04405
CM.rxn:   R00268 R00351 R01067 R01899 R08575
CF.rxn:   R00424 R01218 R02864 R03197 *R04178* R04985 R05000 R06529 *R07463* *R08768*
					R08769 *R09083* R09720 R09737 R10391 R10394 R11329 R11522 *R12067* R12184
					R12423 R12424 R12427 R12428 R12657 *R12658* *R13426* R13439 R13440 R13441
PUR.rxn:  R00331 R01137 R01857 R02014 R02017 R02019 R02020
PYR.rxn:  R00570 R02018 R02022 R02023 R02024 R02331 R06613

*** I will need to assess how important additions are.

FAB.rxn: Synthesis of C16:0 AcylCoA, is not in Wimmers.rxn (except R01626).
GGB.rxn: Synthesis of GGPP (M00900), is not in Wimmers.rxn.
ML.rxn:  Synthesis of PS/PE on diacyl or digeranylgeranyl GolP, is not in Wimmers.rxn

Next steps:
1. Edit reaction lists to remove unnecessary reactions (especially if I added
   them without justification).
2. Create scripts to edit reactions to clarify abstractions notably Pr, Rn and
   Dn as units rather than 'R'

Convert to .csv
AAB: OK
CM:  OK
CF:  OK
FAB: OK
GGB: OK
ML:  OK
PUR: OK
PYR: OK

Remove oxygen requiring reactions ML R00846, AAB, FAB, CM, GGB, PUR, PYR,
CF R04178 R07463 R08768 R09083 R12067 R12658 R13426

Identify compounds with 'R' atoms...
AAB: LysW system C5H8NO4R and Holo Lys2 HSR and their adducts.
CM:  None
CF:  ACP HSR and adducts, Lipoyl carrier protein NH2R and adducts, Glycine cleavage system NH2R
FAB: ACP = HSR and adducts...
GGB: None
ML:  Acyl-CoA acyl group R
PUR: NTP/NDP R for nucleotide.
PYR: None

Delete reactions with NTP/NDP R00331 and R00333 from PYR.rxn (DONE)

Use this command to make the csv files from the rxn files.
`python ../tools/kegg_reactions.py CM.rxn CM.csv`

Used this script to make the modifications.
`for name in AAB CF CM FAB GGB ML PUR PYR; do
awk 'NR==FNR {if ($1 !~ /^#/) {rules[$1]=$2; replacements[$1]=$3} next}
     {for (key in rules) if ($0 ~ key) gsub(rules[key], replacements[key])}
     1' csv_editions.txt $name".csv" | tr "~" " " > $name"2.csv"
done
grep -vF 'Mi_C01641' < CF2.csv > tmp.csv; mv tmp.csv CF2.csv
grep -vF '^Mi_C02987' < CF2.csv > tmp.csv; mv tmp.csv CF2.csv
cat << EOT >>CF2.csv
Mi_C00046;Generic RNA;1;0;RnOH
Mi_C02987;Glutamyl-tRNA;1;0;C5H8O4NRn
EOT`

Next steps
1. make sure all reactions are balanced and fix in csv_editions if not!
    `python ../tools/csv_table_verify.py PYR2.csv` (not currently working -hangs on load)
2. work on `csv_table_verify.py` command line source and flags, writing sbml or csv
		and dealing with missing metabolites.

1/3/25

Unbalanced reactions:
AAB2.csv  OK
CF2.csv   Ri_R01078 Ri_R03231 Ri_R03348 Ri_R04109 Ri_R05220 Ri_R05223 Ri_R05578
					Ri_R07459 Ri_R07461 Ri_R07773 Ri_R08716 Ri_R09153 Ri_R09394 Ri_R10246
					Ri_R10397 Ri_R10404 Ri_R10712 Ri_R10802 Ri_R11580 Ri_R11628 Ri_R12026
					Ri_R12184 Ri_R12423 Ri_R12424 Ri_R13439 Ri_R13440 Ri_R13441
CM2.csv   Ri_R00344 Ri_R00742 Ri_R02164 Ri_R10866
FAB2.csv  Ri_R04534
GGB2.csv  OK
ML2.csv   Ri_R02242 Ri_R06872
PUR2.csv  Ri_R07404
PYR2.csv  R100Ri_R00575 Ri_R01868

Added to csv_editions.txt
Bicarbonate H
Flavodoxin formulae
Quinone formulae (fake)
Acceptor formulae (fake)
C22154 Iron cluster scaffold protein with 4Fe4S 2+ cluster... Fe4S4Pr
Remove H+ from substrates for Ri_R10092
Add H+ to products for Ri_R04534
C22154 Iron cluster scaffold protein with 4Fe4S 2+ cluster... Fe4S4Pr
Electron charge
Sulfur donor formula (fake)

Various incomplete reactions are problematic
Ri_R08716 Mi_C17401;Mi_C11539 added acceptor/donor
Ri_R09153 Mi_C00593;Mi_C03576 added H2S, reductant
Ri_R10397 Mi_C16590;Mi_C16593 added H2S, reductant
Ri_R10404 Mi_C00019;Mi_C00021 methylation by unknown mechanism

As are some metabolites that are "undefined" like acceptor and electron.

Need to keep an eye on these when I do the fba.

Next:
1. verify_table command line (not just stdin and flags for actions)
2. write_sbml from network.
3. read csv in network including extra fields.

then...
* setup environments for fba of different modules...
* find enzyme prosthetic groups and incorporate information... consider how, probably enzyme object.
