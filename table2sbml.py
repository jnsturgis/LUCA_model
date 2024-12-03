# Lets start with some libraries
# Find chemical formulae and charges seem to be missing and CHEBI
"""
This file should generate an initial sbml model based on Table1 and Table3 of Wimmers et al.
"""
from urllib import request
from bs4 import BeautifulSoup
import pandas as pd

col_list=[0,2]
PH = 7
TEMP = 30
col_list.append(int(3+14*(TEMP-25)/5+(PH-1)))

start1 = pd.read_excel("../Wimmers_Tables/Table 1.xlsx", header=4, index_col=0, usecols=[0,1,2,3,6])
start2 = pd.read_excel("../Wimmers_Tables/Table 3.XLSX",
                sheet_name='A', header=1, index_col=0, usecols=col_list)

reactions = pd.concat( [start1, start2], axis=1 )
reactions.rename({"Enzyme Commission number":"EC#", "pH 7.1":"dG0",
                "Reaction equation (Compound name) left":"R name",
                "Reaction equation (Compound name) right":"P name" }, axis=1, inplace=True)

reactions["R1"] = reactions["Reaction equation (KEGG compound)"].str.split("<=>").str[0]
reactions["P1"] = reactions["Reaction equation (KEGG compound)"].str.split("<=>").str[1]
reactions["Reagents"] = reactions["R1"].str.split("+")
reactions["Products"] = reactions["P1"].str.split("+")
reactions["R name"] = reactions["R name"].str.split(" ")
reactions["P name"] = reactions["P name"].str.split(" ")

reactions.drop("Reaction equation (KEGG compound)", axis=1, inplace=True)
reactions.drop("R1", axis=1, inplace=True)
reactions.drop("P1", axis=1, inplace=True)
reactions.drop("R name", axis=1, inplace=True)
reactions.drop("P name", axis=1, inplace=True)

start1.rename({"Enzyme Commission number":"EC#", "pH 7.1":"dG0",
                "Reaction equation (Compound name) left":"R name",
                "Reaction equation (Compound name) right":"P name" }, axis=1, inplace=True)
start1["R1"] = start1["Reaction equation (KEGG compound)"].str.split("<=>").str[0]
start1["P1"] = start1["Reaction equation (KEGG compound)"].str.split("<=>").str[1]

cpds = []
DIGITS = "0123456789"

def add_pair( mycode, myname ):
    """Build the structures"""
    if mycode[0] in DIGITS:
        assert mycode[0] == myname[0]
        mycode = mycode.lstrip(DIGITS)
        myname = myname.lstrip(DIGITS)
        mycode = mycode.strip()
        myname = myname.strip()

    for elem in cpds:
        if elem[0] == mycode and elem[1] == myname:
            return
    cpds.append([mycode, myname])

for ii in range(0,start1.shape[0]):
    codes = start1["R1"][ii].split(" + ")
    names = start1["R name"][ii].split(" + ")
    assert len(codes) == len(names)
    for code, name in zip(codes, names):
        add_pair(code.strip(),name.strip())

for ii in range(0,start1.shape[0]):
    codes = start1["P1"][ii].split(" + ")
    names = start1["P name"][ii].split(" + ")
    assert len(codes) == len(names)
    for code, name in zip(codes, names) :
        add_pair(code.strip(),name.strip())

compounds = pd.DataFrame( cpds, columns=['KEGG ID','Name'] )
print(compounds.shape)
compounds.to_csv("LUCAcompounds.csv")
reactions.to_csv("LUCAreactions.csv")

def fetch_info( species ):
    """Collect information from KEGG - this seems broken"""
    url = f"https://www.kegg.jp/entry/{species}"
    with request.urlopen(url).read().decode('utf8') as html:
        parsed_html = BeautifulSoup(html, "html.parser")
        table = parsed_html.find_all(class_ = "cel")[2]
    return ' '.join(table.stripped_strings)

formulae = []
charges = []
for index, compound in compounds.iterrows() :
    if index % 100 == 0 :
        print(compound)
    formulae.append( fetch_info(compound.to_list()[0]))
    charges.append( '0' )

compounds['Formula'] = formulae
compounds['Charge']  = charges

# Should fix these better
# Compound, Formula, Charge
modifications = [ ('C17023' , 'H2S', '0' ),             # Reduced sulfur donor
                 ('C00003' , 'C21H28N7O14P2', '1' ),    # NAD+
                 ('C00006' , 'C21H29N7O17P3', '1' ),    # NADP+
                 ('C01185' , 'C11H15NO9P', '1' ),       # NMN+
                 ('C00857' , 'C21H27N6O15P2', '1' ),    # deamino NAD+
                 ('C00080' , 'H', '1' ),                # H+
                 ('C01081' , 'C12H18N4O4PS', '1' ),     # TMP+
                 ('C00068' , 'C12H19N4O7P2S', '1' ),    # TPP+
                 ('C19609' , 'Ni', '2' ),               # Ni2+
                 ('C21512' , 'C42H53N6NiO14', '1' ),    # 15,17(3)-Seco-F430-17(3)-acid
                 ('C05777' , 'C42H51N6NiO13', '1' ),    # F430
                 ('C00175' , 'Co', '2' ),               # Co2+
                 ('C05773' , 'C45H59CoN4O14', '1' ),    # Cobyrinate
                 ('C06504' , 'C45H61CoN6O12', '1' ),    # Cob(II)yrinate diamide
                 ('C06506' , 'C55H73CoN11O15', '1' ),   # Adenosyl cobyrinate diamide
                 ('C06507' , 'C55H77CoN15O11', '1' ),   # Adenosyl cobyrinate hexaamide
                 ('C06508' , 'C58H84CoN16O11', '1' ),   # Adenosyl cobinamide
                 ('C06509' , 'C58H85CoN16O14P', '1' ),  # Adenosyl cobinamide phosphate
                 ('C06510' , 'C68H97CoN21O21P2', '1' ), # Adenosine-GDP-cobinamide
                 ('C11542' , 'C44H53CoN4O16', '1' ),    # Cobalt-precorrin 6
                 ('C11543' , 'C44H55CoN4O16', '1' ),    # Cobalt-dihydro-precorrin 6
                 ('C16244' , 'C44H57CoN4O14', '1' ),    # Cobalt-precorrin 7
                 ('C11545' , 'C45H59CoN4O14', '1' ),    # Cobalt-precorrin 8
                 ('C00288' , 'CHO3', '-1'),             # Bicarbonate anion
                 ('C05359' , '', '-1'),                 # Free electron
                 ('C01641' , 'C38H47N15O28P4', "-4"),   # Generic RNA composition PO4- AUCG
                 ('C02987' , 'C43H54N16O31P4', "-4"),   # Gluramyl-tRNA(Glu)
                 ('C00445' , 'C20H22N7O6', "1" ),       # +ve charge 5,10-Methenyltetrahydrofolate
                 ('C00138' , 'Fe2S2', '4'),             # Reduced ferredoxin is less positive than
                 ('C00139' , 'Fe2S2', '5'),             # Oxidised ferredoxin
                 ('C02745' , 'CH6', '0'),               # Reduced flavodoxin 2 H+ and 2 e- more than
                 ('C02869' , 'CH4', '0'),               # Oxidised flavodoxin (protein and FMN?)
                 ('C00030' , 'CH6', '0'),               # Reduced acceptor 2 H+ and 2 e- more than
                 ('C00028' , 'CH4', '0'),               # Oxidized acceptor
                 ('C04425' , 'C15H20N5O6S' , '1')       # S-Adenosyl-4-methylthio-2-oxobutanoate
               ]
for element in modifications:
    compounds.loc[ element[0], 'Formula'] = element[1]
    compounds.loc[ element[0], 'Charge'] = element[2]
# Should fix these better
# Now fix the reactions that need attention for Mass and charge balance mainly protons!!!
# 8716 add a H donor!
# 10397, 09153 add H donor and reduced sulfur donor H2S
# 01078 reduced sulfur donor is H2S
# RMAN4 remove a stray proton

modifications = [('R03231',['C00019', 'C01092', 'C00080'],['C04425', 'C01037']),
                ('R00742',['C00002', 'C00024', 'C00288', 'C00080'],
                ['C00008', 'C00009', 'C00083' ]),
                ('R00344',['C00002', 'C00022', 'C00288', 'C00080'],
                ['C00008', 'C00009', 'C00036' ]),
                ('R05220',['C06505', 'C00002', 'C00080'], ['C06506', 'C00536']),
                ('R07773',['C16243', 'C00019', 'C00080'], ['C11542', 'C00021']),
                ('R00575',['2 C00002', 'C00064', 'C00288', 'C00001', 'C00080'],
                ['2 C00008', 'C00009', 'C00025', 'C00169']),
                ('R05223',['C00194', 'C00144', 'C00080'], ['C06510', 'C05775']),
                ('R03348',['C03722', 'C00119', 'C00080'], ['C01185', 'C00013', 'C00011']),
                ('R07404',['C00002', 'C03373', 'C00288', 'C00080'], ['C00008', 'C00009', 'C15667']),
                ('R10712',['C04752', 'C20247', 'C00080'], ['C01081', 'C00013', 'C00011']),
                ('R12026',['C00003', 'C00037', 'C00283'],
                ['C00153', 'C20784', '3 C00001', 'C00080']),
                ('RMAN4', ['C00011', 'C00282'], ['C00058']),
                ('R11628',['C00002', 'C21511', '3 C00030', 'C00001', 'C00080'],
                ['C00008', 'C00009', 'C21512', '3 C00028']),
                ('RMAN1', ['C21107', 'C00014'], ['C20559', 'C00001']),
                ('R10397',['C17023', 'C16590', 'C00030'], ['C00001', 'C16593', 'C00028']),
                ('R09153',['C17023', 'C00593', 'C00030'], ['C00001', 'C03576', 'C00028']),
                ('R01078',['C01909', 'C17023', '2 C00019'], ['C00120', '2 C00073', '2 C05198']),
                ('R08716',['C17401', 'C00030'], ['C11539', 'C00028'])
                ]
for element in modifications:
#    print (element)
    print (reactions.loc[element[0]])
    new_element = reactions.loc[ element[0] ].to_list()
    new_element[3] = element[1]
    new_element[4] = element[2]
    reactions.loc[ element[0] ] = new_element

def species_id( species ) -> str:
    """Pretty silly function"""
    return f"i_{species}"

def species2sbml( species ) ->str:
    """Produce an sbml species reference"""
    if species[0] in DIGITS:
        number = 0
        i = 0
        while species[i] in DIGITS :
            number += int(species[i])
            i += 1
    else :
        number = 1
    answer = f"          <speciesReference species=\"{species_id(species.lstrip(DIGITS).strip())}\""
    answer += f" stoichiometry=\"{number}\" constant=\"true\"/>\n"
    return answer

def compound2sbml( my_index, my_compound) ->str :
    """Produce and sbml compound declaration"""
    # Consruct id from compartment and kegg ID
    my_id = species_id(my_index)

    output  = f"      <species metaid=\"meta_{my_id}\" id=\"{my_id}\" name=\"{my_compound[0]}\""
    output +=  " compartment=\"i\"\n"
    output +=  "          hasOnlySubstanceUnits=\"false\" boundaryCondition=\"false\""
    output +=    " constant=\"false\"\n"
    output +=  "          initialConcentration=\"1e-3\" substanceUnits=\"mole\"\n"
    output += f"          fbc:charge=\"{my_compound[2]}\" fbc:chemicalFormula=\"{my_compound[1]}\""
    output +=  ">\n"
    # Annotations
    output +=  "        <annotation>\n"
    output +=  "          <rdf:RDF xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\""
    output += " xmlns:dcterms=\"http://purl.org/dc/terms/\" xmlns:vCard=\"http://www.w3.org/2001"
    output += "/vcard-rdf/3.0#\" xmlns:vCard4=\"http://www.w3.org/2006/vcard/ns#\""
    output += " xmlns:bqbiol=\"http://biomodels.net/biology-qualifiers/\" xmlns:bqmodel=\""
    output += "http://biomodels.net/model-qualifiers/\">\n"
    output += f"            <rdf:Description rdf:about=\"#meta_{my_id}\">\n"
    output +=  "              <bqbiol:is>\n"
    output +=  "                <rdf:Bag>\n"
    output +=  "                  <rdf:li rdf:resource=\"https://identifiers.org/kegg/"
    output += f"KEGG:{my_compound[0]}\"/>\n"
    output +=  "                </rdf:Bag>\n"
    output +=  "              </bqbiol:is>\n"
    output +=  "            </rdf:Description>\n"
    output +=  "          </rdf:RDF>\n"
    output +=  "        </annotation>\n"
    output +=  "      </species>\n"
    return output

def reaction2sbml( my_index, my_reaction) ->str :
    """Produce an sbml reaction block"""
    output = f"      <reaction metaid=\"meta_{my_index}\" id=\"{my_index}\" name=\""
    output += f"KEGG {my_index}\""
    output += "          reversible=\"true\" fast=\"false\" fbc:lowerFluxBound=\"cobra_0_bound\""
    output += "fbc:upperFluxBound=\"cobra_default_ub\">\n"
    # Annotations
    output +=  "        <annotation>\n"
    output +=  "          <rdf:RDF xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" "
    output += "xmlns:dcterms=\"http://purl.org/dc/terms/\" xmlns:vCard=\"http://www.w3.org/2001"
    output += "/vcard-rdf/3.0#\" xmlns:vCard4=\"http://www.w3.org/2006/vcard/ns#\" xmlns:bqbiol"
    output += "=\"http://biomodels.net/biology-qualifiers/\" xmlns:bqmodel=\"http://biomodels"
    output += ".net/model-qualifiers/\">\n"
    output += f"            <rdf:Description rdf:about=\"#meta_{my_index}\">\n"
    output +=  "              <bqbiol:is>\n"
    output +=  "                <rdf:Bag>\n"
    output += "                  <rdf:li rdf:resource=\"https://identifiers.org/"
    output += f"dG0/{my_reaction[2]}\"/>\n"
    output += "                  <rdf:li rdf:resource=\"https://identifiers.org/"
    output += f"dG0_uncertainty/{my_reaction[1]}\"/>\n"
    output += "                  <rdf:li rdf:resource=\"https://identifiers.org/"
    output += f"EC #/{my_reaction[0]}\"/>\n"
    output +=  "                </rdf:Bag>\n"
    output +=  "              </bqbiol:is>\n"
    output +=  "            </rdf:Description>\n"
    output +=  "          </rdf:RDF>\n"
    output +=  "        </annotation>\n"
    # Reaction species and stochiometry
    output += "        <listOfReactants>\n"
    for species in my_reaction[3]:
        output += species2sbml( species.strip() )
    output += "        </listOfReactants>\n"
    output += "        <listOfProducts>\n"
    for species in my_reaction[4]:
        output += species2sbml( species.strip() )
    output += "        </listOfProducts>\n"
    output += "      </reaction>\n"
    return output

# pylint: disable="invalid-name"

SBML_HEADER = """<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core"
    xmlns:fbc="http://www.sbml.org/sbml/level3/version1/fbc/version2"
    metaid="meta_ML" id="ML" sboTerm="SBO:0000624" level="3" version="1" fbc:required="false">
  <model metaid="meta_LUCA" id="LUCA" fbc:strict="true">"""

SBML_UNITS = """      <listOfUnitDefinitions>
      <unitDefinition id="mmol_per_gDW_per_hr">
        <listOfUnits>
          <unit kind="mole" exponent="1" scale="-3" multiplier="1"/>
          <unit kind="gram" exponent="-1" scale="0" multiplier="1"/>
          <unit kind="second" exponent="-1" scale="0" multiplier="3600"/>
        </listOfUnits>
      </unitDefinition>
    </listOfUnitDefinitions>\n"""

SBML_COMPARTMENTS = """    <listOfCompartments>
      <compartment id="i" constant="true"/>
      <compartment id="e" constant="true"/>
    </listOfCompartments>"""

sbml_species  = "    <listOfSpecies>\n"
for index, compound in compounds.iterrows() :
    sbml_species += compound2sbml(index, compound.to_list())
sbml_species += "    </listOfSpecies>\n"

SBML_PARAMETERS = """    <listOfParameters>
      <parameter sboTerm="SBO:0000626" id="cobra_default_lb" value="-1000" constant="true"/>
      <parameter sboTerm="SBO:0000626" id="cobra_default_ub" value="1000" constant="true"/>
      <parameter sboTerm="SBO:0000626" id="cobra_0_bound" value="0" constant="true"/>
      <parameter sboTerm="SBO:0000626" id="minus_inf" value="-INF" constant="true"/>
      <parameter sboTerm="SBO:0000626" id="plus_inf" value="INF" constant="true"/>
      <parameter sboTerm="SBO:0000625" id="R_EX_CN_lower_bound" value="-100"
        units="mmol_per_gDW_per_hr" constant="true"/>
    </listOfParameters>\n"""

sbml_reactions  = "    <listOfReactions>"
for index, reaction in reactions.iterrows() :
    sbml_reactions += reaction2sbml(index.strip(), reaction.to_list())
sbml_reactions += "    </listOfReactions>\n"

SBML_OBJECTIVES = """    <fbc:listOfObjectives fbc:activeObjective="obj">
      <fbc:objective fbc:id="obj" fbc:type="maximize">
        <fbc:listOfFluxObjectives>
           <fbc:fluxObjective fbc:reaction="R_Z_BIOMASS" fbc:coefficient="1"/>
        </fbc:listOfFluxObjectives>
      </fbc:objective>
    </fbc:listOfObjectives>\n"""
SBML_FOOTER = """  </model>
</sbml>\n"""

with open('LUCA_v0.0.sbml', 'w+', encoding="utf-8") as fh:
    fh.write(SBML_HEADER)
    fh.write(SBML_UNITS)
    fh.write(SBML_COMPARTMENTS)
    fh.write(sbml_species)
    fh.write(SBML_PARAMETERS)
    fh.write(sbml_reactions)
#    fh.write(SBML_OBJECTIVES)
    fh.write(SBML_FOOTER)
