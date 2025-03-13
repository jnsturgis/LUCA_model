#!/usr/bin/env bash

files="AAB CF CM FAB GGB ML PUR PYR";

# First steps convert .rxn to .csv
for f in $files; do
  echo python tools/kegg_reactions.py "data/"$f".rxn" "data/"$f".csv"
done

# Edit resulting .csv files
for name in files; do
  awk 'NR==FNR {if ($1 !~ /^#/) {rules[$1]=$2; replacements[$1]=$3} next}
     {for (key in rules) if ($0 ~ key) gsub(rules[key], replacements[key])}
     1' data/csv_editions.txt "data/"$name".csv" | tr "~" " " > "data/"$name"2.csv"
  rm "data/"$name".csv"
  mv "data/"$name"2.csv" "data/"$name".csv"
done
grep -vF  'Mi_C01641' < data/CF.csv > tmp.csv; mv tmp.csv data/CF.csv
grep -vF '^Mi_C02987' < data/CF.csv > tmp.csv; mv tmp.csv data/CF.csv
cat << EOT >>data/CF.csv
Mi_C00046;Generic RNA;1;0;RnOH
Mi_C02987;Glutamyl-tRNA;1;0;C5H8O4NRn
EOT

# Next steps verify reactions are balanced

# Merge csv with environments to test fba
 
# Produce sbml models
