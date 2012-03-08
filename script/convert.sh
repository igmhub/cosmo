#!/bin/bash

# Builds ROOT files from a set of text files using text2root.
# Each file is converted <file>.txt => <file>.root[bs]
# Usually the input files are specified with a wildcard, e.g.
#
#   convert.sh results/fix_bao_no_noise_*.txt

for f in $@
do
    base=${f%.txt}
    text2root -i "$base.txt" -o "$base.root[bs]"
done
