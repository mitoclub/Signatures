#!/bin/bash

PATH_TO_MUTATIONS="./data/sample_mutations.tsv"
PATH_TO_MUTSPEC192="./data/sample_mutspec192.tsv"
PATH_TO_MUTSPEC96="./data/sample_mutspec96.tsv"
PATH_TO_SIGNATURES="./data/output"

# activate python environment (need to install 1 GB of libraries)
if [ -e env_signatures ]; then
    sourse env_signatures/bin/activate
else
    python3 -m venv env_signatures
    pip install -r requirements.txt
    sourse env_signatures/bin/activate
fi

if [ -e $PATH_TO_MUTATIONS ]; then
# calculate synonimous mutational spectra using mutations sample
python3 scripts/calculate_mutational_spectrum.py --data $PATH_TO_MUTATIONS --outfile $PATH_TO_MUTSPEC192 --syn

if [ -e $PATH_TO_MUTSPEC192 ]; then
# collapse 192-component mutationals spectrum to 96-component (COSMIC/Signal style) and 
# renormalize using human genome trinucleotides frequencies
python3 scripts/collapse_mutcspec.py $PATH_TO_MUTSPEC192 $PATH_TO_MUTSPEC96

if [ -e $PATH_TO_MUTSPEC96 ]; then
# extract signatures
python3 scripts/mut_signature_extraction.py -t 4 -m 5 $PATH_TO_MUTSPEC96 $PATH_TO_SIGNATURES

fi
fi
fi
