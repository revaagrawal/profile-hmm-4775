#!/usr/bin/env bash
# To run this script, change into the hmmer directory and execute `bash build_models.py` in terminal
# Note: need HMMer software

# FIRST MODEL - default parameters
hmmbuild --amino --seed 42 -o models/default_model.txt models/default_model.hmm msa_train.fasta

# SECOND MODEL - BLOSUM62 for weighting data, rest is default
hmmbuild --amino --seed 42 --wblosum -o models/blosum_model.txt models/blosum_model.hmm msa_train.fasta
 
# THIRD MODEL  - equally weighted data, +1 laplace prior 
hmmbuild --amino --seed 42 --wnone --enone --plaplace -o models/simple_model.txt models/simple_model.hmm msa_train.fasta

# Get the logos of these models
# hmmemit  -o logos/default_seq.fasta default_model.hmm
# weblogo -f logos/default_seq.fasta -o default_logo.png
