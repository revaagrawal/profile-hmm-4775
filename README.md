# profile-hmm-4775


# To Run the Benchmarking:
### HMMer
`brew install hmmer`
- version 3.4

To build the profile HMM (from directory hmmer):
  `bash build_models.sh`

To search the profile HMM (from directory hmmer):
  `bash threshold_tests.sh`

### PSI-BLAST 
  `brew install blast`
  
  `makeblastdb -version`
  - makeblastdb: 2.16.0+

  `psiblast -version`
  - psiblast: 2.16.0+

  To build BLAST-compatible database:
  `makeblastdb -in training_set.fasta -dbtype prot -out target_db`

  To run PSI-BLAST, call:
  `psiblast -query complete_sequences.fasta -db target_db -num_iterations 3 -out [output txt file name] -outfmt 7 -evalue 0.001 -num_threads 4`
  
  - user-defined parameters to adjust: `num_iterations`, `evalue`, `num_threads`

  
# To Run the Reimplemention:
`pip install numpy`

To build the profile HMM,
Run through by the command:

  `python create_hmm.py -f [file containing the training sequences] -sigma 0.01`
- the arg -f is the training sequence (confirmed sequences in the igF family)
- the arg -sigma is the pseudocount

Resulting in a 'package' containing the transition probability matrix, the match emission probabilties, 
and the insertion emission probabilities saved into a file called hmm_model.pkl

To classify a file with query proteins, run:

  `python test_proteins.py -m hmm_model.pkl -f [file containing the query sequences]`
which then saves to results.csv

Alternatively, if the user wishes to manually test a protein sequence, the following can be run:

  `python classify_protein.py -m hmm_model.pkl -q [query protein]`
- the arg -q is the sequence of amino acids to be classified
