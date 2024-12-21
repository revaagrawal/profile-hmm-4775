# profile-hmm-4775


# To Run the Benchmarking:

todo


# To Run the Reimplemention:

To build the profile HMM,
Run through by the command:

  python create_hmm.py -f [file containing the training sequences] -sigma 0.01
- the arg -f is the training sequence (confirmed sequences in the igF family)
- the arg -sigma is the pseudocount

Resulting in a 'package' containing the transition probability matrix, the match emission probabilties, 
and the insertion emission probabilities saved into a file called hmm_model.pkl

To classify a file with query proteins, run:

  python test_proteins.py -m hmm_model.pkl -f [file containing the query sequences]
which then saves to results.csv

Alternatively, if the user wishes to manually test a protein sequence, the following can be run:

  python classify_protein.py -m hmm_model.pkl -q [query protein]
- the arg -q is the sequence of amino acids to be classified
