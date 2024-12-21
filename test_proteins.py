import argparse
import numpy as np
import pickle
import classify_protein

def read_fasta(filename):
    with open(filename, "r") as f:
        output = []
        s = ""
        for l in f.readlines():
            if l.strip()[0] == ">":
                if s == "": continue
                output.append(s)
                s = ""
            else: s += l.strip()
        output.append(s)
        return output

def main():
    parser = argparse.ArgumentParser(description="Evaluate Sequence Membership to Profile HMM")
    parser.add_argument("-m", required=True, type=str, help="Path to saved HMM model")
    parser.add_argument('-f', action="store", dest="f", type=str, default='data/complete_sequences.fasta')
    args = parser.parse_args()

    # opens thre Python object (profile HMM with the transition, emissions, and length) made from 'create_hmm.py'
    with open(args.m, "rb") as f:
        hmm_model = pickle.load(f)

    transition_probs = hmm_model["transition_probs"]
    m_emission_probs = hmm_model["m_emission_probs"]
    i_emission_probs = hmm_model["i_emission_probs"]
    max_length = hmm_model["max_length"]

    sequences = read_fasta(args.f)
    for seq in sequences:
        print(seq)
        score = classify_protein.score_sequence(seq, transition_probs, m_emission_probs, i_emission_probs, max_length)
        print(f"Log-probability of query sequence: {score}")

        # TODO: decide what threshold would be best? not sure if there's a paper we can reference or
        threshold = -45

        if score > threshold: 
            print(f"The queried protein sequence is likely a member of the IgSF family.")
        else:
            print(f"The queried protein sequence is unlikely to be a member of the IgSF family.")
        
        print("-----------------------------------------------")

if __name__ == "__main__":
    main()
