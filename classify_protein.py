import argparse
import numpy as np
import pickle

def score_sequence(query, transition_probs, m_emission_probs, i_emission_probs, max_length):
    """
    Scores a query protein sequence against the HMM model using the Viterbi algorithm
    """
    
    return 0

def main():
    parser = argparse.ArgumentParser(description="Evaluate Sequence Membership to Profile HMM")
    parser.add_argument("-q", required=True, type=str, help="Query protein sequence")
    parser.add_argument("-m", required=True, type=str, help="Path to saved HMM model")
    args = parser.parse_args()

    # opens thre Python object (profile HMM with the transition, emissions, and length) made from 'create_hmm.py'
    with open(args.m, "rb") as f:
        hmm_model = pickle.load(f)

    transition_probs = hmm_model["transition_probs"]
    m_emission_probs = hmm_model["m_emission_probs"]
    i_emission_probs = hmm_model["i_emission_probs"]
    max_length = hmm_model["max_length"]

    score = score_sequence(args.q, transition_probs, m_emission_probs, i_emission_probs, max_length)
    print(f"Log-probability of query sequence: {score}")

    # TODO: decide what threshold would be best? not sure if there's a paper we can reference or
    threshold = 0

    if score > threshold: 
        print(f"The queried protein sequence is likely a member of the IgSF family.")
    else:
        print(f"The queried protein sequence is unlikely to be a member of the IgSF family.")

if __name__ == "__main__":
    main()
