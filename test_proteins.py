import argparse
import numpy as np
import pickle
import classify_protein
import csv
import os

def read_fasta(filename):
    with open(filename, "r") as f:
        output = []
        s = ""
        header = ""
        for l in f.readlines():
            if l.strip()[0] == ">":
                if s != "":
                    output.append((header, s))
                    s = ""
                header = l
            else: s += l.strip()
        output.append((header, s))
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

    transition_probs = classify_protein.preprocess_probabilities(transition_probs)
    m_emission_probs = classify_protein.preprocess_probabilities(m_emission_probs)
    i_emission_probs = classify_protein.preprocess_probabilities(i_emission_probs)

    sequences = read_fasta(args.f)

    output_file = "results.csv"
    if not os.path.exists(output_file):
        with open(output_file, "w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(["header", "sequence", "log", "result", "expected", "match"])

    # stats: (e_p = 0 or 1, match: 1)

    for (header, seq) in sequences:
        try:
            print(header)
            print(seq)

            if len(seq) > 800:
                with open(output_file, "a", newline="") as file:  # Append to the file
                    writer = csv.writer(file)
                    writer.writerow([header, seq, None, None, expected_output, "too long", "too long"])
                continue

            expected_output = 0
            match = None

            expected = (header.strip())[-1]
            if expected == "P":
                expected_output = 1
            
            score = classify_protein.score_sequence(seq, transition_probs, m_emission_probs, i_emission_probs, max_length)
            print(f"Log-probability of query sequence: {score}")

            # TODO: decide what threshold would be best? not sure if there's a paper we can reference or
            threshold = -45

            if score > threshold: 
                print(f"The queried protein sequence is likely a member of the IgSF family.")
                if expected == "P":
                    match = 1
                elif expected == "N":
                    match = 0
            else:
                print(f"The queried protein sequence is unlikely to be a member of the IgSF family.")
                if expected == "P":
                    match = 0
                elif expected == "N":
                    match = 1
            
            with open(output_file, "a", newline="") as file:  # Append to the file
                writer = csv.writer(file)
                writer.writerow([header, seq, score, score > threshold, expected_output, match])
        except Exception as e:
            print(f"Error processing {header}: {e}")
            with open(output_file, "a", newline="") as file:  # Append to the file
                writer = csv.writer(file)
                writer.writerow([header, seq, None, None, "ERROR"])

if __name__ == "__main__":
    main()
