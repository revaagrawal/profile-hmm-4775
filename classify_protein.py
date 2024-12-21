import argparse
import numpy as np
import pickle


def preprocess_probabilities(matrix):
    """
    Converts a matrix to log-scale, ignoring 0 entries (keeps them as 0).
    """
    log_matrix = np.zeros_like(matrix) 
    nonzero_mask = matrix > 0
    log_matrix[nonzero_mask] = np.log(matrix[nonzero_mask]) 
    return log_matrix

def score_sequence(query, transition_probs, m_emission_probs, i_emission_probs, max_length):
    """
    Scores a query protein sequence against the HMM model using the Viterbi algorithm,
    assuming all the probabilities are in log-scale (similar to hw)
    """

    # TODO: fix that it doesn't work when converted to log-scale
    # transition_probs = np.log(transition_probs)
    # m_emission_probs = np.log(m_emission_probs)
    # i_emission_probs = np.log(i_emission_probs)

    transition_probs = preprocess_probabilities(transition_probs)
    m_emission_probs = preprocess_probabilities(m_emission_probs)
    i_emission_probs = preprocess_probabilities(i_emission_probs)

    n = len(query)
    amino_to_index = {
        "A": 0, "R": 1, "N": 2, "D": 3, "C": 4, "Q": 5, "E": 6, "G": 7, 
        "H": 8, "I": 9, "L": 10, "K": 11, "M": 12, "F": 13, "P": 14, 
        "S": 15, "T": 16, "W": 17, "Y": 18, "V": 19
    }

    m = {
        "M": [[float("-inf")] * (n + 1) for _ in range(max_length + 1)],
        "I": [[float("-inf")] * (n + 1) for _ in range(max_length + 1)],
        "D": [[float("-inf")] * (n + 1) for _ in range(max_length + 1)]
    }

    m["M"][0][0] = 0

    # iterate over all the match states
    for i in range(1, max_length + 1):
        # iterate over the pos in the query sequence
        for j in range(1, n + 1):
            amino = query[j - 1]
            amino_index = amino_to_index.get(amino, -1)

            # when a char in the query seq is not a valid amino acid
            if amino_index == -1:
                # TODO: throw an erro maybe?
                continue
            
            # match
            if 0 <= amino_index < m_emission_probs.shape[1]:
                m["M"][i][j] = max(
                    # M -> M, score = prev_match + weight going from prev match to curr match
                    m["M"][i - 1][j - 1] + transition_probs[i - 1][i], 
                    # I -> M
                    m["I"][i - 1][j - 1] + transition_probs[max_length + i - 1][i], 
                    # D -> M
                    m["D"][i - 1][j - 1] + transition_probs[2 * max_length + i - 1][i] 
                ) + m_emission_probs[i][amino_index]

            # insertion
            if 0 <= amino_index < i_emission_probs.shape[1]:
                m["I"][i][j] = max(
                    # M -> I
                    m["M"][i][j - 1] + transition_probs[i][max_length + i],  
                    # I -> I
                    m["I"][i][j - 1] + transition_probs[max_length + i][max_length + i] 
                ) + i_emission_probs[i, amino_index]
            
            # deletion
            m["D"][i][j] = max(
                # M -> D
                m["M"][i - 1][j] + transition_probs[i - 1][2 * max_length + i],
                # D -> D
                m["D"][i - 1][j] + transition_probs[2 * max_length + i - 1][2 * max_length + i] 
            )

    print(f'the viterbi matrix is the following: {m}')
    
    final_score = max(
        m["M"][max_length][n],
        m["I"][max_length][n],
        m["D"][max_length][n]
    )

    return final_score


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
