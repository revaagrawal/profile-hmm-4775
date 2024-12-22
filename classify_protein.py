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

    n = len(query)
    final_max_score = float('-inf')

    for start in range(len(query)):
        new_query = query[start:start+max_length]
        n = min(len(new_query), max_length)

        if len(new_query) < max_length and start != 0:
            break

        # i, j
        m = {
            "M": [[float("-inf")] * (n+1) for _ in range(max_length+1)], 
            "I": [[float("-inf")] * (n+1) for _ in range(max_length+1)],
            "D": [[float("-inf")] * (n+1) for _ in range(max_length+1)]
        }

        tb = {
            "M": [[float("-inf")] * (n+1) for _ in range(max_length+1)],
            "I": [[float("-inf")] * (n+1) for _ in range(max_length+1)],
            "D": [[float("-inf")] * (n+1) for _ in range(max_length+1)]
        }

        # initialize
        m["M"][0][0] = 0

        amino_to_index = {
            "A": 0, "R": 1, "N": 2, "D": 3, "C": 4, "Q": 5, "E": 6, "G": 7, 
            "H": 8, "I": 9, "L": 10, "K": 11, "M": 12, "F": 13, "P": 14, 
            "S": 15, "T": 16, "W": 17, "Y": 18, "V": 19
        }

        for i in range(1, n+1): # i == sequence number, not index (paper)
            query_index = i
            for j in range(max_length+1): # j == actual state num (paper)
                state_index = j
                # MATCH STATE
                if state_index > 0:
                    cur_index = (3 * (j-1)) + 2
                    amino = new_query[query_index-1]
                    amino_acid = amino_to_index[amino]
                    m_emission_p = m_emission_probs[j][amino_acid]
                
                    # from M
                    prev_index = cur_index - 3
                    transition_p = transition_probs[prev_index][cur_index]
                    prev_m = m["M"][state_index-1][query_index-1] + transition_p

                    # from I
                    prev_index = cur_index - 1
                    transition_p = transition_probs[prev_index][cur_index]
                    prev_i = m["I"][state_index-1][query_index-1] + transition_p

                    # from D
                    prev_index = cur_index - 2
                    transition_p = transition_probs[prev_index][cur_index]
                    prev_d = m["D"][state_index-1][query_index-1] + transition_p

                    max_score = max(prev_m, prev_i, prev_d)
                    m["M"][state_index][query_index] = max_score + m_emission_p

                    if max_score == prev_m:
                        tb["M"][state_index][query_index] = "M"
                    elif max_score == prev_i:
                        tb["M"][state_index][query_index] = "I"
                    elif max_score == prev_d:
                        # print(prev_m, prev_i, prev_d)
                        # print(m_emission_p)
                        tb["M"][state_index][query_index] = "D"

                # INSERTION STATE
                cur_index = (3 * j) + 1
                amino = new_query[query_index-1]
                amino_acid = amino_to_index[amino]
                emission_p = i_emission_probs[j][amino_acid]
                
                # from M
                if state_index == 0:
                    prev_m = 0
                else:
                    prev_index = cur_index - 2
                    transition_p = transition_probs[prev_index][cur_index]
                    prev_m = m["M"][state_index][query_index-1] + transition_p

                # from I
                if query_index == 0:
                    prev_i = 0
                else:
                    prev_index = cur_index
                    transition_p = transition_probs[prev_index][cur_index]
                    prev_i = m["I"][state_index][query_index-1] + transition_p

                # from D
                if state_index == 0:
                    prev_d = 0
                else:
                    prev_index = cur_index - 1
                    transition_p = transition_probs[prev_index][cur_index]
                    prev_d = m["D"][state_index][query_index-1] + transition_p

                max_score = max(prev_m, prev_i, prev_d)
                m["I"][state_index][query_index] = max_score + emission_p

                if max_score == prev_m:
                    tb["I"][state_index][query_index] = "M"
                elif max_score == prev_i:
                    tb["I"][state_index][query_index] = "I"
                elif max_score == prev_d:
                    tb["I"][state_index][query_index] = "D"

                # DELETE STATE
                if state_index > 0:
                    cur_index = (3 * (j-1)) + 3

                    # from M
                    prev_index = cur_index - 4
                    transition_p = transition_probs[prev_index][cur_index]
                    prev_m = m["M"][state_index-1][query_index] + transition_p

                    # from I
                    prev_index = cur_index - 2
                    transition_p = transition_probs[prev_index][cur_index]
                    prev_i = m["I"][state_index-1][query_index] + transition_p

                    # from D
                    prev_index = cur_index - 3
                    transition_p = transition_probs[prev_index][cur_index]
                    prev_d = m["D"][state_index-1][query_index] + transition_p

                    max_score = max(prev_m, prev_i, prev_d)
                    m["D"][state_index][query_index] = max_score

                    if max_score == prev_m:
                        tb["D"][state_index][query_index] = "M"
                    elif max_score == prev_i:
                        tb["D"][state_index][query_index] = "I"
                    elif max_score == prev_d:
                        tb["D"][state_index][query_index] = "D"
            
        max_m = m["M"][-1][-1]
        max_i = m["I"][-1][-1]
        max_d = m["D"][-1][-1]

        max_score = max(max_m, max_i, max_d)
        if max_score > final_max_score:
            final_max_score = max_score

    return final_max_score


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
