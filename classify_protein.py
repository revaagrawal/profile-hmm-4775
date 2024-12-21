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
    max_score = float('-inf')

    # for start in range(len(query)):
    for start in range(1):
        new_query = query[start:start+max_length]

        if len(new_query) < max_length and start != 0:
            break

        # i, j
        m = {
            "M": [[float("-inf")] * (n+1) for _ in range(max_length+1)], # 0th column ALWAYS 0
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
                    amino = query[query_index-1]
                    amino_acid = amino_to_index[amino]
                    emission_p = m_emission_probs[j][amino_acid]
                
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
                    m["M"][state_index][query_index] = max_score + emission_p

                # INSERTION STATE
                cur_index = (3 * j) + 1
                amino = query[query_index-1]
                amino_acid = amino_to_index[amino]
                emission_p = i_emission_probs[j][amino_acid]
                
                # from M
                prev_index = cur_index - 2
                transition_p = transition_probs[prev_index][cur_index]
                prev_m = m["M"][state_index][query_index-1] + transition_p

                # from I
                prev_index = cur_index
                transition_p = transition_probs[prev_index][cur_index]
                prev_i = m["I"][state_index][query_index-1] + transition_p

                # from D
                prev_index = cur_index - 1
                transition_p = transition_probs[prev_index][cur_index]
                prev_d = m["D"][state_index][query_index-1] + transition_p
                max_score = max(prev_m, prev_i, prev_d)
                m["I"][state_index][query_index] = max_score + emission_p

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
            
        max_m = m["M"][-1][-1]
        max_i = m["I"][-1][-1]
        max_d = m["D"][-1][-1]

        print(max_m)
        print(max_i)
        print(max_d)
        max_score = max(max_m, max_i, max_d)
        return max_score

        # # initialize
        # for i in range(n):
        #     m["M"][i][0] = 0
        #     m["M"][i][1] = 0

        # for i in range(1, n):
        #     for j in range(2, n):
        #         # match state 
        #         end_index = (3 * (j-1)) + 2

        #         prev_index = (3 * (j-2)) + 2
        #         prev_m = m["M"][i-1][j-1] + transition_probs[prev_index][end_index]

        #         prev_index = (3 * (j-1)) + 1
        #         prev_i = m["I"][i-1][j-1] + transition_probs[prev_index][end_index]

        #         prev_index = (3 * (j-2)) + 3
        #         prev_d = m["D"][i-1][j-1] + transition_probs[prev_index][end_index]

        #         max_score = max(prev_m, prev_i, prev_d)
                
        #         amino = new_query[i]
        #         amino_acid = amino_to_index[amino]
        #         emission_prob = m_emission_probs[j][amino_acid]

        #         m["M"][i][j] = max_score + emission_prob

        #         # insertion state
        #         end_index = (3 * j) + 1

        #         prev_index = (3 * (j-1)) + 2
        #         prev_m = m["M"][i-1][j] + transition_probs[prev_index][end_index]

        #         prev_index = end_index
        #         prev_i = m["I"][i-1][j] + transition_probs[prev_index][end_index]

        #         prev_index = (3 * (j-1)) + 3
        #         prev_d = m["D"][i-1][j] + transition_probs[prev_index][end_index]

        #         max_score = max(prev_m, prev_i, prev_d)

        #         emission_prob = i_emission_probs[j][amino_acid]

        #         m["I"][i][j] = max_score + emission_prob

        #         # delete state
        #         end_index = (3 * (j-1)) + 3
    
        #         prev_index = (3 * (j-2)) + 2
        #         prev_m = m["M"][i][j-1] + transition_probs[prev_index][end_index]

        #         prev_index = (3 * (j-1)) + 1
        #         prev_i = m["I"][i][j-1] + transition_probs[prev_index][end_index]

        #         prev_index = (3 * (j-2)) + 3
        #         prev_d = m["D"][i][j-1] + transition_probs[prev_index][end_index]

        #         max_score = max(prev_m, prev_i, prev_d)

        #         m["D"][i][j] = max_score

        # m_score = m["M"][-1][-1]
        # i_score = m["I"][-1][-1]
        # d_score = m["D"][-1][-1]

        # print(m_score, i_score, d_score)

        # result_score = max(m_score, i_score, d_score)
        # return m_score


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
