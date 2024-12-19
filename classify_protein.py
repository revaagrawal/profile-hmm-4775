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

    print(m_emission_probs[15])

    n = max_length
    max_score = float('-inf')
    start_index = None

    for start in range(len(query)):
        new_query = query[start:start+max_length]

        if len(new_query) < max_length and start != 0:
            break

        m = {
            "M": [float("-inf")] * (n + 1),
            "I": [float("-inf")] * (n + 1),
            "D": [float("-inf")] * (n + 1)
        }

        amino_to_index = {
            "A": 0, "R": 1, "N": 2, "D": 3, "C": 4, "Q": 5, "E": 6, "G": 7, 
            "H": 8, "I": 9, "L": 10, "K": 11, "M": 12, "F": 13, "P": 14, 
            "S": 15, "T": 16, "W": 17, "Y": 18, "V": 19
        }

        m["M"][0] = 0

        i = 0
        state = 0
        while i < min(max_length, len(query)):
            if i == 0:
                # start with match state
                i += 1
                state += 1
                continue

            # find next end_state
            for end_state in ["M", "I", "D"]:
                for prev_state in ["M", "I", "D"]:
                    # find transition probability
                    if end_state == "I":
                        end_index = (3 * state) + 1
                    elif end_state == "M":
                        end_index = (3 * state) + 2
                    elif end_state == "D":
                        end_index = (3 * state) + 3

                    if prev_state == "M":
                        prev_index = (3 * (state-1)) + 2
                    elif prev_state == "D":
                        prev_index = (3 * (state-1)) + 3
                    elif prev_state == "I":
                        prev_index = (3 * (state)) + 1
                    
                    transition_prob = transition_probs[prev_index][end_index]

                    # find emission probability
                    amino = query[i]
                    amino_acid = amino_to_index[amino]
                    if end_state == "M":
                        emission_prob = m_emission_probs[state+1][amino_acid]
                    elif end_state == "I":
                        emission_prob = i_emission_probs[state][amino_acid]
                    
                    prob = m[prev_state][i-1] + transition_prob + emission_prob

                    if prob > m[end_state][i]:
                        m[end_state][i] = prob
            i += 1

        i = len(m["M"]) - 1
        while m["M"][i] == float('-inf') and i > -1:
            i -= 1
        m_max = m["M"][i]

        i = len(m["I"]) - 1
        while m["I"][i] == float('-inf') and i > -1:
            i -= 1
        i_max = m["I"][i]

        i = len(m["D"]) - 1
        while m["D"][i] == float('-inf') and i > -1:
            i -= 1
        d_max = m["D"][i]

        # print(m_max, i_max, d_max)
        score = max(m_max, i_max, d_max)

        if score > max_score:
            max_score = score
            start_index = start
        
        # if start == 230:
        #     print(score)
        #     print(new_query)
        
    print("start", start_index)
    return max_score

    # m = {
    #     "M": [float("-inf")] * (n + 1),
    #     "I": [float("-inf")] * (n + 1),
    #     "D": [float("-inf")] * (n + 1)
    # }

    # amino_to_index = {
    #     "A": 0, "R": 1, "N": 2, "D": 3, "C": 4, "Q": 5, "E": 6, "G": 7, 
    #     "H": 8, "I": 9, "L": 10, "K": 11, "M": 12, "F": 13, "P": 14, 
    #     "S": 15, "T": 16, "W": 17, "Y": 18, "V": 19
    # }

    # m["M"][0] = 0

    # i = 0
    # state = 0
    # while i < min(max_length, n):
    #     if i == 0:
    #         # start with match state
    #         i += 1
    #         state += 1
    #     cur_state = None
    #     # find next end_state
    #     for end_state in ["M", "I", "D"]:
    #         for prev_state in ["M", "I", "D"]:
    #             # find transition probability
    #             if end_state == "I":
    #                 end_index = (3 * state) + 1
    #             elif end_state == "M":
    #                 end_index = (3 * state) + 2
    #             elif end_state == "D":
    #                 end_index = (3 * state) + 3

    #             if prev_state == "M":
    #                 prev_index = (3 * (state-1)) + 2
    #             elif prev_state == "D":
    #                 prev_index = (3 * (state-1)) + 3
    #             elif prev_state == "I":
    #                 prev_index = (3 * (state)) + 1
                
    #             transition_prob = transition_probs[prev_index][end_index]

    #             # find emission probability
    #             amino = query[i]
    #             amino_acid = amino_to_index[amino]
    #             if end_state == "M":
    #                 emission_prob = m_emission_probs[state+1][amino_acid]
    #             elif end_state == "I":
    #                 emission_prob = i_emission_probs[state][amino_acid]
                
    #             prob = m[prev_state][i-1] + transition_prob + emission_prob

    #             if prob > m[end_state][i]:
    #                 cur_state = end_state
    #                 m[end_state][i] = prob
    #     if cur_state != "D":
    #         i += 1

    # i = len(m["M"]) - 1
    # while m["M"][i] == float('-inf') and i > -1:
    #     i -= 1
    # m_max = m["M"][i]

    # i = len(m["I"]) - 1
    # while m["I"][i] == float('-inf') and i > -1:
    #     i -= 1
    # i_max = m["I"][i]

    # i = len(m["D"]) - 1
    # while m["D"][i] == float('-inf') and i > -1:
    #     i -= 1
    # d_max = m["D"][i]
    # score = max(m_max, i_max, d_max)
    # return score

    # m = {
    #     "M": [[float("-inf")] * (n + 1)],
    #     "I": [[float("-inf")] * (n + 1)],
    #     "D": [[float("-inf")] * (n + 1)]
    # }

    # m["M"][0] = 0
    # m["I"][0] = 0

    # i = 0
    # state = 1
    # while i < len(query):
    #     amino = query[i]
    #     amino_index = amino_to_index.get(amino, -1)

    #     if amino_index == -1:
    #         # TODO: throw an erro maybe?
    #         continue

    #     i_emission_p = i_emission_probs[state-1][amino_index]
    #     m_emission_p = m_emission_probs[state][amino_index]
    #     if i == 0: # initial state
    #         m["I"][state-1] = i_emission_p
    #         m["M"][state] = m_emission_p
    #     else:
            
    #     i += 1

    # # iterate over all the match states
    # for i in range(1, max_length + 1):
    #     # iterate over the pos in the query sequence
    #     for j in range(1, n + 1):
    #         amino = query[j - 1]
    #         amino_index = amino_to_index.get(amino, -1)

    #         # when a char in the query seq is not a valid amino acid
    #         if amino_index == -1:
    #             # TODO: throw an erro maybe?
    #             continue
            
    #         # match
    #         if 0 <= amino_index < m_emission_probs.shape[1]:
    #             m["M"][i][j] = max(
    #                 # M -> M, score = prev_match + weight going from prev match to curr match
    #                 m["M"][i - 1][j - 1] + transition_probs[i - 1][i], 
    #                 # I -> M
    #                 m["I"][i - 1][j - 1] + transition_probs[max_length + i - 1][i], 
    #                 # D -> M
    #                 m["D"][i - 1][j - 1] + transition_probs[2 * max_length + i - 1][i] 
    #             ) + m_emission_probs[i][amino_index]

    #         # insertion
    #         if 0 <= amino_index < i_emission_probs.shape[1]:
    #             m["I"][i][j] = max(
    #                 # M -> I
    #                 m["M"][i][j - 1] + transition_probs[i][max_length + i],  
    #                 # I -> I
    #                 m["I"][i][j - 1] + transition_probs[max_length + i][max_length + i] 
    #             ) + i_emission_probs[i, amino_index]
            
    #         # deletion
    #         m["D"][i][j] = max(
    #             # M -> D
    #             m["M"][i - 1][j] + transition_probs[i - 1][2 * max_length + i],
    #             # D -> D
    #             m["D"][i - 1][j] + transition_probs[2 * max_length + i - 1][2 * max_length + i] 
    #         )

    # print(f'the viterbi matrix is the following: {m}')
    
    # final_score = max(
    #     m["M"][max_length][n],
    #     m["I"][max_length][n],
    #     m["D"][max_length][n]
    # )

    # return final_score


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
