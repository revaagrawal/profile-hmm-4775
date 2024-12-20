import argparse
import numpy as np
import pickle

"""
run command 
python classify_protein.py -m hmm_model.pkl -q PSTKLQEGGSVTMTCSSEEEEEEGLPAPEIFWSKKLDNGNLQHLGGGGGGGGGGGGGGGGGGGGGGGGGSGNAAAAAAAAATLTLIAMRMEDSGGIYVGCEGVNLIGKNRKEV
"""


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
    assuming all the probabilities are not in log-scale (this function converts to log-scale)

    Arguments:
        query: observed sequence of emitted states (string of a protein sequence)
        transition_probs: transition probabilities (2d np.array)
            note: matrix[row][col] where row is the previous state and col is the next state 
            and indexing is S, I_0, M_1, D_1, I_1, M_2, D_2, I_2
        m_emission_probs: match emission log-probabilities (2d np.array)
            each row corresponds to a match state and each column corresponds to an amino acid
            thus, m[1][A] gives the probability of emitting A in match state 1
            note: first col is all 0 bc M_0 does nto exist

        i_emission_probs: insertion emission log-probabilities (2d np.array)
        max_length: the length of the optimal sequence (hidden state sequence) produced by viterbi
    Returns:
        l: list of most likely hidden states at each position
            (list of hidden states)
        p: log-probability of the returned hidden state sequence
    """
    trans_probs = preprocess_probabilities(transition_probs)
    m_emiss_probs = preprocess_probabilities(m_emission_probs)
    i_emiss_probs = preprocess_probabilities(i_emission_probs)


    # print(f'trans: {trans_probs}')
    # print(f'm_emiss: {m_emiss_probs}')
    # print(f'i_emiss: {i_emiss_probs}')


    # print(f'dimensions trans: {trans_probs.shape}')
    # print(f'dimensions m_emiss: {m_emiss_probs.shape}')
    # print(f'dimensions i_emiss: {i_emiss_probs.shape}')


    n = len(query)
    amino_to_index = {
        "A": 0, "R": 1, "N": 2, "D": 3, "C": 4, "Q": 5, "E": 6, "G": 7, 
        "H": 8, "I": 9, "L": 10, "K": 11, "M": 12, "F": 13, "P": 14, 
        "S": 15, "T": 16, "W": 17, "Y": 18, "V": 19
    }

    # initialize, m represents the entire viterbi matrix (the best prob to get to the current state)
    # m[_][i] is the probability to observe a certain amino acid emitted when it is _ state at position i
    # m["I"][2] is prob that the query emitted the a.a. at position 2 due to insertion
    # the +1 is because there is this represents the start state S
    m = {
        "M": [float("-inf")] * (n + 1),
        "I": [float("-inf")] * (n + 1),
        "D": [float("-inf")] * (n + 1)
    }
    # print(f'm: {(m)}')


    # 0 represents the start state, which all have states have prob ln(1) = 0 because all begin at start
    for state in m:
        m[state][0] = 0

    # the first amino acid observed converted to its number form
    x_0 =amino_to_index[query[0]]

    # prob of starting with an insertion (I_0) and emitting x_0 (first amino acid observed in query)
    # prob = t(S -> I_0) + e(emitting x_0 due to insertion state I_0)
    # m["I"][1] = trans_probs[0][1] + i_emiss_probs[0][x_0]
    m["I"][1] = trans_probs[0][1] + i_emiss_probs[0][x_0]
    # print(f' t(S -> I_0): {trans_probs[0][1]}')
    # print(f' e(emitting x_0 due to insertion state I_0): {i_emiss_probs[0][x_0]}\n')


    # prob of starting with a match (M_1) and emitting x_0
    # prob = t(S -> M_1) + e(emitting x_0 due to match state M_1)
    m["M"][1] = trans_probs[0][2] + m_emiss_probs[1][x_0]
    print(f' t(S -> M_1): {trans_probs[0][2]}')
    print(f' e(emitting x_0 due to state M_1): {m_emiss_probs[1][x_0]}\n')
    print(f'm["M"][1] = {m["M"][1]}')


    # prob of starting with an deletion (D_1) and NO emission
    # prob = t(S -> D_1) 
    m["D"][1] = trans_probs[0][3] 


    # iterate through each amino acid in the query sequence (besides the very first one w/ index 0)
    for i in range(2, n+1):
        # the current emission/amino acid converted ot number form
        print(f'curr time step is t = {i}')
        x_i = amino_to_index[query[i-1]]
        
        # index for seq pos pos i (current) at each states
        curr_match_index = 2*i + 1
        curr_del_index = 2*i + 2
        curr_ins_index = 2*i + 3

        # index for seq i-1 (previous)at each states
        if (i <= 2):
            prev_match_index = 2
            prev_del_index = 3

            # when i-1 = 0-> I_0 has index 1, and when i-1=1 -> I_1 has index 4
            prev_ins_index = (i-1)*3 + 1


        # TODO: need to actually calc instead of making all same freq
        q_aa = np.log(1/20)

        # match
        # current emission term relative to background (in log term so subtration!)
        # print(f'M_{i}, x_i: {x_i} when the shape is {m_emiss_probs.shape}')
        if (m_emiss_probs[i][x_i] < q_aa):
            curr_match_emiss = m_emiss_probs[i][x_i] - q_aa
            print(f'match emission @ {i}: {m_emiss_probs[i][x_i]} + {q_aa}')
        else:
            curr_match_emiss = m_emiss_probs[i][x_i]
        # curr_match_emiss = m_emiss_probs[i][x_i]
        # print(f'curr_match_index: {curr_match_index}\n')

        v_m_to_m =  m["M"][i-1]  + trans_probs[prev_match_index][curr_match_index]
        print(f'viterbi at {i-1} at m: { m["M"][i-1]}')
        print(f'm{i-1} to m_{i} ({prev_match_index} to {curr_match_index}) : {trans_probs[prev_match_index][curr_match_index]}\n')
        print(f'm to m : { m["M"][i-1]} + { trans_probs[prev_match_index][curr_match_index]}\n')
        # print(f'm to m_{i}: { m["M"][i-1]  + trans_probs[prev_match_index][curr_match_index]}\n')
        v_i_to_m =  m["I"][i-1] + trans_probs[prev_ins_index][curr_match_index]
        print(f'i to m : { m["I"][i-1]} + { trans_probs[prev_ins_index][curr_match_index]}\n')
        v_d_to_m =  m["D"][i-1] +  trans_probs[prev_del_index][curr_match_index]
        print(f'd to m : {  m["D"][i-1]} + {trans_probs[prev_del_index][curr_match_index]}\n')


        print(f'prev {i-1} to match at {i}: {max(v_m_to_m, v_i_to_m, v_d_to_m)}\n')


        m["M"][i] = curr_match_emiss + max(v_m_to_m, v_i_to_m, v_d_to_m)


        # insertion. NOTE: that if you look at the graph, states that GO TO INSERTIONS
        # all have the same time step indexing as the insertion. but if you look at the other two
        # sections TO GO TO matches and deletions, the curr state is based on the previous time step

        if (i_emiss_probs[i][x_i] < q_aa):
            curr_ins_emiss = m_emiss_probs[i][x_i] - q_aa
        else:
            curr_ins_emiss = m_emiss_probs[i][x_i]

        v_m_to_i = m["M"][i-1] + trans_probs[curr_match_index][curr_ins_index]
        v_i_to_i = m["I"][i-1] + trans_probs[curr_ins_index][curr_ins_index]
        v_d_to_i = m["D"][i-1] + trans_probs[curr_del_index][curr_ins_index]

        m["I"][i] = curr_ins_emiss + max(v_m_to_i, v_i_to_i, v_d_to_i)


        # deletion, has no emission! does not have emission temr
        v_m_to_d = m["M"][i-1] + trans_probs[prev_match_index][curr_del_index]
        v_i_to_d = m["I"][i-1] + trans_probs[prev_ins_index][curr_del_index]
        v_d_to_d = m["D"][i-1] + trans_probs[prev_del_index][curr_del_index]

        m["D"][i] = max( v_m_to_d, v_i_to_d, v_d_to_d)

        print(f'm at time {i} is {m}\n')


    print(f'the viterbi matrix is the following: {m}')
    
    final_score = max(
        m["M"][-1],
        m["I"][-1],
        m["D"][-1]
    )
    # print(f)

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
    # print(f'transition: {transition_probs}')

    m_emission_probs = hmm_model["m_emission_probs"]
    # print(f'm emis: {len(m_emission_probs)}')

    i_emission_probs = hmm_model["i_emission_probs"]
    # print(f'i emis: {len(i_emission_probs)}')
    max_length = hmm_model["max_length"]

    score = score_sequence(args.q, transition_probs, m_emission_probs, i_emission_probs, max_length)
    print(f"Log-probability of query sequence: {score}")

    # TODO: decide what threshold would be best? not sure if there's a paper we can reference or 
    threshold = np.log(0.5)

    if score < threshold: 
        print(f"The queried protein sequence is likely a member of the IgSF family.")
    else:
        print(f"The queried protein sequence is unlikely to be a member of the IgSF family.")

if __name__ == "__main__":
    main()
