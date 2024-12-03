import argparse
import numpy as np

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
    
def get_insertion_columns(sequences):
    insertions = dict()
    for s in sequences:
        for i in range(len(s)):
            if s[i] == ".":
                if i in insertions:
                    insertions[i] += 1
                else:
                    insertions[i] = 1

    return insertions

def get_deletion_columns(sequences):
    deletions = dict()
    for s in sequences:
        for i in range(len(s)):
            if s[i] == "-":
                if i in deletions:
                    deletions[i] += 1
                else:
                    deletions[i] = 1

    return deletions

def remove_columns(cols, theta, total):
    col_keys = list(cols.keys())
    for col in col_keys:
        if cols[col] / total < theta:
            del cols[col]
    return cols

def get_state(base):
    state = None
    if base == "-":
        state = "D"
    else:
        state = "M"
    
    return state

def get_transition_probabilities(sequences, insertions, deletions, max_length):
    num = (max_length * 3) + 3
    print(num)
    table = np.zeros((num, num))

    # add pseudocounts
    for i in range(max_length):
        start_idx = (i*3) + 2
        for j in range(3):
            left = start_idx + j
            end_idx = ((i+1)*3)+1
            for end in range(end_idx, min(len(table), end_idx+3)):
                table[left][end] = 1

    for seq in sequences:
        prev_state = None
        start = 0
        for i in range(len(seq)):
            # find current state
            cur_state = None
            if i in deletions and seq[i] == ".":
                cur_state = "D"
            elif i in insertions and seq[i] != ".":
                cur_state = "I"
            elif seq[i] == ".":
                continue
            else:
                cur_state = "M"

            
            # determine transition
            if prev_state is None: # first
                if cur_state == "I":
                    table[0][1] += 1
                elif cur_state == "M":
                    table[0][2] += 1
                elif cur_state == "D":
                    table[0][3] += 1
            else:
                start_index = (3 * (start-1)) + 1
                end_index = (3 * start) + 1
            
                if prev_state == "I":
                    end_index -= 3
                    if cur_state == "M":
                        end_index += 1
                    elif cur_state == "D":
                        end_index += 2
                    elif cur_state == "I":
                        start -= 1

                elif prev_state == "M":
                    start_index += 1
                    if cur_state == "M":
                        end_index += 1
                    elif cur_state == "D":
                        end_index += 2
                    elif cur_state == "I":
                        start -= 1

                elif prev_state == "D":
                    start_index += 2
                    if cur_state == "M":
                        end_index += 1
                    elif cur_state == "D":
                        end_index += 2
    
                # print(seq, start_index, end_index, i)
                # print(seq[:i])
                # print(prev_state, cur_state)
                # print(start)
                table[start_index][end_index] += 1

            start += 1
            prev_state = cur_state

            if i == len(seq) - 1:
                start_index = (3 * (start-1)) + 1
                if cur_state == "M":
                    start_index += 1
                elif cur_state == "D":
                    start_index += 2
                elif cur_state == "I":
                    start_index += 3
                table[start_index][-1] += 1
        
    return table

def get_emission_probabilities(sequences, insertions, deletions, max_length):
    amino_acids = {
        "A": 0,
        "R": 1,
        "N": 2,
        "D": 3,
        "C": 4,
        "Q": 5,
        "E": 6,
        "G": 7,
        "H": 8,
        "I": 9,
        "L": 10,
        "K": 11,
        "M": 12,
        "F": 13,
        "P": 14,
        "S": 15,
        "T": 16,
        "W": 17,
        "Y": 18,
        "V": 19,
    }
    
    m_emissions = np.ones((max_length+1,20))
    m_emissions[:,0] = 0
    i_emissions = np.ones((max_length+1,20))

    for seq in sequences:
        state = 0
        for i in range(len(seq)):
            # find current state
            cur_state = None
            if i in deletions and seq[i] == ".":
                cur_state = "D"
            elif i in insertions and seq[i] != ".":
                cur_state = "I"
            elif seq[i] == ".":
                continue
            else:
                cur_state = "M"

            aa = seq[i]
            state += 1
            # update emission count
            if cur_state == "I":
                state -= 1
                i_emissions[state][amino_acids[aa]] += 1
            elif cur_state == "M":
                m_emissions[state][amino_acids[aa]] += 1
        
    return m_emissions, i_emissions 
    
def create_hmm():
    return
    
def main():
    parser = argparse.ArgumentParser(description='description')
    parser.add_argument('-f', action="store", dest="f", type=str, default='motif1.fasta')

    args = parser.parse_args()
    sequences = read_fasta(args.f)
    for s in sequences:
        print(s)

    total = len(sequences)

    dots = dict()
    for s in sequences:
        for i in range(len(s)):
            if s[i] == ".":
                if i in dots:
                    dots[i] += 1
                else:
                    dots[i] = 1
    insertions = set()
    deletions = set()
    for k in dots.keys():
        if dots[k] / total > 0.5:
            insertions.add(k)
        else:
            deletions.add(k)
    print(insertions)
    print(deletions)

    # TODO: already removed for "seed"? (lecture notes)
    # removed_cols = remove_columns(insertion_cols, 0.1, total)
    # print(removed_cols.keys())

    # filtered_sequences = filter_sequences(sequences, removed_cols)

    max_length = len(sequences[0]) - len(insertions)
    transition_matrix = get_transition_probabilities(sequences, insertions, deletions, max_length)

    row_sums = transition_matrix.sum(axis=1, keepdims=True)
    transition_probs = transition_matrix / row_sums

    m_emissions, i_emissions = get_emission_probabilities(sequences, insertions, deletions, max_length)
    row_sums = m_emissions.sum(axis=1, keepdims=True)
    m_emission_probs = m_emissions / row_sums

    row_sums = i_emissions.sum(axis=1, keepdims=True)
    i_emission_probs = i_emissions / row_sums

    for i in range(len(m_emission_probs)):
        print(i)
        print(i_emission_probs[i])

if __name__ == '__main__':
    main()
