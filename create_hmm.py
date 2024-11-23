import argparse


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

def remove_columns(cols, theta, total):
    col_keys = list(cols.keys())
    for col in col_keys:
        if cols[col] / total < theta:
            del cols[col]
    return cols
    
def create_hmm():
    return
    
def main():
    parser = argparse.ArgumentParser(description='description')
    parser.add_argument('-f', action="store", dest="f", type=str, default='motif1.fasta')
    # parser.add_argument('-k', action="store", dest="k", type=int, default=10)
    # parser.add_argument('-epsilon', action="store", dest="epsilon", type=float, default=1.0)

    args = parser.parse_args()
    sequences = read_fasta(args.f)
    for s in sequences:
        print(s)

    total = len(sequences)
    insertion_cols = get_insertion_columns(sequences)
    print(insertion_cols.keys())

    removed_cols = remove_columns(insertion_cols, 0.1, total)
    print(removed_cols.keys())
    # k = args.k
    # epsilon = args.epsilon

    # starts, pi = gibbs_sampling(sequences, k, epsilon)
    # print("The starting positions for the motif (starts) in each sequence are:")
    # print(starts)
    # print("============================")
    # print("The motif model (pi) is:")
    # print(np.matrix(pi))
    # motif_sequences = get_motif_sequences(sequences, starts, k)
    # print("============================")
    # print("The motifs in each sequence (motif_sequences) are:")
    # print('\n'.join([("Sequence %d: " % i) + m for i, m in enumerate(motif_sequences)]))
    # print("============================")
    # print("The consensus motif is: %s" % majority(motif_sequences))
    # print("============================")
    # # Change pi to log scale
    # for j in range(k):
    #     for b in range(4):
    #         pi[j][b] = np.log(pi[j][b])
    # print("The log likelihood of the motif model is:")
    # print(get_log_liklihood(sequences, starts, pi, k))

if __name__ == '__main__':
    main()
