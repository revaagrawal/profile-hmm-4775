def process_fasta(input_files, output_file):
    with open(output_file, 'w') as out_file:
        for i, input_file in enumerate(input_files):
            # positive (P) for first file, negative (N) for second and third file
            label = 'P' if i == 0 else ('NF' if i == 1 else 'NK')
            
            with open(input_file, 'r') as in_file:
                lines = in_file.readlines()
                for line in lines:
                    if line.startswith('>'):
                        header = line.strip() + f"|{label}\n"
                        out_file.write(header)
                    else:
                        out_file.write(line)

if __name__ == "__main__":
    input_files = ['test_sequences.fasta', 'sampled_fibropectin3.fasta', 'sampled_kinase.fasta']
    output_file = 'complete_sequences.fasta'

    process_fasta(input_files, output_file)
    print(f"Concatenated sequences with labels saved to {output_file}")