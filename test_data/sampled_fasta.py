from Bio import SeqIO
import random

def sample_fasta(input_fasta, output_fasta, num_samples):
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    
    sampled_sequences = random.sample(sequences, num_samples)
    
    SeqIO.write(sampled_sequences, output_fasta, "fasta")
    print(f"Sampled {num_samples} sequences and saved to {output_fasta}")

input_fasta = "protein-matching-PF00069.fasta"
output_fasta = "sampled_kinase.fasta"
num_samples = 111

sample_fasta(input_fasta, output_fasta, num_samples)