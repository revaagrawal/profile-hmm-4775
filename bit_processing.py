def filter_psiblast_hits(input_file, output_file_prefix, thresholds):
    """
    Filters PSI-BLAST hits by bit score thresholds.

    Args:
        input_file (str): Path to the PSI-BLAST output file in outfmt 7 format.
        output_file_prefix (str): Prefix for the filtered output files.
        thresholds (list): List of bit score thresholds to filter on.

    Outputs:
        Files with hits filtered by each threshold.
    """
    for threshold in thresholds:
        with open(input_file, "r") as infile, open(f"{output_file_prefix}_bit_{threshold}.txt", "w") as outfile:
            for line in infile:
                if line.startswith("#") or not line.strip():
                    # Write comments and empty lines as-is
                    outfile.write(line)
                    continue
                
                # Parse the tabular data
                parts = line.split()
                if len(parts) < 12:
                    continue  # Skip malformed lines
                
                try:
                    bit_score = float(parts[11])  # Bit score is the 12th column (index 11)
                    if bit_score >= threshold:
                        outfile.write(line)
                except ValueError:
                    continue  # Skip lines with non-numeric bit score

# Example usage
input_file = "psi_3_0.001.txt"
output_file_prefix = "psi_3_0.001"
thresholds = [0, 10, 50, 100]

filter_psiblast_hits(input_file, output_file_prefix, thresholds)