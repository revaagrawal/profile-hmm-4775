def filter_psiblast_hits(input_file, output_file_prefix, thresholds):
    for threshold in thresholds:
        with open(input_file, "r") as infile, open(f"{output_file_prefix}_bit_{threshold}.txt", "w") as outfile:
            for line in infile:
                if line.startswith("#") or not line.strip():
                    outfile.write(line)
                    continue
                
                parts = line.split()
                if len(parts) < 12:
                    continue
                
                try:
                    bit_score = float(parts[11]) 
                    if bit_score >= threshold:
                        outfile.write(line)
                except ValueError:
                    continue

input_file = "psi_3_0.001.txt"
output_file_prefix = "psi_3_0.001"
thresholds = [0, 10, 50, 100]

filter_psiblast_hits(input_file, output_file_prefix, thresholds)