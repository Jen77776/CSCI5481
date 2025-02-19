def calculate_conservation_rates(msa):
    """
    Calculate the conservation rate at each position in a multiple sequence alignment.
    
    Args:
        msa (list): A list of strings, each string representing a sequence in the MSA.
    
    Returns: 
        list: A list of conservation rates, one per position in the alignment.
    """
    num_sequences = len(msa)
    alignment_length = len(msa[0])
    
    conservation_rates = []
    
    # Iterate over each position in the alignment
    for position in range(alignment_length):
        # Get all bases at this position
        bases_at_position = [seq[position].upper() for seq in msa]
        
        # Count the occurrences of each valid base (ACGT, ignoring gaps and invalid bases)
        valid_bases = "ACGT"
        base_counts = {base: bases_at_position.count(base) for base in valid_bases}
        
        # Find the count of the most common valid base
        max_base_count = max(base_counts.values()) if base_counts else 0
        
        # Calculate conservation rate (denominator includes all characters, including gaps)
        conservation_rate = max_base_count / num_sequences
        conservation_rates.append(conservation_rate)
    
    return conservation_rates


def main():
    # Read in the alignment file 
    alignment_file = "Homework4-seqs-with-primers.fna"
    with open(alignment_file, "r") as file:
        lines = file.readlines()

    # Extract sequences from the lines
    sequences = []
    current_sequence = ""
    for line in lines:
        if line.startswith(">"):
            if current_sequence:
                sequences.append(current_sequence)
                current_sequence = ""
        else:
            current_sequence += line.strip()
    if current_sequence:
        sequences.append(current_sequence)

    # Calculate conservation rates
    conservation_rates = calculate_conservation_rates(sequences)

    # Write conservation rates to file
    output_file = "solution-problem-1.txt"
    with open(output_file, "w") as file:
        for rate in conservation_rates:
            file.write(f"{rate:.6f}\n")


if __name__ == "__main__":
    main()