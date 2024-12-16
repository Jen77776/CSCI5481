import sys
import csv
from collections import Counter

def count_codons(input_file):
    codon_counts = Counter()
    stop_codons = ['TAA', 'TAG', 'TGA']
    stop_codon_positions = []
    first_atg_position = -1
    position = 1
    with open(input_file, 'r') as file:
        sequence = ""
        for line in file:
            if line.startswith('>'):
                if sequence:
                    codons = [sequence[i:i+3] for i in range(0, len(sequence) - len(sequence) % 3, 3)]
                    codon_counts.update(codons)
                    for i, codon in enumerate(codons):
                        current_position = position + i * 3
                        if codon in stop_codons and len(stop_codon_positions) < 10:
                                stop_codon_positions.append(current_position)
                        if codon == 'ATG' and first_atg_position == -1:
                            first_atg_position = current_position
                    position += len(sequence)                            
                    sequence = ""
            else:
                sequence += line.strip().upper()
        
        if sequence:
            codons = [sequence[i:i+3] for i in range(0, len(sequence) - len(sequence) % 3, 3)]
            codon_counts.update(codons)
            for i, codon in enumerate(codons):
                current_position = position + i * 3
                if codon in stop_codons and len(stop_codon_positions) < 10:
                        stop_codon_positions.append(current_position)
                if codon == 'ATG' and first_atg_position == -1:
                    first_atg_position = current_position
        print(f"Positions of the first 10 stop codons: {stop_codon_positions}")
        print(f"Position of the first ATG codon: {first_atg_position}")
    
    return codon_counts

def write_to_csv(codon_counts, output_file):
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(codon_counts.items())

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    codon_counts = count_codons(input_file)
    write_to_csv(codon_counts, output_file)