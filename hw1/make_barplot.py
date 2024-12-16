import csv
import matplotlib.pyplot as plt

def read_csv(file_path):
    with open(file_path, 'r') as file:
        reader = csv.reader(file)
        return {row[0]: int(row[1]) for row in reader}

def combine_data(whole_genome_counts, separate_genes_counts):
    codons = set(whole_genome_counts.keys()) | set(separate_genes_counts.keys())
    combined_data = {codon: (whole_genome_counts.get(codon, 0), separate_genes_counts.get(codon, 0)) for codon in codons}
    return sorted(combined_data.items(), key=lambda x: x[1][1], reverse=True)

def create_barplot(combined_data):
    codons = [item[0] for item in combined_data]
    whole_genome_counts = [item[1][0] for item in combined_data]
    coding_sequences_counts = [item[1][1] for item in combined_data]

    fig, ax = plt.subplots(figsize=(20, 6))
    bar_width = 0.35
    x = range(len(codons))

    ax.bar([i - bar_width/2 for i in x], coding_sequences_counts, width=bar_width, label='Coding sequences (correct frame shift)', color='#8884d8')
    ax.bar([i + bar_width/2 for i in x], whole_genome_counts, width=bar_width, label='Whole genome (random frame shift)', color='#82ca9d')

    ax.set_ylabel('Frequency')
    ax.set_xticks(x)
    ax.set_xticklabels(codons, rotation=90)
    ax.legend(loc='upper right')

    plt.show()
    
if __name__ == "__main__":
    whole_genome_counts = read_csv('SARS-CoV-2_whole_genome_output.csv')
    separate_genes_counts = read_csv('SARS-CoV-2_separate_genes_output.csv')
    combined_data = combine_data(whole_genome_counts, separate_genes_counts)
    create_barplot(combined_data)
