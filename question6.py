import csv
from collections import defaultdict

genetic_code = {
    'AAA':'Lys', 'AAC':'Asn', 'AAG':'Lys', 'AAT':'Asn',
    'ACA':'Thr', 'ACC':'Thr', 'ACG':'Thr', 'ACT':'Thr',
    'AGA':'Arg', 'AGC':'Ser', 'AGG':'Arg', 'AGT':'Ser', 
    'ATA':'Ile', 'ATC':'Ile', 'ATG':'Met', 'ATT':'Ile',
    'CAA':'Gln', 'CAC':'His', 'CAG':'Gln', 'CAT':'His',
    'CCA':'Pro', 'CCC':'Pro', 'CCG':'Pro', 'CCT':'Pro',
    'CGA':'Arg', 'CGC':'Arg', 'CGG':'Arg', 'CGT':'Arg',
    'CTA':'Leu', 'CTC':'Leu', 'CTG':'Leu', 'CTT':'Leu',
    'GAA':'Glu', 'GAC':'Asp', 'GAG':'Glu', 'GAT':'Asp',
    'GCA':'Ala', 'GCC':'Ala', 'GCG':'Ala', 'GCT':'Ala',
    'GGA':'Gly', 'GGC':'Gly', 'GGG':'Gly', 'GGT':'Gly',
    'GTA':'Val', 'GTC':'Val', 'GTG':'Val', 'GTT':'Val',
    'TAA':'Stp', 'TAC':'Tyr', 'TAG':'Stp', 'TAT':'Tyr',
    'TCA':'Ser', 'TCC':'Ser', 'TCG':'Ser', 'TCT':'Ser',
    'TGA':'Stp', 'TGC':'Cys', 'TGG':'Trp', 'TGT':'Cys',
    'TTA':'Leu', 'TTC':'Phe', 'TTG':'Leu', 'TTT':'Phe'
}

def count_amino_acids(csv_file):
    amino_acid_counts = defaultdict(int)
    
    with open(csv_file, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            codon, count = row
            amino_acid = genetic_code[codon]
            amino_acid_counts[amino_acid] += int(count)
    
    return amino_acid_counts


whole_genome_counts = count_amino_acids('SARS-CoV-2_whole_genome_output.csv')
separate_genes_counts = count_amino_acids('SARS-CoV-2_separate_genes_output.csv')

import matplotlib.pyplot as plt
def create_barplot(whole_genome_counts, separate_genes_counts):
    amino_acids = sorted(separate_genes_counts, key=separate_genes_counts.get, reverse=True)
    whole_genome_freqs = [whole_genome_counts.get(aa, 0) for aa in amino_acids]
    separate_genes_freqs = [separate_genes_counts[aa] for aa in amino_acids]
    
    fig, ax = plt.subplots(figsize=(12, 8))
    bar_width = 0.35
    x = range(len(amino_acids))
    
    ax.bar([i - bar_width/2 for i in x], separate_genes_freqs, width=bar_width, label='Separate Genes', color='#1f77b4')
    ax.bar([i + bar_width/2 for i in x], whole_genome_freqs, width=bar_width, label='Whole Genome', color='#63a8e2')
    
    ax.set_ylabel('Frequency')
    ax.set_title('Amino Acid Frequency Comparison')
    ax.set_xticks(x)
    
    ax.set_xticklabels(amino_acids, rotation='vertical')
    
    ax.legend()
    
    plt.tight_layout()
    plt.savefig('amino_acid_frequency_barplot.png')
    plt.show()
    
create_barplot(whole_genome_counts, separate_genes_counts)