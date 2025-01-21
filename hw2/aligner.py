import argparse

def read_fasta(file_path):
    with open(file_path, 'r') as f:
        header = f.readline().strip()
        sequence = ''.join(line.strip() for line in f)
    return header, sequence

def needleman_wunsch(seq1, seq2, match_score, mismatch_penalty, gap_penalty, ignore_outer_gaps=False):
    m, n = len(seq1), len(seq2)
    score_matrix = [[0] * (n + 1) for _ in range(m + 1)]
    
    # Initialize first row and column
    if not ignore_outer_gaps:
        for i in range(1, m + 1):
            score_matrix[i][0] = score_matrix[i-1][0] + gap_penalty
        for j in range(1, n + 1):
            score_matrix[0][j] = score_matrix[0][j-1] + gap_penalty
    
    # Fill the score matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score_matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty)
            delete = score_matrix[i-1][j] + (0 if ignore_outer_gaps and j == n else gap_penalty)
            insert = score_matrix[i][j-1] + (0 if ignore_outer_gaps and i == m else gap_penalty)
            score_matrix[i][j] = max(match, delete, insert)
    
    # Traceback
    align1, align2 = "", ""
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and (j == 0 or score_matrix[i][j] == score_matrix[i-1][j] + (0 if ignore_outer_gaps and j == n else gap_penalty)):
            align1 = seq1[i-1] + align1
            align2 = "_" + align2
            i -= 1
        elif j > 0 and (i == 0 or score_matrix[i][j] == score_matrix[i][j-1] + (0 if ignore_outer_gaps and i == m else gap_penalty)):
            align1 = "_" + align1
            align2 = seq2[j-1] + align2
            j -= 1
        else:
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
    
    return score_matrix[m][n], align1, align2

def format_alignment(align1, align2):
    return ''.join('|' if a == b and a != '_' else 'x' if a != '_' and b != '_' else ' ' for a, b in zip(align1, align2))

def main():
    parser = argparse.ArgumentParser(description='Needleman-Wunsch Global Alignment')
    parser.add_argument('-q', '--query', required=True, help='Query sequence file')
    parser.add_argument('-r', '--reference', required=True, help='Reference sequence file')
    parser.add_argument('-o', '--output', required=True, help='Output file')
    parser.add_argument('-g', '--gap_penalty', type=int, required=True, help='Gap penalty')
    parser.add_argument('-p', '--mismatch_penalty', type=int, required=True, help='Mismatch penalty')
    parser.add_argument('-m', '--match_score', type=int, required=True, help='Match score')
    parser.add_argument('--ignore_outer_gaps', action='store_true', help='Ignore gaps at the start and end of sequences')
    
    args = parser.parse_args()
    
    ref_header, seq1 = read_fasta(args.reference)
    query_header, seq2 = read_fasta(args.query)
    
    score, align1, align2 = needleman_wunsch(seq1, seq2, args.match_score, args.mismatch_penalty, args.gap_penalty, args.ignore_outer_gaps)
    
    with open(args.output, 'w') as f:
        f.write(f"{score}\n")
        f.write(f"{ref_header}\n")
        f.write(f"{align1}\n")
        f.write(f"{format_alignment(align1, align2)}\n")
        f.write(f"{align2}\n")
        f.write(f"{query_header}\n")

if __name__ == "__main__":
    main()