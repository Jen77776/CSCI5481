from typing import Dict, List, Tuple
from question1and2 import (
    PhylogeneticAnalyzer,
    SequenceRepository,
    Sequence,
    GeneticDistanceStrategy,
    DistanceMatrix,
    TreeBuilder
)
import random
import os
import sys

class BootstrapAnalyzer:
    """
    Bootstrap Analyzer: Responsible for generating bootstrap samples and constructing corresponding phylogenetic trees
    """
    def __init__(self):
        """Initialize the bootstrap analyzer"""
        self._repository = SequenceRepository()
        self._distance_strategy = GeneticDistanceStrategy()
        
    def _generate_bootstrap_sequences(self, sequences: dict[str, Sequence]) -> dict[str, Sequence]:
        """
        Generate a single bootstrap sample
        
        Parameters:
            sequences: Original sequence dictionary
            
        Returns:
            Dict[str, Sequence]: Resampled sequence dictionary
        """
        if not sequences:
            return {}
            
        # Get the sequence length (assuming all sequences have the same length)
        seq_length = next(iter(sequences.values())).length
        
        # Randomly sample column indices (with replacement)
        columns = list(range(seq_length))
        sampled_columns = random.choices(columns, k=seq_length)
        
        # Create new bootstrap samples for each sequence
        bootstrap_seqs = {}
        for seq_id, sequence in sequences.items():
            # Convert the sequence to a list of characters for column-wise sampling
            seq_chars = list(sequence.content)
            # Create a new sequence based on the sampled column indices
            new_content = ''.join(seq_chars[i] for i in sampled_columns)
            bootstrap_seqs[seq_id] = Sequence(seq_id, new_content)
            
        return bootstrap_seqs
    
    def generate_bootstrap_trees(self, fasta_file: str, n_bootstrap: int, output_folder: str) -> None:
        """
        Generate a specified number of bootstrap trees
        
        Parameters:
            fasta_file: Path to the input FASTA file
            n_bootstrap: Number of bootstrap samples
            output_folder: Path to the output folder
        """
        # Create the output folder
        os.makedirs(output_folder, exist_ok=True)
        
        # Load the original sequence data
        self._repository.load_from_fasta(fasta_file)
        original_sequences = self._repository.sequences
        
        # Generate a tree for each bootstrap sample
        for i in range(n_bootstrap):
            # Create a new repository to store the bootstrap sample
            bootstrap_repo = SequenceRepository()
            
            # Generate a bootstrap sample
            bootstrap_seqs = self._generate_bootstrap_sequences(original_sequences)
            
            # Add the bootstrap sample to the repository
            for seq_id, sequence in bootstrap_seqs.items():
                bootstrap_repo._sequences[seq_id] = sequence
            
            # Calculate the distance matrix
            matrix = DistanceMatrix(bootstrap_repo, self._distance_strategy)
            
            # Build the tree
            builder = TreeBuilder(matrix.matrix, bootstrap_repo.sequence_ids)
            tree = builder.build()
            
            # Save the tree to a file
            output_file = os.path.join(output_folder, f'bootstrap_tree_{i+1}.txt')
            self._save_tree(tree, output_file)
    
    def _save_tree(self, tree, filepath: str) -> None:
        """
        Save the tree structure to a file
        
        Parameters:
            tree: List of tree edges
            filepath: Path to the output file
        """
        with open(filepath, 'w') as f:
            for ancestor, descendant, length in tree:
                f.write(f"{ancestor}\t{descendant}\t{length}\n")

def main():
    """Main function"""
    if len(sys.argv) != 4:
        print("Usage: python question6.py <fasta_file> <n_bootstrap> <output_folder>")
        sys.exit(1)
    
    try:
        fasta_file = sys.argv[1]
        n_bootstrap = int(sys.argv[2])
        output_folder = sys.argv[3]
        
        analyzer = BootstrapAnalyzer()
        analyzer.generate_bootstrap_trees(fasta_file, n_bootstrap, output_folder)
        
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()