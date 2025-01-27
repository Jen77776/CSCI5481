from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Dict, List, Tuple
import numpy as np
import sys

@dataclass
class Sequence:
    """
    Sequence dataclass: Represents a biological sequence
    
    Attributes:
        id (str): The unique identifier of the sequence
        content (str): The content string of the sequence
    """
    id: str
    content: str
    
    @property
    def length(self) -> int:
        """Returns the length of the sequence"""
        return len(self.content)

class SequenceRepository:
    """
    Sequence repository class: Responsible for reading, storing and managing sequence data
    
    Implements the repository pattern, centrally manages all sequence data, provides a unified access interface
    """
    def __init__(self):
        """Initialize an empty sequence dictionary"""
        self._sequences: Dict[str, Sequence] = {}
    
    def load_from_fasta(self, filepath: str) -> None:
        """
        Load sequence data from FASTA format file
        
        Parameters:
            filepath: The path of the FASTA file
            
        Exceptions:
            FileNotFoundError: Raised when the file does not exist
        """
        try:
            with open(filepath) as file:
                buffer = []
                current_id = None
                # Read FASTA file line by line
                for line in file:
                    if line.startswith('>'): # Start of a new sequence
                        if current_id: # Save the previous sequence
                            self._sequences[current_id] = Sequence(
                                current_id, 
                                ''.join(buffer)
                            )
                            buffer.clear()
                        current_id = line[1:].strip()# Remove '>' and store ID
                    else:
                        buffer.append(line.strip())# Collect sequence content
                # Save the last sequence
                if current_id and buffer:
                    self._sequences[current_id] = Sequence(
                        current_id,
                        ''.join(buffer)
                    )
        except FileNotFoundError:
            raise FileNotFoundError(f"File not found: {filepath}")
    
    @property
    def sequences(self) -> Dict[str, Sequence]:
        """Returns a copy of the sequence dictionary to prevent external modification"""
        return self._sequences.copy()
    
    @property
    def sequence_ids(self) -> List[str]:
        """Returns a list of all sequence IDs"""
        return list(self._sequences.keys())

class DistanceStrategy(ABC):
    """
    Distance calculation strategy interface
    
    Defines the abstract interface for sequence distance calculation, allows for different distance calculation methods to be implemented
    """
    @abstractmethod
    def calculate(self, seq1: str, seq2: str) -> float:
        """
        Calculate the distance between two sequences
        
        Parameters:
            seq1: The first sequence
            seq2: The second sequence
        Returns:
            float: The distance value between the two sequences
        """
        pass

class GeneticDistanceStrategy(DistanceStrategy):
    """
    Genetic distance calculation strategy implementation
    
    Calculates genetic distance based on the proportion of mismatched sites in the sequences
    """
    def calculate(self, seq1: str, seq2: str) -> float:
        """
        Calculate the genetic distance between two sequences
        
        Parameters:
            seq1: The first sequence
            seq2: The second sequence
        Returns:
            float: The genetic distance value (between 0-1)
        """
        if not (seq1 and seq2):
            return 0.0
            
        min_len = min(len(seq1), len(seq2))
        # Calculate the number of mismatched sites
        mismatches = sum(
            1 for a, b in zip(seq1, seq2)
            if a != b and not (a == '-' and b == '-') 
            and not (a == 'N' and b == 'N')
        )
        
        return mismatches / min_len

class DistanceMatrix:
    """
    Distance matrix class: Responsible for calculating and managing the distance matrix between sequences
    """
    def __init__(self, repository: SequenceRepository, 
                 strategy: DistanceStrategy):
        """
        Initialize the distance matrix
        
        Parameters:
            repository: The sequence repository instance
            strategy: The distance calculation strategy instance
        """
        self._repository = repository
        self._strategy = strategy
        self._matrix = None
        self._compute_matrix()
    
    def _compute_matrix(self) -> None:
        """Calculate the distance matrix between all sequence pairs"""
        sequences = self._repository.sequences
        ids = self._repository.sequence_ids
        size = len(ids)
        
        self._matrix = np.zeros((size, size))
         # Calculate the distance between each pair of sequences
        for i, id1 in enumerate(ids):
            for j, id2 in enumerate(ids):
                self._matrix[i, j] = self._strategy.calculate(
                    sequences[id1].content,
                    sequences[id2].content
                )
    
    def save(self, filepath: str) -> None:
        """
        Save the distance matrix to a file
        
        Parameters:
            filepath: The output file path
        """
        ids = self._repository.sequence_ids
        with open(filepath, 'w') as f:
            # Write the header
            f.write('\t' + '\t'.join(ids) + '\n')
            # Write the matrix data
            for i, id1 in enumerate(ids):
                row = [id1] + [f"{x}" for x in self._matrix[i]]
                f.write('\t'.join(row) + '\n')
    
    @property
    def matrix(self) -> np.ndarray:
        """Returns a copy of the distance matrix"""
        return self._matrix.copy()

class TreeBuilder:
    """
    Phylogenetic tree builder: Implements the Neighbor-Joining algorithm to construct the phylogenetic tree
    """
    def __init__(self, distance_matrix: np.ndarray, 
                 sequence_ids: List[str]):
        """
        Initialize the tree builder
        
        Parameters:
            distance_matrix: The distance matrix
            sequence_ids: The list of sequence IDs
        """
        self._matrix = distance_matrix.copy()
        self._names = sequence_ids.copy()
        self._n = len(sequence_ids) # Number of sequences
        self._size = self._n # Used for node numbering
        
    def build(self) -> List[Tuple[int, int, float]]:
        """
        Build the phylogenetic tree
        
        Use the Neighbor-Joining algorithm to construct an unrooted tree
        
        Returns:
            List[Tuple[int, int, float]]: The edge list of the tree, each edge contains (ancestor node, descendant node, branch length)
        """
        tree = []
        nodes = list(range(1, self._n + 1)) # Initial node list
        # Main loop: Continue until only two nodes remain
        while len(nodes) > 2:
            # Calculate the Q matrix
            row_sums = self._matrix.sum(axis=1)
            q_matrix = ((len(nodes) - 2) * self._matrix - 
                       row_sums[:, None] - row_sums)
            
            # Find the sequence pair corresponding to the minimum Q value
            i, j = self._find_min_pair(q_matrix)
            
           # Calculate the new internal node
            new_node = self._n + (self._size - len(nodes)) + 1
            limb_i, limb_j = self._calculate_limb_lengths(
                i, j, row_sums, len(nodes)
            )
            
            # Record the new branches
            tree.append((new_node, nodes[i], limb_i))
            tree.append((new_node, nodes[j], limb_j))
            
            # Update the distance matrix
            new_distances = np.zeros(len(nodes) - 2)
             # Update the node list
            remaining_indices = [k for k in range(len(nodes)) 
                               if k != i and k != j]
            
            for k, orig_k in enumerate(remaining_indices):
                new_distances[k] = 0.5 * (
                    self._matrix[i, orig_k] + 
                    self._matrix[j, orig_k] - 
                    self._matrix[i, j]
                )
            
            # Create a new distance matrix
            new_size = len(nodes) - 1
            new_matrix = np.zeros((new_size, new_size))
            
            # Copy the remaining distances
            for k1, orig_k1 in enumerate(remaining_indices):
                for k2, orig_k2 in enumerate(remaining_indices):
                    new_matrix[k1, k2] = self._matrix[orig_k1, orig_k2]
            
            # Add the distances for the new node
            for k in range(len(remaining_indices)):
                new_matrix[k, -1] = new_distances[k]
                new_matrix[-1, k] = new_distances[k]
            
            self._matrix = new_matrix
            
            # Update the node list
            nodes = [n for k, n in enumerate(nodes) 
                    if k not in (i, j)] + [new_node]
        
        # Add the last edge
        if len(nodes) == 2:
            tree.append((nodes[1], nodes[0], self._matrix[0, 1]))
        
        return tree
    
    def _find_min_pair(self, matrix: np.ndarray) -> Tuple[int, int]:
        """
        Find the row and column indices corresponding to the minimum value in the matrix
        
        Parameters:
            matrix: The matrix to search
        Returns:
            Tuple[int, int]: The (row, column) indices corresponding to the minimum value
        """
        n = matrix.shape[0]
        min_i, min_j = 0, 1
        min_val = float('inf')
        
        for i in range(n):
            for j in range(i + 1, n):
                if matrix[i, j] < min_val:
                    min_val = matrix[i, j]
                    min_i, min_j = i, j
        
        return min_i, min_j
    
    def _calculate_limb_lengths(self, i: int, j: int, 
                              row_sums: np.ndarray, 
                              n: int) -> Tuple[float, float]:
        """
        Calculate the branch lengths from the new node to its child nodes
        
        Parameters:
            i, j: The indices of the nodes to be merged
            row_sums: The row sums of the distance matrix
            n: The current number of nodes
        Returns:
            Tuple[float, float]: (length to i, length to j)
        """
        limb_i = 0.5 * (self._matrix[i, j] + 
                        (row_sums[i] - row_sums[j]) / 
                        (n - 2))
        limb_j = self._matrix[i, j] - limb_i
        return limb_i, limb_j

class PhylogeneticAnalyzer:
    """
    Phylogenetic analyzer: Integrates all components to perform a complete phylogenetic analysis
    """
    def __init__(self):
        self._repository = SequenceRepository()
        self._distance_strategy = GeneticDistanceStrategy()
    
    def analyze(self, fasta_file: str, 
               distance_output: str, 
               tree_output: str) -> None:
        """
        Execute the complete phylogenetic analysis workflow
        
        Parameters:
            fasta_file: The input FASTA file path
            distance_output: The output file path for the distance matrix
            tree_output: The output file path for the phylogenetic tree
        """
        # Load sequence data
        self._repository.load_from_fasta(fasta_file)
        
        # Calculate the genetic distance matrix
        matrix = DistanceMatrix(self._repository, self._distance_strategy)
        matrix.save(distance_output)
        
        # Build the phylogenetic tree
        builder = TreeBuilder(matrix.matrix, 
                            self._repository.sequence_ids)
        tree = builder.build()
        
        # Save the tree structure
        self._save_tree(tree, tree_output)
    
    def _save_tree(self, tree: List[Tuple[int, int, float]], 
                  filepath: str) -> None:
        """
        Save the tree structure to a file
        
        Parameters:
            tree: The edge list of the tree
            filepath: The output file path
        """
        with open(filepath, 'w') as f:
            for ancestor, descendant, length in tree:
                f.write(f"{ancestor}\t{descendant}\t{length}\n")

def main():
    """Main function"""
    if len(sys.argv) != 2:
        print("Usage: python3 phylogenetic_analysis.py <fasta_file>")
        sys.exit(1)
    
    try:
        analyzer = PhylogeneticAnalyzer()
        analyzer.analyze(
            sys.argv[1],
            "genetic-distances.txt",
            "edges.txt"
        )
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()