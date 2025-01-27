import os
import sys
from collections import defaultdict
from typing import Dict, List, Set, Tuple

class TreeAnalyzer:
    """Tree Analyzer: Calculate bootstrap support"""
    
    def __init__(self):
        """Initialize the tree analyzer"""
        self.node_to_leaves = {}  # Store leaf nodes under each internal node
    
    def _get_leaves_below_node(self, edges: List[Tuple[int, int, float]], 
                              node: int) -> Set[int]:
        """
        Get all leaf nodes below the specified node
        
        Parameters:
            edges: List of tree edges [(ancestor, descendant, length),...]
            node: Node ID to search for
            
        Returns:
            Set of all leaf nodes below the specified node
        """
        leaves = set()
        # Find direct child nodes
        children = [(d, l) for a, d, l in edges if a == node]
        
        for child, _ in children:
            # If the child node is a leaf node (ID <= number of leaves)
            if not any(a == child for a, _, _ in edges):
                leaves.add(child)
            else:
                # Recursively get leaf nodes of the subtree
                leaves.update(self._get_leaves_below_node(edges, child))
        
        return leaves
    
    def _read_tree_file(self, filepath: str) -> List[Tuple[int, int, float]]:
        """Read tree file and return the list of edges"""
        edges = []
        with open(filepath) as f:
            for line in f:
                ancestor, descendant, length = line.strip().split('\t')
                edges.append((int(ancestor), int(descendant), float(length)))
        return edges
    
    def _get_node_partitions(self, edges: List[Tuple[int, int, float]]) -> Dict[int, Tuple[Set[int], Set[int]]]:
        """
        Get leaf node partition for each internal node
        
        Returns:
            Dictionary {internal node ID: (left subtree leaf set, right subtree leaf set)}
        """
        partitions = {}
        # Get all internal nodes
        internal_nodes = set(a for a, _, _ in edges)
        # Remove leaf nodes
        internal_nodes = set(n for n in internal_nodes 
                           if any(a == n for a, _, _ in edges))
        
        for node in internal_nodes:
            # Get direct child nodes
            children = [(d, l) for a, d, l in edges if a == node]
            if len(children) >= 2:  # Ensure at least two child nodes
                # Get leaf nodes for each subtree
                child_leaves = [self._get_leaves_below_node(edges, child) 
                              for child, _ in children[:2]]  # Take first two child nodes
                partitions[node] = (child_leaves[0], child_leaves[1])
        
        return partitions
    
    def calculate_support(self, true_tree_file: str, bootstrap_folder: str) -> None:
        """
        Calculate bootstrap support
        
        Parameters:
            true_tree_file: Path to the original tree file
            bootstrap_folder: Path to the bootstrap tree folder
        """
        # Read the original tree
        true_edges = self._read_tree_file(true_tree_file)
        # Get partitions of the original tree
        true_partitions = self._get_node_partitions(true_edges)
        
        # Counter for the number of occurrences of each partition
        partition_counts = defaultdict(int)
        
        # Read all bootstrap tree files
        bootstrap_files = [f for f in os.listdir(bootstrap_folder) 
                          if f.startswith('bootstrap_tree_')]
        total_bootstrap_trees = len(bootstrap_files)
        
        # Analyze each bootstrap tree
        for bootstrap_file in bootstrap_files:
            filepath = os.path.join(bootstrap_folder, bootstrap_file)
            bootstrap_edges = self._read_tree_file(filepath)
            bootstrap_partitions = self._get_node_partitions(bootstrap_edges)
            
            # Check if each partition of the original tree appears in this bootstrap tree
            for node, (left, right) in true_partitions.items():
                # Check if this partition appears in any node of the bootstrap tree
                for _, (boot_left, boot_right) in bootstrap_partitions.items():
                    if (left == boot_left and right == boot_right) or \
                       (left == boot_right and right == boot_left):
                        partition_counts[node] += 1
                        break
        
        # Calculate and save support values
        with open('bootstrap.txt', 'w') as f:
            for node in sorted(true_partitions.keys()):
                support = partition_counts[node] / total_bootstrap_trees
                f.write(f"{node}\t{support:.2f}\n")

def main():
    if len(sys.argv) != 3:
        print("Usage: python question7.py edges.txt bootstrap_trees_folder")
        sys.exit(1)
    
    try:
        true_tree_file = sys.argv[1]
        bootstrap_folder = sys.argv[2]
        
        analyzer = TreeAnalyzer()
        analyzer.calculate_support(true_tree_file, bootstrap_folder)
        
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()