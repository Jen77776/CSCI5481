import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def smooth_conservation_scores(input_file, window_size=61):
    """
    Smooth the conservation scores using a sliding window.

    Args:
        input_file (str): Path to the input file containing conservation scores.
        window_size (int): Size of the sliding window, must be odd. Default is 61.

    Returns:
        pandas.Series: Smoothed conservation scores.
    """
    # Read conservation scores data
    conservation_scores = pd.read_csv(input_file, header=None, names=['score'])
    
    # Apply sliding window smoothing to the scores
    smoothed_scores = conservation_scores['score'].rolling(window=window_size, center=True).mean()
    
    return smoothed_scores

def plot_conservation_scores(scores, output_file):
    """
    Plot the conservation scores as a line chart.

    Args:
        scores (pandas.Series): Series containing conservation scores.
        output_file (str): Path to the output image file.
    """
    # Create figure and axis objects
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Plot conservation scores as a line chart
    ax.plot(scores, color='darkblue', linewidth=0.8, alpha=0.9)
    
    # Set plot title and axis labels
    ax.set_title('Conservation Scores')
    ax.set_xlabel('Position in Gapped Alignment')
    ax.set_ylabel('Conservation Score')
    
    # Adjust plot layout and save to file
    fig.tight_layout()
    plt.savefig(output_file)
    plt.close()

def main():
    # Input and output file paths
    input_file = 'solution-problem-1.txt'
    output_file = 'solution-problem-2.pdf'
    
    # Smooth conservation scores
    smoothed_scores = smooth_conservation_scores(input_file)
    
    # Plot smoothed conservation scores
    plot_conservation_scores(smoothed_scores, output_file)

if __name__ == '__main__':
    main()