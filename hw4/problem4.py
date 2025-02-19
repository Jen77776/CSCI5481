import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plot_conservation_rates(data_file, regions_file, output_file):
    """
    Plot the smoothed conservation rates and highlight variable regions.

    Parameters:
    data_file (str): Path to the file containing conservation rates data.
    regions_file (str): Path to the file containing variable regions data.
    output_file (str): Path to save the generated plot.
    """
    # Read and smooth conservation rates data
    conservation_data = pd.read_csv(data_file, header=None, names=['conservation_rate'])
    window_size = 61
    smoothed_rates = conservation_data['conservation_rate'].rolling(window=window_size, center=True).mean()
    
    # Read variable regions data
    variable_regions = []
    with open(regions_file, 'r') as file:
        for line in file:
            start, end = map(int, line.strip().split('\t'))
            variable_regions.append((start, end))
    
    # Create a figure with specified size
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Plot the smoothed conservation rates curve
    ax.plot(smoothed_rates, color='black', linewidth=0.8, label='Sliding Window Average')
    
    # Highlight variable regions
    for start, end in variable_regions:
        ax.hlines(y=0.5, xmin=start, xmax=end, color='red', linewidth=2,
                  label='Variable Region' if start == variable_regions[0][0] else "")
    
    # Set labels and legend
    ax.set_xlabel('Gene Position')
    ax.set_ylabel('Conservation Rate')
    ax.legend()
    
    # Set y-axis limits
    ax.set_ylim(0.2, 1.0)
    
    # Save the plot to a file
    fig.savefig(output_file, bbox_inches='tight')
    plt.close(fig)

def main():
    # Specify input and output file paths
    conservation_data_file = 'solution-problem-1.txt'
    variable_regions_file = 'solution-problem-3.txt'
    output_plot_file = 'solution-problem-4.pdf'
    
    # Plot conservation rates and highlight variable regions
    plot_conservation_rates(conservation_data_file, variable_regions_file, output_plot_file)
    
    print(f"Conservation rates plot saved to: {output_plot_file}")

if __name__ == '__main__':
    main()