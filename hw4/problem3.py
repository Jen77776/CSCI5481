import numpy as np
import pandas as pd

def identify_variable_regions(data_file, threshold=0.75, min_length=30):
    """
    Identify variable regions in the conservation rates data.

    Parameters:
    data_file (str): Path to the file containing conservation rates data.
    threshold (float): Threshold value for determining variable regions (default: 0.75).
    min_length (int): Minimum length of a region to be considered variable (default: 30).

    Returns:
    list: A list of tuples representing variable regions, where each tuple contains the start and end positions.
    """
    # Read conservation rates data from file
    conservation_data = pd.read_csv(data_file, header=None, names=['conservation_rate'])
    
    # Apply rolling window smoothing to the conservation rates
    window_size = 61
    smoothed_rates = conservation_data['conservation_rate'].rolling(window=window_size, center=True).mean()
    
    # Find positions with conservation rates below the threshold
    variable_positions = []
    region_start = None
    
    # Iterate over the smoothed conservation rates to identify variable regions
    for position, rate in enumerate(smoothed_rates):
        if pd.isna(rate):
            # Skip NaN values at the beginning and end of the smoothed data
            continue
        
        if region_start is None and rate < threshold:
            # Start a new variable region
            region_start = position
        elif region_start is not None and rate >= threshold:
            # End the current variable region
            if position - region_start >= min_length:
                variable_positions.append((region_start, position))
            region_start = None
    
    # Handle the last variable region if it extends to the end of the data
    if region_start is not None and len(smoothed_rates) - region_start >= min_length:
        variable_positions.append((region_start, len(smoothed_rates) - 1))
    
    return variable_positions

def main():
    # Specify the input data file
    data_file = 'solution-problem-1.txt'
    
    # Identify variable regions
    variable_regions = identify_variable_regions(data_file)
    
    # Save the results to a file
    output_file = 'solution-problem-3.txt'
    with open(output_file, 'w') as file:
        for start, end in variable_regions:
            file.write(f"{start}\t{end}\n")
    
    # Print the identified variable regions
    num_regions = len(variable_regions)
    print(f"Identified {num_regions} variable regions:")
    for i, (start, end) in enumerate(variable_regions, 1):
        region_length = end - start + 1
        print(f"Region {i}: {start}-{end} (Length: {region_length})")

if __name__ == '__main__':
    main()