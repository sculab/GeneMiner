import numpy as np
import argparse

def estimate_k(p, v, l):
    """
    Estimate the most likely value of k given p, v, and l without using scipy.

    Parameters:
    p (float): A real number between 0 and 1.
    v (float): A real number between 0 and 1.
    l (int): An integer between 50 and 500.

    Returns:
    float: The estimated value of k, rounded to the nearest integer.
    """
    # Define the function to minimize
    def equation_diff(k):
        return abs(p - (1 - (1 - (1-v)**k)**(l - k + 1)))

    # Initialize minimum difference and corresponding k value
    min_diff = float('inf')
    min_k = 9

    # Brute-force search in the range of k
    for k in range(9, l):
        diff = equation_diff(k)
        if diff < min_diff:
            min_diff = diff
            min_k = k

    # Return the k value that minimizes the difference
    return min_k

# Set up argument parsing
parser = argparse.ArgumentParser(description="Estimate the most likely value of k given p, v, and l. For instance: python estimate_k.py 0.99 0.1 150")
parser.add_argument("p", type=float, help="Expected recall rate of reads. A real number between 0 and 1.")
parser.add_argument("v", type=float, help="Mutation rate of the reference sequence relative to the target sequence. A real number between 0 and 1.")
parser.add_argument("l", type=int, help="Reads length. An integer between 50 and 500.")

# Parse arguments
args = parser.parse_args()

# Call the function with provided arguments
k_estimated = estimate_k(args.p, args.v, args.l)
print("Estimated k:", k_estimated)
