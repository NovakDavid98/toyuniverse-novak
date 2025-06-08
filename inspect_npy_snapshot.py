import numpy as np
import os
import argparse

def inspect_snapshot(npy_file_path):
    """
    Loads a .npy file, prints its basic statistics, and a small slice.
    """
    if not os.path.exists(npy_file_path):
        print(f"Error: File not found at {npy_file_path}")
        return

    try:
        data = np.load(npy_file_path)
        print(f"Successfully loaded: {npy_file_path}")
        print(f"Data type: {data.dtype}")
        print(f"Shape: {data.shape}")

        if data.size == 0:
            print("The array is empty.")
            return

        print(f"\n--- Basic Statistics ---")
        print(f"Minimum value: {np.min(data):.6e}")
        print(f"Maximum value: {np.max(data):.6e}")
        print(f"Mean value:    {np.mean(data):.6e}")
        print(f"Std deviation: {np.std(data):.6e}")

        print(f"\n--- NaN/Inf Check ---")
        num_nan = np.isnan(data).sum()
        num_inf = np.isinf(data).sum()
        print(f"Number of NaN values: {num_nan}")
        print(f"Number of Inf values: {num_inf}")

        print(f"\n--- Array Slice (top-left 5x5 or smaller) ---")
        slice_rows = min(5, data.shape[0])
        slice_cols = min(5, data.shape[1])
        print(data[:slice_rows, :slice_cols])

    except Exception as e:
        print(f"An error occurred while processing {npy_file_path}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Inspect a .npy snapshot file.")
    parser.add_argument("npy_file_path", type=str, help="Path to the .npy file to inspect.")
    args = parser.parse_args()
    inspect_snapshot(args.npy_file_path)
