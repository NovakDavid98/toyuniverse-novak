import numpy as np
import os

def inspect_simulation_times(npy_file_path, target_times):
    """
    Loads a .npy file containing time values for simulation frames
    and prints information about it, including the closest frames
    for a list of target times.
    """
    if not os.path.exists(npy_file_path):
        print(f"Error: File not found at {npy_file_path}")
        return

    try:
        t_values = np.load(npy_file_path)
        print(f"Successfully loaded: {npy_file_path}")
        print(f"Data type: {t_values.dtype}")
        print(f"Shape: {t_values.shape}")
        num_entries = len(t_values)
        print(f"Number of time entries: {num_entries}")

        if num_entries == 0:
            print("The file is empty.")
            return

        print(f"First 5 time values: {t_values[:5]}")
        if num_entries > 5:
            print(f"Last 5 time values: {t_values[-5:]}")
        else:
            print(f"All time values: {t_values}")

        min_time = np.min(t_values)
        max_time = np.max(t_values)
        print(f"Minimum time: {min_time:.2f}")
        print(f"Maximum time: {max_time:.2f}")

        print("\n--- Frames closest to target times ---")
        # Assuming frame filenames are 0-indexed, e.g., frame_0000.png, frame_0001.png, ...
        # The index from t_values directly corresponds to the number in the filename.
        for target_t in target_times:
            if target_t < min_time and not np.isclose(target_t, min_time):
                 print(f"Target time t={target_t:.2f} is below the range of available times (starting from {min_time:.2f}).")
                 continue
            if target_t > max_time and not np.isclose(target_t, max_time):
                 print(f"Target time t={target_t:.2f} is above the range of available times (up to {max_time:.2f}).")
                 continue
            closest_idx = np.argmin(np.abs(t_values - target_t))
            actual_t = t_values[closest_idx]
            print(f"Target t={target_t:.2f}: Closest frame index is {closest_idx} (corresponds to frame_{closest_idx:04d}.png), actual t={actual_t:.2f}")

    except Exception as e:
        print(f"An error occurred while processing {npy_file_path}: {e}")

if __name__ == "__main__":
    file_path = "/teamspace/studios/this_studio/toyuniverse-novak/simulationdatajun808pm/t_values.npy"
    targets = [20.0, 50.0, 100.0, 200.0, 300.0]
    inspect_simulation_times(file_path, targets)
