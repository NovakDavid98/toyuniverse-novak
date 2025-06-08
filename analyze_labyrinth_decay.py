import numpy as np
import matplotlib.pyplot as plt
import os
import re

def get_time_from_filename(filename):
    """Extracts simulation time from snapshot filename."""
    match = re.search(r'_t(\d+\.\d+)\.npy$', filename)
    if match:
        return float(match.group(1))
    return None

def analyze_spatial_variance():
    """
    Calculates and plots the spatial variance of Phi_field over time for labyrinth snapshots.
    """
    base_path = "/teamspace/studios/this_studio/toyuniverse-novak/"
    snapshots_dir = os.path.join(base_path, "run_labyrinth_snapshots_v2", "phi_snapshots_npy")
    output_plot_filename = os.path.join(base_path, "labyrinth_variance_decay.png")

    if not os.path.isdir(snapshots_dir):
        print(f"Error: Snapshot directory not found at {snapshots_dir}")
        return

    snapshot_files = sorted(
        [f for f in os.listdir(snapshots_dir) if f.startswith('phi_snapshot_step') and f.endswith('.npy')],
        key=lambda f: get_time_from_filename(f) if get_time_from_filename(f) is not None else float('inf')
    )

    times = []
    variances = []
    max_valid_time = -1.0

    print(f"Analyzing snapshots in: {snapshots_dir}")

    for i, filename in enumerate(snapshot_files):
        filepath = os.path.join(snapshots_dir, filename)
        sim_time = get_time_from_filename(filename)

        if sim_time is None:
            print(f"Could not extract time from {filename}, skipping.")
            continue
        
        try:
            data = np.load(filepath)
            if np.isnan(data).any():
                print(f"NaN values found in {filename} (t={sim_time:.2f}). Stopping variance calculation at this point.")
                break 
            
            spatial_variance = np.var(data)
            times.append(sim_time)
            variances.append(spatial_variance)
            max_valid_time = sim_time
            if i % 10 == 0: # Print progress every 10 files
                 print(f"Processed {filename} (t={sim_time:.2f}), variance: {spatial_variance:.4e}")

        except Exception as e:
            print(f"Error processing {filename}: {e}")
            # Optionally, decide if to break or continue on other errors
            break

    if not times:
        print("No valid snapshots found or processed to plot variance.")
        return

    plt.figure(figsize=(10, 6))
    plt.plot(times, variances, marker='o', linestyle='-')
    plt.xlabel("Time")
    plt.ylabel("Spatial Variance of Phi_field")
    plt.title("Labyrinth Pattern Decay: Spatial Variance over Time")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_plot_filename)
    print(f"\nVariance decay plot saved to: {output_plot_filename}")
    print(f"Analysis included snapshots up to t={max_valid_time:.2f}.")

if __name__ == '__main__':
    analyze_spatial_variance()
