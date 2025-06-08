import numpy as np
import matplotlib.pyplot as plt
import os

def create_qualitative_comparison_plot():
    base_path = "/teamspace/studios/this_studio/toyuniverse-novak/"
    spot_run_dir = os.path.join(base_path, "run_spot_snapshots_v2", "phi_snapshots_npy")
    labyrinth_run_dir = os.path.join(base_path, "run_labyrinth_snapshots_v2", "phi_snapshots_npy")

    target_times_info = [
        {"t": 20.00, "step": "0002000"},
        {"t": 50.00, "step": "0005000"},
        {"t": 100.00, "step": "0010000"},
        {"t": 200.00, "step": "0020000"},
        {"t": 295.00, "step": "0029500"} # Closest to t=300
    ]

    num_times = len(target_times_info)
    fig, axes = plt.subplots(num_times, 2, figsize=(10, 5 * num_times))
    fig.suptitle('Qualitative Spatial Comparison: Spot vs. Labyrinth Patterns', fontsize=16)

    for i, time_info in enumerate(target_times_info):
        t_val = time_info["t"]
        step_str = time_info["step"]

        spot_snapshot_file = f"phi_snapshot_step{step_str}_t{t_val:.2f}.npy"
        labyrinth_snapshot_file = f"phi_snapshot_step{step_str}_t{t_val:.2f}.npy"

        spot_path = os.path.join(spot_run_dir, spot_snapshot_file)
        labyrinth_path = os.path.join(labyrinth_run_dir, labyrinth_snapshot_file)

        # Load and plot Spot Pattern
        try:
            spot_data = np.load(spot_path)
            ax_spot = axes[i, 0]
            im_spot = ax_spot.imshow(spot_data, cmap='hot', interpolation='nearest')
            ax_spot.set_title(f'Spot Pattern (t={t_val:.2f})')
            ax_spot.set_xticks([])
            ax_spot.set_yticks([])
            # fig.colorbar(im_spot, ax=ax_spot, fraction=0.046, pad=0.04) # Optional: add colorbar
        except FileNotFoundError:
            print(f"Error: Spot snapshot not found at {spot_path}")
            ax_spot = axes[i, 0]
            ax_spot.text(0.5, 0.5, 'Spot data not found', ha='center', va='center')
            ax_spot.set_title(f'Spot Pattern (t={t_val:.2f}) - Error')
            ax_spot.set_xticks([])
            ax_spot.set_yticks([])
        except Exception as e:
            print(f"Error loading spot data {spot_path}: {e}")
            ax_spot = axes[i, 0]
            ax_spot.text(0.5, 0.5, f'Error loading spot data:\n{e}', ha='center', va='center', fontsize=8)
            ax_spot.set_title(f'Spot Pattern (t={t_val:.2f}) - Error')
            ax_spot.set_xticks([])
            ax_spot.set_yticks([])

        # Load and plot Labyrinth Pattern
        try:
            labyrinth_data = np.load(labyrinth_path)
            ax_labyrinth = axes[i, 1]
            labyrinth_mean = np.mean(labyrinth_data)
            labyrinth_std = np.std(labyrinth_data)
            vmin_lab = labyrinth_mean - 2.5 * labyrinth_std
            vmax_lab = labyrinth_mean + 2.5 * labyrinth_std
            
            # Ensure vmin is not greater than vmax, can happen if std is zero or very small
            if vmin_lab >= vmax_lab:
                vmin_lab = np.min(labyrinth_data) # Fallback to min/max if std is problematic
                vmax_lab = np.max(labyrinth_data)
                if vmin_lab == vmax_lab: # If still equal (e.g. constant array), add small epsilon for imshow
                    vmin_lab -= 0.001
                    vmax_lab += 0.001
            
            im_labyrinth = ax_labyrinth.imshow(labyrinth_data, cmap='hot', interpolation='nearest', vmin=vmin_lab, vmax=vmax_lab)
            ax_labyrinth.set_title(f'Labyrinth Pattern (t={t_val:.2f})\n(vmin={vmin_lab:.2e}, vmax={vmax_lab:.2e})')
            ax_labyrinth.set_xticks([])
            ax_labyrinth.set_yticks([])
            # fig.colorbar(im_labyrinth, ax=ax_labyrinth, fraction=0.046, pad=0.04) # Optional: add colorbar
        except FileNotFoundError:
            print(f"Error: Labyrinth snapshot not found at {labyrinth_path}")
            ax_labyrinth = axes[i, 1]
            ax_labyrinth.text(0.5, 0.5, 'Labyrinth data not found', ha='center', va='center')
            ax_labyrinth.set_title(f'Labyrinth Pattern (t={t_val:.2f}) - Error')
            ax_labyrinth.set_xticks([])
            ax_labyrinth.set_yticks([])
        except Exception as e:
            print(f"Error loading labyrinth data {labyrinth_path}: {e}")
            ax_labyrinth = axes[i, 1]
            ax_labyrinth.text(0.5, 0.5, f'Error loading labyrinth data:\n{e}', ha='center', va='center', fontsize=8)
            ax_labyrinth.set_title(f'Labyrinth Pattern (t={t_val:.2f}) - Error')
            ax_labyrinth.set_xticks([])
            ax_labyrinth.set_yticks([])

    plt.tight_layout(rect=[0, 0, 1, 0.96]) # Adjust layout to make space for suptitle
    output_filename = os.path.join(base_path, 'qualitative_comparison_snapshots.png')
    plt.savefig(output_filename, dpi=300)
    print(f"Comparison plot saved to {output_filename}")

if __name__ == '__main__':
    create_qualitative_comparison_plot()
