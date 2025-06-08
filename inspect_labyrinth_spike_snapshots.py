import numpy as np
import matplotlib.pyplot as plt
import os
import re

def get_step_from_time(time_val):
    """Converts time to step string, assuming dt=0.01 and snapshot_factor=1 (for step in filename)."""
    # This needs to align with how steps were actually saved if not directly from time.
    # Assuming step = time * 100 for these files based on common pattern.
    return int(time_val * 100)

def plot_specific_labyrinth_snapshots():
    base_path = "/teamspace/studios/this_studio/toyuniverse-novak/"
    snapshots_dir = os.path.join(base_path, "run_labyrinth_snapshots_v2", "phi_snapshots_npy")
    output_plot_filename = os.path.join(base_path, "labyrinth_spike_snapshots.png")

    target_times = [100.00, 105.00, 110.00, 115.00, 120.00, 125.00]
    
    num_times = len(target_times)
    # Adjust layout based on number of times, e.g., 2 rows, 3 cols for 6 images
    ncols = 3
    nrows = (num_times + ncols - 1) // ncols 
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 5, nrows * 5), squeeze=False)
    fig.suptitle('Labyrinth Snapshots During High Variance Period', fontsize=16)

    plot_idx = 0
    for t_val in target_times:
        step_int = get_step_from_time(t_val) # Correctly calculate step based on your simulation's dt and snapshot saving frequency
        # We need to find the exact filename. The step might not be simply t*100 if snapshot_storage_skip_factor was used differently.
        # Let's search for the file by time string, as that's more robust from previous scripts.
        expected_filename_pattern = f"phi_snapshot_step.*?_t{t_val:.2f}.npy"
        
        found_file = None
        for fname in os.listdir(snapshots_dir):
            if re.fullmatch(expected_filename_pattern, fname):
                found_file = fname
                break
        
        ax = axes[plot_idx // ncols, plot_idx % ncols]

        if found_file:
            filepath = os.path.join(snapshots_dir, found_file)
            try:
                data = np.load(filepath)
                if np.isnan(data).any():
                    ax.text(0.5, 0.5, 'Data contains NaNs',
                            horizontalalignment='center', verticalalignment='center',
                            transform=ax.transAxes, color='red')
                    ax.set_title(f'Labyrinth Pattern (t={t_val:.2f}) - NaNs')
                else:
                    data_mean = np.mean(data)
                    data_std = np.std(data)
                    vmin = data_mean - 2.5 * data_std
                    vmax = data_mean + 2.5 * data_std
                    if vmin >= vmax: # Fallback
                        vmin = np.min(data)
                        vmax = np.max(data)
                        if vmin == vmax:
                            vmin -= 0.001
                            vmax += 0.001
                    
                    im = ax.imshow(data, cmap='hot', interpolation='nearest', vmin=vmin, vmax=vmax)
                    ax.set_title(f'Labyrinth (t={t_val:.2f})\n(vmin={vmin:.2e}, vmax={vmax:.2e})')
                    # fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04) # Optional

            except FileNotFoundError:
                ax.text(0.5, 0.5, 'File not found',
                        horizontalalignment='center', verticalalignment='center',
                        transform=ax.transAxes)
                ax.set_title(f'Labyrinth Pattern (t={t_val:.2f}) - Not Found')
            except Exception as e:
                ax.text(0.5, 0.5, f'Error: {e}',
                        horizontalalignment='center', verticalalignment='center',
                        transform=ax.transAxes, color='red', fontsize=8)
                ax.set_title(f'Labyrinth Pattern (t={t_val:.2f}) - Error')
        else:
            ax.text(0.5, 0.5, 'File not found by pattern',
                    horizontalalignment='center', verticalalignment='center',
                    transform=ax.transAxes)
            ax.set_title(f'Labyrinth Pattern (t={t_val:.2f}) - Not Found')

        ax.set_xticks([])
        ax.set_yticks([])
        plot_idx += 1

    # Hide any unused subplots if num_times is not a perfect multiple for grid
    for i in range(plot_idx, nrows * ncols):
        fig.delaxes(axes[i // ncols, i % ncols])

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output_plot_filename, dpi=300)
    print(f"Specific snapshots plot saved to {output_plot_filename}")

if __name__ == '__main__':
    plot_specific_labyrinth_snapshots()
