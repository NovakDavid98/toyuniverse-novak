import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os

# Define output directory for plots
PLOTS_OUTPUT_DIR = "results1"
os.makedirs(PLOTS_OUTPUT_DIR, exist_ok=True)

def load_data(regime_suffix):
    """Loads time-series data for a given regime (e.g., '_spot' or '_labyrinth')."""
    print(f"Loading data for regime: {regime_suffix}...")
    data = {}
    try:
        data['t'] = np.load(f't_values{regime_suffix}.npy')
        data['psi'] = np.load(f'psi_id14_values{regime_suffix}.npy')
        data['i'] = np.load(f'i_id14_values{regime_suffix}.npy')
        data['r'] = np.load(f'r_id14_values{regime_suffix}.npy')
        data['avg_phi'] = np.load(f'avg_phi_field_values{regime_suffix}.npy')
        print(f"Successfully loaded all .npy files for {regime_suffix}.")
    except FileNotFoundError as e:
        print(f"Error loading data for {regime_suffix}: {e}")
        return None
    return data

def plot_temporal_comparison(spot_data, labyrinth_data):
    """Plots temporal evolution of ID14 variables and avg_Phi_field for both regimes."""
    print("Plotting temporal comparison...")
    if not spot_data or not labyrinth_data:
        print("Missing data for one or both regimes. Skipping temporal comparison plot.")
        return

    variables_to_plot = {
        'psi': ('Psi_ID14', 'Psi_ID14'),
        'i': ('I_ID14', 'I_ID14'),
        'r': ('R_ID14', 'R_ID14'),
        'avg_phi': ('Average Phi Field', 'avg_Phi_field')
    }

    fig, axs = plt.subplots(len(variables_to_plot), 1, figsize=(12, 4 * len(variables_to_plot)), sharex=True)
    fig.suptitle('Temporal Evolution Comparison: Spot vs. Labyrinthine', fontsize=16)

    for idx, (var_key, (title, ylabel)) in enumerate(variables_to_plot.items()):
        ax = axs[idx]
        ax.plot(spot_data['t'], spot_data[var_key], label=f'Spot ({var_key})', color='blue')
        ax.plot(labyrinth_data['t'], labyrinth_data[var_key], label=f'Labyrinth ({var_key})', color='red')
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.legend()
        ax.grid(True)

    axs[-1].set_xlabel('Time')
    plt.tight_layout(rect=[0, 0, 1, 0.96]) # Adjust layout to make space for suptitle
    plot_filename = os.path.join(PLOTS_OUTPUT_DIR, 'temporal_comparison.png')
    plt.savefig(plot_filename)
    print(f"Saved temporal comparison plot to {plot_filename}")
    plt.close(fig)

def plot_phase_space(data, regime_suffix, var1_key, var2_key):
    """Plots a phase space diagram for two variables of a given regime."""
    print(f"Plotting phase space for {regime_suffix}: {var1_key} vs {var2_key}...")
    if not data:
        print(f"Missing data for {regime_suffix}. Skipping phase space plot.")
        return

    plt.figure(figsize=(8, 6))
    plt.plot(data[var1_key], data[var2_key])
    plt.title(f'Phase Space: {var1_key.upper()} vs {var2_key.upper()} ({regime_suffix.strip("_")})')
    plt.xlabel(var1_key.upper())
    plt.ylabel(var2_key.upper())
    plt.grid(True)
    plot_filename = os.path.join(PLOTS_OUTPUT_DIR, f'phase_space_{var1_key}_vs_{var2_key}{regime_suffix}.png')
    plt.savefig(plot_filename)
    print(f"Saved phase space plot to {plot_filename}")
    plt.close()

def plot_side_by_side_snapshots(time_points_for_snapshots):
    """Plots side-by-side snapshots for spot and labyrinthine regimes at given time points."""
    print("Plotting side-by-side snapshots...")
    num_times = len(time_points_for_snapshots)
    fig, axs = plt.subplots(num_times, 2, figsize=(10, 5 * num_times))
    fig.suptitle('Comparative Snapshots: Spot (left) vs. Labyrinthine (right)', fontsize=16)

    for i, t_val in enumerate(time_points_for_snapshots):
        # Frame numbers are 1-indexed and correspond to t_val if t_max=300 and 300 frames saved
        frame_num = int(t_val) # Assuming t_val directly maps to frame number for these specific times
        frame_filename_template = f"frame_{frame_num:04d}.png"

        spot_frame_path = os.path.join('seriesofimages_spot', frame_filename_template)
        labyrinth_frame_path = os.path.join('seriesofimages_labyrinth', frame_filename_template)

        # Spot regime plot
        ax_spot = axs[i, 0] if num_times > 1 else axs[0]
        try:
            img_spot = mpimg.imread(spot_frame_path)
            ax_spot.imshow(img_spot, cmap=ANIMATION_CMAP_ANALYSIS if 'ANIMATION_CMAP_ANALYSIS' in globals() else 'hot')
            ax_spot.set_title(f'Spot Regime: t = {t_val}')
            ax_spot.axis('off')
        except FileNotFoundError:
            ax_spot.text(0.5, 0.5, f'Spot frame\n{frame_filename_template}\nnot found', ha='center', va='center')
            ax_spot.set_title(f'Spot Regime: t = {t_val}')
            ax_spot.axis('off')

        # Labyrinthine regime plot
        ax_labyrinth = axs[i, 1] if num_times > 1 else axs[1]
        try:
            img_labyrinth = mpimg.imread(labyrinth_frame_path)
            ax_labyrinth.imshow(img_labyrinth, cmap=ANIMATION_CMAP_ANALYSIS if 'ANIMATION_CMAP_ANALYSIS' in globals() else 'hot')
            ax_labyrinth.set_title(f'Labyrinthine Regime: t = {t_val}')
            ax_labyrinth.axis('off')
        except FileNotFoundError:
            ax_labyrinth.text(0.5, 0.5, f'Labyrinth frame\n{frame_filename_template}\nnot found', ha='center', va='center')
            ax_labyrinth.set_title(f'Labyrinthine Regime: t = {t_val}')
            ax_labyrinth.axis('off')

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plot_filename = os.path.join(PLOTS_OUTPUT_DIR, 'comparative_snapshots.png')
    plt.savefig(plot_filename)
    print(f"Saved comparative snapshots to {plot_filename}")
    plt.close(fig)

# Define a colormap consistent with simulation if possible, otherwise default to 'hot'
# This should ideally be read from the simulation script or set consistently.
ANIMATION_CMAP_ANALYSIS = 'hot' 

if __name__ == '__main__':
    print("Starting analysis script...")
    # Load data for both regimes
    spot_data = load_data('_spot')
    labyrinth_data = load_data('_labyrinth')

    # Generate temporal comparison plot
    plot_temporal_comparison(spot_data, labyrinth_data)

    # Generate phase space plots for each regime
    phase_space_pairs = [('psi', 'i'), ('psi', 'r'), ('i', 'r')]

    if spot_data:
        for var1, var2 in phase_space_pairs:
            plot_phase_space(spot_data, '_spot', var1, var2)
    
    if labyrinth_data:
        for var1, var2 in phase_space_pairs:
            plot_phase_space(labyrinth_data, '_labyrinth', var1, var2)

    # Generate side-by-side comparative snapshots
    # Time points suggested by scientist: t=20, 50, 100, 200, 300
    # Assuming 300 frames were saved for t_max=300, so frame number = t_value
    snapshot_time_points = [20, 50, 100, 200, 300]
    plot_side_by_side_snapshots(snapshot_time_points)

    print("Analysis script finished.")
