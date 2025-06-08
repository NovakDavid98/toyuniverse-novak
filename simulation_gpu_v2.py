#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GPU-accelerated simulation of the Toy Universe model (Version 2).

This script simulates a PDE-ODE model, with optional GPU acceleration using CuPy.
It includes features for:
- Comprehensive command-line argument parsing.
- GPU/CPU agnostic array operations.
- Live animation of the PDE field.
- Saving individual animation frames as PNGs.
- Saving final state and time-series data as NPY files.
- Profiling (cProfile) and memory usage tracking (tracemalloc).
"""

import numpy as np
import os
import argparse
import time
import sys
import gc
import tracemalloc # For memory profiling

# Plotting and Animation
import matplotlib
matplotlib.use('Agg') # Use Agg backend for non-interactive plotting
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Profiling
import cProfile
import pstats
import io

# Attempt to import CuPy for GPU acceleration
cupy_available = False
try:
    import cupy
    cupy_available = True
    print("INFO: CuPy found, GPU acceleration is available.")
except ImportError:
    print("INFO: CuPy not found, GPU acceleration is NOT available. Falling back to NumPy.")

def parse_arguments():
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description="Toy Universe PDE-ODE Simulation (Version 2)")

    # GPU and Backend
    parser.add_argument('--use_gpu', action='store_true', help="Enable GPU acceleration with CuPy.")

    # Simulation Parameters
    parser.add_argument('--t_max', type=float, default=100.0, help="Maximum simulation time.")
    parser.add_argument('--dt', type=float, default=0.01, help="Time step.")
    parser.add_argument('--Nx', type=int, default=128, help="Number of grid points in x.")
    parser.add_argument('--Ny', type=int, default=128, help="Number of grid points in y.")
    parser.add_argument('--Lx', type=float, default=10.0, help="Length of domain in x.")
    parser.add_argument('--Ly', type=float, default=10.0, help="Length of domain in y.")

    # PDE Coefficients (Phi_field)
    parser.add_argument('--gamma_Phi_PDE', type=float, default=0.01, help='gamma_Phi_PDE parameter (UNUSED in current classic PDE model, D_Phi_PDE is for diffusion)')
    parser.add_argument('--alpha_Phi_PDE', type=float, default=2.1, help='PDE: alpha_Phi_PDE, coefficient for Psi_c driven growth term1 (default 2.1)')
    parser.add_argument('--Phi0_saturation_alpha', type=float, default=0.5, help='PDE: Phi0_saturation_alpha, saturation scale in tanh for alpha_Phi_PDE term (default 0.5)')
    parser.add_argument('--epsilon_Phi_PDE', type=float, default=2.1, help='PDE: epsilon_Phi_PDE, coefficient for Psi_d/rho_elec driven feedback term4 (default 2.1)')
    parser.add_argument('--Phi0_saturation_epsilon_PDE', type=float, default=0.5, help='PDE: Phi0_saturation_epsilon_PDE, saturation scale in tanh for epsilon_Phi_PDE term (default 0.5)')
    parser.add_argument('--D_Phi_PDE', type=float, default=0.01, help="PDE: D_Phi_PDE, diffusion coefficient for Phi_field (default 0.01 for spots, try 0.001 for labyrinths)")
    parser.add_argument('--delta_Phi_PDE', type=float, default=1.0, help="PDE: delta_Phi_PDE (UNUSED in current classic PDE model)")
    parser.add_argument('--beta_Phi_PDE', type=float, default=0.0001, help="PDE: beta_Phi_PDE, coefficient for R_rad_curr damping term2 (default 0.0001 for patterns)")

    # New parameters for intermediate field calculations in classic PDE model
    parser.add_argument('--Psi0_scaling_A', type=float, default=1.0, help='PDE: Psi0_scaling_A for Psi_c_s = tanh(Psi_ID14 / Psi0_scaling_A)')
    parser.add_argument('--Psi0_scaling_B', type=float, default=1.0, help='PDE: Psi0_scaling_B for Psi_d_s = tanh(Psi_ID14 * Psi0_scaling_B)')
    parser.add_argument('--k_v_p', type=float, default=1.0, help='PDE: k_v_p for V_field calculation in R_rad_curr')
    parser.add_argument('--rho_eps_p', type=float, default=1e-6, help='PDE: rho_eps_p, small constant for V_field denominator in R_rad_curr')
    parser.add_argument('--kappa_R_p', type=float, default=1.0, help='PDE: kappa_R_p for R_rad_curr calculation')
    parser.add_argument('--k_R_int_p', type=float, default=1.0, help='PDE: k_R_int_p for R_rad_curr calculation')

    # ODE Coefficients (Psi_ID14, I_ID14, R_ID14)
    parser.add_argument('--alpha_Psi_ODE', type=float, default=1.0, help="ODE: alpha_Psi coefficient")
    parser.add_argument('--beta_Psi_ODE', type=float, default=1.0, help="ODE: beta_Psi coefficient")
    parser.add_argument('--gamma_Psi_ODE', type=float, default=1.0, help="ODE: gamma_Psi coefficient")
    parser.add_argument('--delta_Psi_ODE', type=float, default=1.0, help="ODE: delta_Psi coefficient")
    parser.add_argument('--epsilon_Psi_ODE', type=float, default=1.0, help="ODE: epsilon_Psi coefficient")
    parser.add_argument('--zeta_Psi_ODE', type=float, default=1.0, help="ODE: zeta_Psi coefficient")
    parser.add_argument('--eta_Psi_ODE', type=float, default=1.0, help="ODE: eta_Psi coefficient")
    parser.add_argument('--theta_Psi_ODE', type=float, default=1.0, help="ODE: theta_Psi coefficient")
    parser.add_argument('--kappa_Psi_ODE', type=float, default=1.0, help="ODE: kappa_Psi coefficient")
    parser.add_argument('--lambda_Psi_ODE', type=float, default=1.0, help="ODE: lambda_Psi coefficient")
    parser.add_argument('--mu_Psi_ODE', type=float, default=1.0, help="ODE: mu_Psi coefficient")
    parser.add_argument('--nu_Psi_ODE', type=float, default=1.0, help="ODE: nu_Psi coefficient")
    parser.add_argument('--xi_Psi_ODE', type=float, default=1.0, help="ODE: xi_Psi coefficient")
    parser.add_argument('--rho_Psi_ODE', type=float, default=1.0, help="ODE: rho_Psi coefficient")
    parser.add_argument('--sigma_Psi_ODE', type=float, default=1.0, help="ODE: sigma_Psi coefficient")

    # Initial Conditions
    parser.add_argument('--Phi0_noise_amplitude_INIT', type=float, default=0.01, help="Initial noise amplitude for Phi_field.")
    parser.add_argument('--Psi0_INIT', type=float, default=0.1, help="Initial value for Psi_ID14.")
    parser.add_argument('--I0_INIT', type=float, default=0.1, help="Initial value for I_ID14.")
    parser.add_argument('--R0_INIT', type=float, default=0.1, help="Initial value for R_ID14.")

    # Output Control
    parser.add_argument('--storage_skip_factor', type=int, default=10, help="Store data every Nth time step.")
    parser.add_argument('--MAIN_OUTPUT_DIR', type=str, default="simulation_outputs_v2", help="Main directory for simulation outputs.")
    parser.add_argument('--ANIMATION_CMAP', type=str, default='hot', help="Colormap for the animation GIF.")
    parser.add_argument('--FINAL_STATE_CMAP', type=str, default='viridis', help="Colormap for the final state image.")

    # NPY Snapshot Output Control
    parser.add_argument('--phi_snapshot_dir_name', type=str, default="phi_snapshots_npy", \
                        help="Directory name for 2D Phi_field NPY snapshots, relative to MAIN_OUTPUT_DIR. Empty string or snapshot_storage_skip_factor <= 0 disables snapshots.")
    parser.add_argument('--snapshot_storage_skip_factor', type=int, default=500, \
                        help="Store 2D Phi_field NPY snapshot every Nth time step. If 0 or negative, snapshots are disabled. (e.g., 500 for ~60 snapshots if t_max=300, dt=0.01)")

    return parser.parse_args()

# --- Global Variables & Helper Functions ---
XP_MODULE = np # Default to NumPy
SERIES_OF_IMAGES_DIR = ""

def print_memory_usage(tag=""):
    """Prints current and peak memory usage using tracemalloc."""
    gc.collect() # Run garbage collection first
    try:
        current, peak = tracemalloc.get_traced_memory()
        print(f"MEMORY_USAGE [{tag}]: Current memory usage: {current / 10**6:.2f}MB; Peak: {peak / 10**6:.2f}MB")
    except Exception as e:
        print(f"MEMORY_USAGE [{tag}]: Could not get tracemalloc stats (is it started?): {e}")
    sys.stdout.flush()

# --- Helper Functions for Simulation ---
def get_laplacian(field, dx, dy, xp):
    """Calculates the 2D Laplacian of a field using finite differences and periodic boundary conditions."""
    # field_xp = xp.asarray(field) # Ensure it's an xp array, though it should be already
    # Using np.roll for periodic boundary conditions, works for both numpy and cupy
    term_x = (xp.roll(field, -1, axis=0) - 2 * field + xp.roll(field, 1, axis=0)) / (dx**2)
    term_y = (xp.roll(field, -1, axis=1) - 2 * field + xp.roll(field, 1, axis=1)) / (dy**2)
    return term_x + term_y

def pde_model(Phi_field, Psi_c_s, Psi_d_s, rho_elec_field, R_rad_field, args, dx, dy, xp):
    """Calculates dPhi_dt based on the classic pattern-forming PDE model.
    This structure is based on the user's reconstruction of the previously successful model.
    Inputs:
        Phi_field: Current 2D Phi_field array.
        Psi_c_s: Scalar modulating factor from Psi_ID14 for term1.
        Psi_d_s: Scalar modulating factor from Psi_ID14 for term4.
        rho_elec_field: 2D electric charge density field.
        R_rad_field: 2D radiative radius/flux field.
        args: Namespace containing parameters.
        dx, dy: Grid spacings.
        xp: NumPy or CuPy module.
    """
    laplacian_Phi = get_laplacian(Phi_field, dx, dy, xp)

    # Term 1: Psi_c driven growth, with Phi self-saturation
    # effective_Phi_for_alpha_term = Phi0_saturation_alpha * tanh(Phi_curr / Phi0_saturation_alpha)
    effective_Phi_for_alpha_term = args.Phi0_saturation_alpha * xp.tanh(Phi_field / (args.Phi0_saturation_alpha + 1e-9))
    term1_growth_psi_c = args.alpha_Phi_PDE * effective_Phi_for_alpha_term * Psi_c_s

    # Term 2: Damping related to R_radiative_field
    # term2_damping_R_rad = -beta_Phi_PDE * R_rad_curr
    term2_damping_R_rad = -args.beta_Phi_PDE * R_rad_field
    
    # Term 3: Diffusion
    # term3_diffusion = D_Phi_PDE * laplacian_Phi_curr
    term3_diffusion = args.D_Phi_PDE * laplacian_Phi

    # Term 4: Psi_d and rho_electric driven feedback, with Phi self-saturation
    # effective_Phi_for_epsilon_term = Phi0_saturation_epsilon_PDE * tanh(Phi_curr / Phi0_saturation_epsilon_PDE)
    effective_Phi_for_epsilon_term = args.Phi0_saturation_epsilon_PDE * xp.tanh(Phi_field / (args.Phi0_saturation_epsilon_PDE + 1e-9))
    term4_feedback_psi_d_rho = args.epsilon_Phi_PDE * effective_Phi_for_epsilon_term * Psi_d_s * rho_elec_field

    # Summing the terms
    dPhi_dt = term1_growth_psi_c + term2_damping_R_rad + term3_diffusion + term4_feedback_psi_d_rho
    
    return dPhi_dt

def ode_model(Psi_ID14, I_ID14, R_ID14, Phi_avg, args, xp):
    """Calculates dPsi_dt, dI_dt, dR_dt based on the ODE model equations (ID14)."""
    # Ensure inputs are xp arrays if they are scalars, for consistent operations
    # However, they should already be xp.array from initialization

    # dPsi_dt equation
    dPsi_dt = (
        args.alpha_Psi_ODE * Psi_ID14 * (1.0 - Psi_ID14 / args.beta_Psi_ODE)
        - args.gamma_Psi_ODE * Psi_ID14 * I_ID14 / (1.0 + args.delta_Psi_ODE * I_ID14**2)
        - args.epsilon_Psi_ODE * Psi_ID14 * R_ID14
        + args.zeta_Psi_ODE * Phi_avg
    )

    # dI_dt equation
    dI_dt = (
        args.eta_Psi_ODE * I_ID14 * (args.theta_Psi_ODE - I_ID14) * (I_ID14 - args.kappa_Psi_ODE)
        - args.lambda_Psi_ODE * Psi_ID14 * I_ID14
        - args.mu_Psi_ODE * I_ID14 * R_ID14
    )

    # dR_dt equation
    dR_dt = (
        args.nu_Psi_ODE * R_ID14 * (args.xi_Psi_ODE - R_ID14)
        - args.rho_Psi_ODE * Psi_ID14 * R_ID14
        - args.sigma_Psi_ODE * I_ID14 * R_ID14
    )
    return dPsi_dt, dI_dt, dR_dt

def calculate_classic_intermediate_fields(Phi_field, Psi_ID14_scalar, args, xp):
    """Calculates intermediate fields for the classic PDE model.
    
    Psi_c_s, Psi_d_s: Modulating factors derived from Psi_ID14_scalar.
    rho_elec_field: Electric charge density field (simplified version).
    R_rad_field: Radiative radius/flux field.
    """
    # Ensure Psi_ID14_scalar is an xp array for consistent operations if it's a plain float
    # (though it should be an xp.array from the ODE solver output already)
    psi_val = xp.asarray(Psi_ID14_scalar) 

    # Calculate Psi_c_s and Psi_d_s
    Psi_c_s = xp.tanh(psi_val / (args.Psi0_scaling_A + 1e-9)) # Add epsilon for safety in denominator
    Psi_d_s = xp.tanh(psi_val * args.Psi0_scaling_B)

    # Calculate rho_elec_field (simpler form based on user's note)
    rho_elec_field = 0.1 * Phi_field

    # Calculate R_rad_field
    # V_field = k_v_p * (1 / (rho_electric_field + rho_eps_p)) * (6 / np.pi) * Psi_m_s
    # (where Psi_m_s was often just Psi_ID14_current)
    pi_val = xp.pi if hasattr(xp, 'pi') else np.pi # Use xp.pi if available (CuPy), else np.pi
    
    # Add small epsilon to rho_elec_field in denominator to prevent division by zero if rho_elec_field can be zero
    # and to args.rho_eps_p itself if it's zero.
    V_field_denominator = rho_elec_field + args.rho_eps_p + 1e-9 
    V_field = args.k_v_p * (1.0 / V_field_denominator) * (6.0 / pi_val) * psi_val
    
    R_rad_field_intermediate = args.kappa_R_p - args.k_R_int_p * Phi_field * Psi_c_s * V_field
    R_rad_field = xp.maximum(R_rad_field_intermediate, 0.0) # Ensure non-negative

    return Psi_c_s, Psi_d_s, rho_elec_field, R_rad_field

def initialize_fields(args, xp):
    """Initializes Phi_field, Psi_ID14, I_ID14, R_ID14 on the chosen device (xp)."""
    print_memory_usage("Start of initialize_fields")

    # Initialize Phi_field with random noise
    # Ensure noise is generated on the correct device (CPU for np, GPU for cp)
    if xp == np:
        noise = np.random.rand(args.Nx, args.Ny) * args.Phi0_noise_amplitude_INIT
        Phi_field = noise - np.mean(noise) # Centered noise
    else: # xp == cupy
        noise = xp.random.rand(args.Nx, args.Ny, dtype=xp.float64) * args.Phi0_noise_amplitude_INIT
        Phi_field = noise - xp.mean(noise) # Centered noise
    print(f"INFO: Phi_field initialized with shape {Phi_field.shape} on {type(xp).__name__}")

    # Initialize ODE variables (as xp arrays for consistency, even if scalar)
    Psi_ID14 = xp.array(args.Psi0_INIT, dtype=xp.float64)
    I_ID14 = xp.array(args.I0_INIT, dtype=xp.float64)
    R_ID14 = xp.array(args.R0_INIT, dtype=xp.float64)
    print(f"INFO: Psi_ID14, I_ID14, R_ID14 initialized on {type(xp).__name__}")

    print_memory_usage("End of initialize_fields")
    return Phi_field, Psi_ID14, I_ID14, R_ID14

# --- Main Simulation Logic ---
def run_simulation(args):
    """Main function to run the simulation."""
    global XP_MODULE, SERIES_OF_IMAGES_DIR

    print_memory_usage("Start of run_simulation")

    # Determine backend (NumPy or CuPy)
    if args.use_gpu and cupy_available:
        XP_MODULE = cupy
        print("INFO: Using CuPy for GPU acceleration.")
        try:
            device = XP_MODULE.cuda.Device(0)
            device.use() # Select GPU 0
            try:
                # Attempt to get device name using runtime API
                device_props = XP_MODULE.cuda.runtime.getDeviceProperties(device.id)
                # device_props is a dict. The name is usually under the key 'name' and is in bytes.
                raw_name = device_props.get('name')
                if isinstance(raw_name, bytes):
                    device_name = raw_name.decode('utf-8')
                elif isinstance(raw_name, str):
                    device_name = raw_name # Should ideally be bytes, but handle if already str
                else:
                    # If 'name' key isn't found or not bytes/str, this will be caught by the outer except
                    raise ValueError(f"Device name not found or in unexpected format: {type(raw_name)}")
                print(f"INFO: CuPy is using GPU: {device_name} (ID: {device.id})")
            except Exception as e_get_name:
                print(f"INFO: CuPy is using GPU ID: {device.id} (Could not get name: {e_get_name})")
        except Exception as e:
            print(f"ERROR: Could not set CuPy device. {e}")
            print("INFO: Falling back to NumPy.")
            XP_MODULE = np
    else:
        XP_MODULE = np
        if args.use_gpu and not cupy_available:
            print("WARNING: --use_gpu was specified, but CuPy is not available. Using NumPy.")
        else:
            print("INFO: Using NumPy for CPU execution.")

    # Setup output directories
    MAIN_OUTPUT_DIR = args.MAIN_OUTPUT_DIR
    SERIES_OF_IMAGES_DIR = os.path.join(MAIN_OUTPUT_DIR, "series_of_images")
    try:
        os.makedirs(MAIN_OUTPUT_DIR, exist_ok=True)
        os.makedirs(SERIES_OF_IMAGES_DIR, exist_ok=True)
        print(f"INFO: Output will be saved to: {MAIN_OUTPUT_DIR}")
        print(f"INFO: Individual frames will be saved to: {SERIES_OF_IMAGES_DIR}")
    except OSError as e:
        print(f"ERROR: Could not create output directories: {e}", file=sys.stderr)
        sys.exit(1)

    PHI_SNAPSHOT_DIR = "" # Initialize
    if args.phi_snapshot_dir_name and args.snapshot_storage_skip_factor > 0:
        PHI_SNAPSHOT_DIR = os.path.join(MAIN_OUTPUT_DIR, args.phi_snapshot_dir_name)
        try:
            os.makedirs(PHI_SNAPSHOT_DIR, exist_ok=True)
            print(f"INFO: 2D Phi_field NPY snapshots will be saved to: {PHI_SNAPSHOT_DIR}")
        except Exception as e_mkdir_snapshot: # Catch generic Exception for snapshot dir
            print(f"WARNING: Could not create snapshot directory {PHI_SNAPSHOT_DIR}: {e_mkdir_snapshot}. NPY snapshots will be disabled.")
            PHI_SNAPSHOT_DIR = "" # Disable if creation fails
    else:
        print("INFO: 2D Phi_field NPY snapshots are disabled by configuration.")


    # Placeholder for further simulation setup (dx, dy, initial fields, etc.)
    print("INFO: Simulation setup (xp module, directories) complete.")

    # Calculate grid spacing
    dx = args.Lx / args.Nx
    dy = args.Ly / args.Ny
    print(f"INFO: dx = {dx}, dy = {dy}")

    # Initialize fields
    Phi_field, Psi_ID14, I_ID14, R_ID14 = initialize_fields(args, XP_MODULE)

    # Prepare lists for storing time-series data (on CPU)
    t_values_out = []
    psi_id14_values_out = []
    i_id14_values_out = []
    r_id14_values_out = []
    avg_phi_field_values_out = [] # To store average of Phi_field over time

    # Setup for Matplotlib animation
    fig, ax = plt.subplots(figsize=(8, 6))
    # Initial plot (empty or initial state)
    if XP_MODULE == np:
        initial_phi_display = Phi_field
    else:
        initial_phi_display = XP_MODULE.asnumpy(Phi_field)
    im = ax.imshow(initial_phi_display, animated=True, cmap=args.ANIMATION_CMAP, origin='lower',
                   extent=[0, args.Lx, 0, args.Ly])
    plt.colorbar(im, ax=ax, label='Phi_field')
    ax.set_title(f"Phi_field at t=0.00")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.tight_layout()

    frames_for_gif = [] # To store frames for GIF animation
    frame_idx = 0

    print_memory_usage("After full setup in run_simulation")

    # --- Main Time Loop ---
    print("\nINFO: Starting main simulation loop...")
    print_memory_usage("Before Main Loop")
    start_time_loop = time.time()

    num_steps = int(args.t_max / args.dt)
    for i in range(num_steps):
        t_current = i * args.dt

        # Calculate average of Phi_field (ensure it's a scalar on the correct device then CPU)
        Phi_avg_gpu = XP_MODULE.mean(Phi_field)
        if XP_MODULE != np:
            Phi_avg_cpu = float(Phi_avg_gpu.get()) # Get scalar from GPU
        else:
            Phi_avg_cpu = float(Phi_avg_gpu) # Already a CPU scalar

        # Get derivatives
        # Calculate intermediate fields for the classic PDE model
        Psi_c_s, Psi_d_s, rho_elec_field, R_rad_field = calculate_classic_intermediate_fields(Phi_field, Psi_ID14, args, XP_MODULE)

        # Calculate dPhi_dt using the classic PDE model
        dPhi_dt = pde_model(Phi_field, Psi_c_s, Psi_d_s, rho_elec_field, R_rad_field, args, dx, dy, XP_MODULE)
        dPsi_dt, dI_dt, dR_dt = ode_model(Psi_ID14, I_ID14, R_ID14, Phi_avg_gpu, args, XP_MODULE) # Pass GPU version of Phi_avg

        # Update fields using Euler integration
        Phi_field += dPhi_dt * args.dt
        if XP_MODULE.isnan(Phi_field).any() or XP_MODULE.isinf(Phi_field).any():
            print(f"WARNING: NaN or Inf detected in Phi_field at t={t_current:.2f} (step {i})")
            # Consider breaking or special handling if NaNs/Infs occur
            # For now, just print a warning and continue
        Psi_ID14 += dPsi_dt * args.dt
        I_ID14 += dI_dt * args.dt
        R_ID14 += dR_dt * args.dt

        # --- Save 2D Phi_field NPY snapshot (periodically) ---
        if PHI_SNAPSHOT_DIR and i % args.snapshot_storage_skip_factor == 0: # Check PHI_SNAPSHOT_DIR as it's set to "" if disabled
            try:
                # Ensure Phi_field is on CPU before saving
                phi_to_save_np = XP_MODULE.asnumpy(Phi_field) if XP_MODULE != np else Phi_field.copy()
                snapshot_filename = os.path.join(PHI_SNAPSHOT_DIR, f"phi_snapshot_step{i:07d}_t{t_current:.2f}.npy")
                np.save(snapshot_filename, phi_to_save_np)
                # print(f"DEBUG: Saved NPY snapshot {snapshot_filename}") # Uncomment for verbose debug
            except Exception as e_save_snapshot:
                print(f"WARNING: Could not save NPY snapshot at t={t_current:.2f} (step {i}): {e_save_snapshot}")

        # Data storage and animation update (periodically)
        if i % args.storage_skip_factor == 0:
            # Store scalar values (convert from xp.array to Python float)
            t_values_out.append(t_current)
            psi_id14_values_out.append(float(Psi_ID14.get()) if XP_MODULE != np else float(Psi_ID14))
            i_id14_values_out.append(float(I_ID14.get()) if XP_MODULE != np else float(I_ID14))
            r_id14_values_out.append(float(R_ID14.get()) if XP_MODULE != np else float(R_ID14))
            avg_phi_field_values_out.append(Phi_avg_cpu)

            # Update animation plot
            Phi_field_np = XP_MODULE.asnumpy(Phi_field) if XP_MODULE != np else Phi_field
            phi_min_np = np.min(Phi_field_np)
            phi_max_np = np.max(Phi_field_np)
            print(f"DEBUG_FRAME_DATA: t={t_current:.2f}, Phi_field_np min={phi_min_np:.4e}, max={phi_max_np:.4e}")

            im.set_data(Phi_field_np)
            if phi_min_np == phi_max_np:
                # Handle case where the field is flat (e.g., all NaNs, or uniform)
                im.set_clim(vmin=phi_min_np - 1e-5, vmax=phi_max_np + 1e-5) # Add a tiny epsilon if flat
            else:
                im.set_clim(vmin=phi_min_np, vmax=phi_max_np) # Dynamic scale for animation
            ax.set_title(f"Phi_field at t={t_current:.2f}")
            fig.canvas.draw_idle() # Update the plot
            # plt.pause(0.001) # Small pause to allow plot to update if running interactively (not with Agg)

            # Save individual frame as PNG
            frame_filename = os.path.join(SERIES_OF_IMAGES_DIR, f"frame_{frame_idx:04d}.png")
            try:
                plt.savefig(frame_filename)
                # print(f"DEBUG: Saved frame {frame_filename}") # Optional debug print
            except Exception as e_save_frame:
                print(f"WARNING: Could not save frame {frame_filename}: {e_save_frame}")

            # Save frame for GIF
            # Ensure data is on CPU and is a copy
            frames_for_gif.append(Phi_field_np.copy()) # Store a copy for the GIF
            frame_idx += 1

            # Progress print
            if i % (args.storage_skip_factor * 10) == 0: # Print less frequently than saving
                print(f"Progress: t = {t_current:.2f} / {args.t_max:.2f} (Step {i}/{num_steps})")
                sys.stdout.flush()

    end_time_loop = time.time()
    print_memory_usage("After Main Loop")
    print(f"INFO: Main simulation loop completed in {end_time_loop - start_time_loop:.2f} seconds.")


    # --- Post-Loop Processing ---
    print("\nINFO: Starting post-loop processing...")
    print_memory_usage("Before Post-Loop Processing")

    # 1. Save final state of Phi_field as an image
    final_state_filename = os.path.join(MAIN_OUTPUT_DIR, "phi_field_final_state.png")
    try:
        plt.figure(figsize=(8, 6))
        final_phi_display = XP_MODULE.asnumpy(Phi_field) if XP_MODULE != np else Phi_field
        plt.imshow(final_phi_display, cmap=args.FINAL_STATE_CMAP, origin='lower', extent=[0, args.Lx, 0, args.Ly])
        plt.colorbar(label='Phi_field')
        plt.title(f"Final State of Phi_field at t={args.t_max:.2f}")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.tight_layout()
        plt.savefig(final_state_filename)
        plt.close() # Close the figure to free memory
        print(f"INFO: Final state image saved to {final_state_filename}")
    except Exception as e_save_final:
        print(f"WARNING: Could not save final state image: {e_save_final}")

    # 2. Create and save GIF animation
    if frames_for_gif:
        animation_filename = os.path.join(MAIN_OUTPUT_DIR, "phi_field_animation.gif")
        try:
            # Reuse the existing figure and image object for animation
            # The fig and im are already set up from the live plotting part
            # We just need to re-create the animation object with all frames
            print(f"INFO: Creating animation with {len(frames_for_gif)} frames...")
            # Ensure figure is clean for animation generation
            ax.clear()
            im_anim = ax.imshow(frames_for_gif[0], animated=True, cmap=args.ANIMATION_CMAP, origin='lower',
                                extent=[0, args.Lx, 0, args.Ly])
            # Add colorbar again if cleared
            # plt.colorbar(im_anim, ax=ax, label='Phi_field') 
            # ^ This might add multiple colorbars if not handled carefully. 
            # For simplicity, we'll rely on the individual frames having consistent scaling.

            def update_anim(frame_num):
                im_anim.set_array(frames_for_gif[frame_num])
                ax.set_title(f"Phi_field at t={(frame_num * args.storage_skip_factor * args.dt):.2f}")
                return im_anim,

            anim = animation.FuncAnimation(fig, update_anim, frames=len(frames_for_gif), interval=100, blit=True)
            anim.save(animation_filename, writer='imagemagick', fps=10)
            plt.close(fig) # Close the animation figure
            print(f"INFO: Animation GIF saved to {animation_filename}")
        except Exception as e_save_gif:
            print(f"WARNING: Could not save animation GIF: {e_save_gif}")
            if "imagemagick" in str(e_save_gif).lower():
                print("HINT: Ensure ImageMagick is installed and in your system's PATH.")
    else:
        print("INFO: No frames collected for GIF animation.")

    # 3. Save time-series data to .npy files
    print("\nINFO: Saving time-series data to .npy files...")
    print_memory_usage("Before NPY Save Loop")

    required_lists_map = {
        "t_values": t_values_out,
        "psi_id14_values": psi_id14_values_out,
        "i_id14_values": i_id14_values_out,
        "r_id14_values": r_id14_values_out,
        "avg_phi_field_values": avg_phi_field_values_out
    }

    save_error_npy = False
    for name, data_list in required_lists_map.items():
        npy_filename = os.path.join(MAIN_OUTPUT_DIR, f"{name}.npy")
        try:
            if not isinstance(data_list, list):
                print(f"ERROR_NPY_SAVE [{name}]: Data for {name} is not a list (type: {type(data_list)}). Skipping.")
                save_error_npy = True
                continue
            
            # Convert list of Python floats/numbers to a NumPy array before saving
            np_array_to_save = np.array(data_list, dtype=np.float64)
            
            print(f"DEBUG_NPY_SAVE [{name}]: Attempting to save {name} (length: {len(data_list)}, type: {type(np_array_to_save)}, dtype: {np_array_to_save.dtype}) to {npy_filename}")
            np.save(npy_filename, np_array_to_save)
            print(f"INFO_NPY_SAVE [{name}]: Successfully saved {name} to {npy_filename}")
        except Exception as e_npy:
            print(f"ERROR_NPY_SAVE [{name}]: Failed to save {name} to {npy_filename}. Error: {e_npy}")
            save_error_npy = True

    if not save_error_npy:
        print("INFO: Successfully saved all time-series data to .npy files.")
    else:
        print("WARNING: One or more errors occurred during .npy file saving. Check logs above.")
    print_memory_usage("After NPY Save Loop")


    print_memory_usage("End of run_simulation")

if __name__ == '__main__':
    args = parse_arguments()
    print("--- Simulation Start (Version 2) ---")
    print(f"Arguments: {args}")

    # Placeholder for the rest of the script logic
    # Start tracemalloc for memory profiling
    tracemalloc.start()
    print_memory_usage("Initial after tracemalloc start")

    # Profiling setup
    profiler = cProfile.Profile()
    profiler.enable()

    try:
        run_simulation(args) # Call the main simulation function
    except Exception as e:
        print(f"FATAL ERROR in run_simulation: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
    finally:
        profiler.disable()
        print_memory_usage("After run_simulation call")

        # Print cProfile stats
        s = io.StringIO()
        ps = pstats.Stats(profiler, stream=s).sort_stats('cumulative')
        print("\n--- cProfile Stats (Top 20 by cumulative time) ---")
        ps.print_stats(20)
        print(s.getvalue())

        # Print final memory usage
        print("\n--- Final Memory Usage ---")
        print_memory_usage("Final before tracemalloc stop")
        tracemalloc.stop()

    print("\n--- Script End (Version 2) ---")
