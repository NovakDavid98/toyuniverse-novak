import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy.fft import rfft, rfftfreq
from scipy.signal import correlate, correlation_lags
import sys
import faulthandler
faulthandler.enable()

print("Starting simulation (stdout)...", flush=True)
sys.stderr.write("Starting simulation (stderr)...\n")
sys.stderr.flush()
try:
    with open("cascade_debug_output.txt", "w") as f_debug:
        f_debug.write("Script started, attempting to write to debug file.\n")
except Exception as e:
    pass # Silently pass as console might be broken

# --- Parameters ---

# ID14 Core Parameters (from toy_universe_simulation.py)
U_E = 86.4  # Fundamental Energy Unit
omega_image = 0.5184  # Characteristic Angular Frequency (2*pi / T_period where T_period = 12.122)
R_char = 8.0065  # Characteristic Radius
I_char = 52.20  # Characteristic Current

# ID14 Derived & Stabilization Constants (from toy_universe_simulation.py)
alpha_psi_ID14 = np.pi**2 / 6  # for dPsi/dt
beta_psi_ID14 = 2.0             # Was U_E / R_char (~10.79), for dPsi/dt
gamma_Psi_ID14 = 1.0          # Damping for Psi_ID14 stabilization (was 5.0)

# ID14 Coupling Constants (NEW - for coupling with ID26)
k_I_phi = 0.1     # Coupling from Phi to I_ID14 (was 1e-2)
k_R_phi_prime = 0.01 # Coupling from dPhi/dt to R_ID14 (was 1e-2)

# ID26 Parameters (NEW - all placeholders)
O_octave = 1.0 # Octave number
# For Phi_comp_field = (k_g1 * Phi + k_g2 * Psi_d) * 2^(O-1)
k_g1 = 1.0
k_g2 = 0.001  # Was 0.01, Coefficient for Phi^3 term in Phi_compression, reduced
k_rho = 1.0
# For V = k_v * (1 / rho_electric) * (6/pi) * Psi_m
k_v = 1.0
rho_epsilon = 1e-6 # Small value to prevent division by zero in V calculation
# For R_radiative = kappa_R * 2^(O-1) - k_R_int * Phi * Psi_c * V
kappa_R = 1.0  # Was 10.0
k_R_int = 1.0  # Interaction term for R_radiative, (re-enabled from 0.0)
# For PDE: dPhi/dt = alpha_Phi*Phi*Psi_c - beta_Phi*R_rad - gamma_Phi_PDE*Laplacian(Phi) + epsilon_Phi*Phi*Psi_d*rho_elec
alpha_Phi_PDE = 2.1    # Reduced from 2.5, increased from 2.0
beta_Phi_PDE = 0.0001  # Was 0.1, significantly reduced to lessen negative drive on Phi
gamma_Phi_PDE = 0.01 # Diffusion coefficient for Laplacian term
epsilon_Phi_PDE = 2.1   # Reduced from 2.5, increased from 2.0
Phi0_saturation_alpha = 0.2 # Saturation scale for Phi_field in alpha term (reduced from 1.0)
Phi0_saturation_epsilon = 0.2 # Saturation scale for Phi_field in epsilon term (reduced from 1.0)
Phi_coupling_saturation_limit = 0.5 # Max value for Phi_field average passed to ID14

# Grid Parameters (NEW)
N_x = 100  # Number of spatial points
x_min = 0.0
x_max = 10.0
dx = (x_max - x_min) / (N_x - 1)
x_points = np.linspace(x_min, x_max, N_x)

# Time Parameters
T_period_ID14 = 2 * np.pi / omega_image
t_max = 10 * T_period_ID14  # Simulate for 10 characteristic periods of ID14
dt = 0.001  # Time step for PDE and overall simulation. Needs to be small for PDE stability.
num_steps = int(t_max / dt)
print(f"DEBUG: num_steps = {num_steps}", flush=True)
storage_skip_factor = 100  # Store data every Nth actual simulation step
num_storage_points = (num_steps // storage_skip_factor) + 1

# --- Initialization ---

# ID14 Initial Conditions
I0 = 0.0
R0 = 0.0
Psi0_ID14 = 0.5 # Matches previous Psi0

I_current = I0
R_current = R0
Psi_ID14_current = Psi0_ID14

# ID26 Field Initial Conditions
# Phi_field = np.zeros(N_x)  # Initial electric potential field (e.g., all zeros)
Phi_field = 0.01 * np.sin(np.pi * x_points / x_max) # Example: a small sine wave perturbation

# Storage for results (downsampled)
t_values_out = np.zeros(num_storage_points)
I_values_out = np.zeros(num_storage_points)
R_values_out = np.zeros(num_storage_points)
Psi_ID14_values_out = np.zeros(num_storage_points)
Phi_field_out = np.zeros((num_storage_points, N_x))

t_values_out[0] = 0
I_values_out[0] = I_current
R_values_out[0] = R_current
Psi_ID14_values_out[0] = Psi_ID14_current
Phi_field_out[0, :] = Phi_field.copy()

# --- Function Definitions ---

def id14_model(t, y, Phi_from_id26_scalar, id14_core_params, id14_coupling_params):
    # print(f"t={t:.2f}, I={y[0]:.2e}, R={y[1]:.2e}, Psi={y[2]:.2e}, Phi_scalar_in={Phi_from_id26_scalar:.2e}") # DEBUG inside model
    I, R, Psi_ID14 = y
    U_E_p, omega_image_p, R_char_p, I_char_p, alpha_psi_p, beta_psi_p, gamma_Psi_p = id14_core_params
    k_I_phi_p, k_R_phi_prime_p = id14_coupling_params

    # dI/dt = k_I_phi * Phi_electric_ID14(t) - (U_E / I_char) * I * (1 + cos(omega_image * t))
    dIdt = k_I_phi_p * Phi_from_id26_scalar - (U_E_p / I_char_p) * I * (1 + np.cos(omega_image_p * t))

    # dR/dt = k_R_phi_prime * Phi_electric_ID14(t) * I * sin(omega_image * t) - (2 * U_E / R_char) * R
    dRdt = k_R_phi_prime_p * Phi_from_id26_scalar * I * np.sin(omega_image_p * t) - (2 * U_E_p / R_char_p) * R

    # dPsi_ID14/dt = alpha * I - beta * R - gamma_Psi * Psi_ID14
    dPsidt = alpha_psi_p * I - beta_psi_p * R - gamma_Psi_p * Psi_ID14

    return [dIdt, dRdt, dPsidt]

def calculate_laplacian(field, dx_val):
    laplacian = np.zeros_like(field)
    laplacian[1:-1] = (field[:-2] - 2 * field[1:-1] + field[2:]) / (dx_val**2)
    # Boundary conditions for laplacian (e.g., zero at boundaries if field is zero there)
    # Or handle Neumann by setting laplacian[0] and laplacian[-1] based on field values
    laplacian[0] = (field[1] - 2 * field[0] + field[0]) / (dx_val**2) # Assuming Phi(x<0) = Phi(0) (zero gradient approx for first point)
    laplacian[-1] = (field[-2] - 2 * field[-1] + field[-1]) / (dx_val**2) # Assuming Phi(x>x_max) = Phi(x_max)
    return laplacian

def calculate_id26_intermediate_fields(Phi_field_current, Psi_scalar_current, id26_p, grid_p):
    O, k_g1_p, k_g2_p, k_rho_p, k_v_p, rho_eps_p, kappa_R_p, k_R_int_p = id26_p
    # Assuming Psi_c = Psi_d = Psi_m = Psi_scalar_current for simplicity. Apply tanh saturation.
    # Psi_c = Psi_ID14 / Psi0_ID14
    # Psi_d = Psi_ID14 * Psi0_ID14
    Psi_c_s = np.tanh(Psi_scalar_current / Psi0_ID14) # Using global Psi0_ID14 (defined as 0.5)
    Psi_d_s = np.tanh(Psi_scalar_current * Psi0_ID14) # Using global Psi0_ID14 (defined as 0.5)
    Psi_m_s = Psi_scalar_current # Psi_m is not typically scaled by Psi0 or tanh'd here

    # A. Compression Potential Component
    Phi_comp_field = (k_g1_p * Phi_field_current + k_g2_p * Psi_d_s) * (2**(O-1))

    # B. Electric Charge Density
    rho_electric_field = k_rho_p * Phi_comp_field * Psi_c_s

    # C. Velocity (handle potential division by zero)
    V_field = k_v_p * (1 / (rho_electric_field + rho_eps_p)) * (6 / np.pi) * Psi_m_s

    # D. Radiative Radius/Flux
    R_radiative_field = kappa_R_p * (2**(O-1)) - k_R_int_p * Phi_field_current * Psi_c_s * V_field
    
    return Phi_comp_field, rho_electric_field, V_field, R_radiative_field

def id26_pde_rhs(Phi_curr, Psi_c_s, Psi_d_s, rho_elec_curr, R_rad_curr, laplacian_Phi_curr, id26_pde_consts_extended):
    alpha_p, beta_p, gamma_p, epsilon_p, phi0_sat_a, phi0_sat_e = id26_pde_consts_extended

    # Apply saturation to Phi_field's contribution to its own growth terms
    effective_Phi_alpha = phi0_sat_a * np.tanh(Phi_curr / phi0_sat_a)
    term1 = alpha_p * effective_Phi_alpha * Psi_c_s

    term2 = -beta_p * R_rad_curr
    term3 = gamma_p * laplacian_Phi_curr  # Diffusion term

    effective_Phi_epsilon = phi0_sat_e * np.tanh(Phi_curr / phi0_sat_e)
    term4 = epsilon_p * effective_Phi_epsilon * Psi_d_s * rho_elec_curr
    
    return term1 + term2 + term3 + term4

# --- Main Simulation Loop ---
print("Starting simulation...", flush=True)

id14_params_core = (U_E, omega_image, R_char, I_char, alpha_psi_ID14, beta_psi_ID14, gamma_Psi_ID14)
id14_params_coupling = (k_I_phi, k_R_phi_prime)
id26_params_intermediate = (O_octave, k_g1, k_g2, k_rho, k_v, rho_epsilon, kappa_R, k_R_int)
id26_params_pde = (alpha_Phi_PDE, beta_Phi_PDE, gamma_Phi_PDE, epsilon_Phi_PDE, Phi0_saturation_alpha, Phi0_saturation_epsilon)

print("DEBUG: Functions defined, loop params set. Pre-loop.", flush=True)

for i in range(num_steps):
    print(f"DEBUG: Starting loop iteration i={i} (step {i+1})", flush=True)
    current_t = i * dt

    # 1. Extract Phi_for_ID14 from Phi_field (e.g., center value or average)
    # Phi_for_ID14_scalar = Phi_field[N_x // 2] 
    Phi_for_ID14_raw = np.mean(Phi_field)
    Phi_for_ID14_scalar = Phi_coupling_saturation_limit * np.tanh(Phi_for_ID14_raw / Phi_coupling_saturation_limit)
    # print(f"Loop {i}: t={current_t:.2f}, Pre-ID14: I={I_current:.2e}, R={R_current:.2e}, Psi={Psi_ID14_current:.2e}, Phi_to_ID14={Phi_for_ID14_scalar:.2e}") # DEBUG

    # 2. Solve ID14 ODEs for one step dt
    # Using solve_ivp for one step. More accurate but slower than Euler.
    print("DEBUG: Before solve_ivp", flush=True)
    sol = solve_ivp(id14_model, [current_t, current_t + dt], 
                    [I_current, R_current, Psi_ID14_current], 
                    args=(Phi_for_ID14_scalar, id14_params_core, id14_params_coupling),
                    dense_output=False, t_eval=[current_t + dt])
    print(f"DEBUG: After solve_ivp, sol.success={sol.success}", flush=True)
    
    if not sol.success:
        print(f"ERROR: ID14 solve_ivp failed at t={current_t:.2f}! Status: {sol.status}, Message: {sol.message}")
        print(f"       Input y0: I={I_current:.2e}, R={R_current:.2e}, Psi={Psi_ID14_current:.2e}")
        print(f"       Input Phi_scalar: {Phi_for_ID14_scalar:.2e}")
        # Potentially break or handle error
        I_current, R_current, Psi_ID14_current = np.nan, np.nan, np.nan # Propagate failure
    else:
        I_current = sol.y[0, -1]
        R_current = sol.y[1, -1]
        Psi_ID14_current = sol.y[2, -1]

    if np.isnan(I_current) or np.isnan(R_current) or np.isnan(Psi_ID14_current):
        print(f"ERROR: NaN detected in ID14 solution at t={current_t+dt:.2f}. Stopping.")
        break

    # Update Psi_c, Psi_d, Psi_m (assuming they are all = Psi_ID14_current)
    # Apply tanh saturation to Psi_c and Psi_d before passing to ID26 calculations.
    # Psi_c = Psi_ID14 / Psi0_ID14
    # Psi_d = Psi_ID14 * Psi0_ID14
    current_Psi_c_scalar = np.tanh(Psi_ID14_current / Psi0_ID14) # Using global Psi0_ID14
    current_Psi_d_scalar = np.tanh(Psi_ID14_current * Psi0_ID14) # Using global Psi0_ID14
    # current_Psi_m_scalar = Psi_ID14_current # V_field uses this, Psi_m not usually tanh'd here

    # 4. Calculate ID26 intermediate fields
    _, current_rho_electric_field, _, current_R_radiative_field = \
        calculate_id26_intermediate_fields(Phi_field, current_Psi_c_scalar, 
                                           id26_params_intermediate, None) # grid_params not used yet

    # 5. Calculate Laplacian of Phi_field
    laplacian_Phi = calculate_laplacian(Phi_field, dx)

    # 6. Calculate RHS of ID26 PDE for all x
    pde_rhs_values = id26_pde_rhs(Phi_field, current_Psi_c_scalar, current_Psi_d_scalar, 
                                  current_rho_electric_field, current_R_radiative_field, 
                                  laplacian_Phi, id26_params_pde)

    # 7. Update Phi_field using Euler step
    Phi_field_new = Phi_field + dt * pde_rhs_values

    # 8. Apply Boundary Conditions to Phi_field_new
    # The calculate_laplacian function implements Neumann (zero-flux) conditions.
    # Forcing Dirichlet here (Phi_field_new[0]=0, Phi_field_new[-1]=0) overrides that
    # and is inconsistent with previously observed stable results. Commenting out.
    # Phi_field_new[0] = 0.0  # Fixed at x_min
    # Phi_field_new[-1] = 0.0 # Fixed at x_max
    Phi_field = Phi_field_new

    # 9. Store results (downsampled)
    # Note: i is 0-indexed for steps 0 to num_steps-1. (i+1) is the current step number.
    if (i + 1) % storage_skip_factor == 0:
        storage_idx = (i + 1) // storage_skip_factor
        if storage_idx < num_storage_points: # Ensure we don't write out of bounds if num_steps isn't perfectly divisible
            t_values_out[storage_idx] = current_t + dt
            I_values_out[storage_idx] = I_current
            R_values_out[storage_idx] = R_current
            Psi_ID14_values_out[storage_idx] = Psi_ID14_current
            Phi_field_out[storage_idx, :] = Phi_field.copy()

    if (i+1) % (num_steps // 100) == 0:
        print(f"DEBUG: Entered progress print block for step {i+1}.", flush=True)
        avg_phi_val_str = "ErrorInAvgPhi" # Default in case of issues
        # Check for numerical instability before calling np.mean in the progress print
        if np.isnan(Phi_field).any() or np.isinf(Phi_field).any():
            print(f"WARNING: NaN or Inf detected in Phi_field at t={current_t + dt:.2f} (step {i+1}) *immediately before* avg(Phi) calculation in progress print.", flush=True)
            avg_phi_val_str = "NaN/Inf"
        else:
            try:
                avg_phi_val_str = f"{np.mean(Phi_field):.2e}"
            except Exception as e_mean:
                print(f"ERROR calculating np.mean(Phi_field) at step {i+1}: {e_mean}", flush=True)
        
        print(f"Progress: {((i+1)/num_steps*100):.0f}%, t={current_t + dt:.2f}, Psi_ID14={Psi_ID14_current:.6e}, avg(Phi)={avg_phi_val_str}", flush=True)
        print(f"DEBUG: Completed progress print statement for step {i+1}.", flush=True)

print("DEBUG: Main simulation loop finished.", flush=True)

print("Simulation complete.")

# --- Quantitative Analysis of Stable Dynamics ---
print("\n--- Quantitative Analysis of Stable Dynamics ---")

t_start_active_phase = 80.0  # Define start of active phase for temporal analysis
dt_main_loop = dt * storage_skip_factor # Adjusted time step for downsampled data

# Helper function to get active phase data
def get_active_phase_data(t_values, data_values, t_start):
    active_indices = t_values >= t_start
    return t_values[active_indices], data_values[active_indices]

# 1. Temporal Analysis: Dominant Frequency and Amplitude
def analyze_temporal_oscillations(t_values, data_values, series_name, t_start, dt_sim):
    print(f"\nAttempting Temporal Analysis for {series_name}...")
    try:
        t_active, data_active = get_active_phase_data(t_values, data_values, t_start)

        if len(data_active) < 2:
            print(f"  Not enough data in active phase for {series_name} analysis (found {len(data_active)} points).")
            return

        # Amplitude
        peak_to_trough_amplitude = np.max(data_active) - np.min(data_active)
        print(f"  Peak-to-Trough Amplitude ({series_name}, active phase, t>{t_start:.1f}): {peak_to_trough_amplitude:.4e}")

        # Dominant Frequency (FFT)
        N_active = len(data_active)
        if N_active > 1 and dt_sim > 0:
            yf = rfft(data_active - np.mean(data_active))  # Remove DC offset
            xf = rfftfreq(N_active, dt_sim)
            if len(xf) > 1 and len(yf) > 1:
                dominant_freq_idx = np.argmax(np.abs(yf[1:])) + 1  # Exclude DC component
                if dominant_freq_idx < len(xf):
                    dominant_freq = xf[dominant_freq_idx]
                    print(f"  Dominant Frequency ({series_name}, active phase, t>{t_start:.1f}): {dominant_freq:.4f} Hz (approx.)")
                else:
                    print(f"  FFT Error ({series_name}): Dominant frequency index out of bounds.")
            else:
                print(f"  Not enough data points for meaningful FFT frequency analysis after filtering ({series_name}). len(xf)={len(xf)}, len(yf)={len(yf)}.")
        else:
            print(f"  Cannot perform FFT ({series_name}): Not enough data (N_active={N_active}) or invalid dt_sim (dt_sim={dt_sim}).")
        print(f"Completed Temporal Analysis for {series_name}.")
    except Exception as e:
        print(f"ERROR during temporal analysis for {series_name}: {e}")

# 2. Temporal Analysis: Phase Relationships (Cross-Correlation)
def analyze_phase_relationship(t_values, data1, data2, name1, name2, t_start, dt_sim):
    print(f"\nAttempting Phase Relationship Analysis between {name1} and {name2}...")
    try:
        t_active, d1_active = get_active_phase_data(t_values, data1, t_start)
        _, d2_active = get_active_phase_data(t_values, data2, t_start)

        if len(d1_active) < 2 or len(d2_active) < 2:
            print(f"  Not enough data in active phase for cross-correlation ({name1} vs {name2}). Found {len(d1_active)}, {len(d2_active)} points.")
            return

        # Normalize data for cross-correlation
        std_d1 = np.std(d1_active)
        std_d2 = np.std(d2_active)

        if std_d1 < 1e-9 or std_d2 < 1e-9: # Check for near-zero std dev
            print(f"  Cannot normalize data for cross-correlation ({name1} vs {name2}) - std is near zero for one or both series in active phase. Std1: {std_d1:.2e}, Std2: {std_d2:.2e}")
            return
            
        d1_norm = (d1_active - np.mean(d1_active)) / std_d1
        d2_norm = (d2_active - np.mean(d2_active)) / std_d2
        
        correlation = correlate(d1_norm, d2_norm, mode='full')
        lags_samples = correlation_lags(len(d1_norm), len(d2_norm), mode='full')
        
        lag_at_max_corr_samples = lags_samples[np.argmax(np.abs(correlation))] # Use abs for max magnitude of correlation
        lag_at_max_corr_time = lag_at_max_corr_samples * dt_sim

        print(f"  Lag at max correlation magnitude ({name1} vs {name2}): {lag_at_max_corr_time:.4f} time units.")
        if lag_at_max_corr_time < -1e-9: # Comparing floats, allow for small epsilon
            print(f"    Interpretation: {name1} lags {name2} by approx. {abs(lag_at_max_corr_time):.4f} time units.")
        elif lag_at_max_corr_time > 1e-9:
            print(f"    Interpretation: {name1} leads {name2} by approx. {lag_at_max_corr_time:.4f} time units.")
        else:
            print(f"    Interpretation: {name1} and {name2} are approximately in phase.")
        print(f"Completed Phase Relationship Analysis for {name1} vs {name2}.")
    except Exception as e:
        print(f"ERROR during phase relationship analysis for {name1} vs {name2}: {e}")

# 3. Spatial Analysis of Phi_field
def analyze_phi_spatial_profile(x_coords, phi_profile, profile_time):
    print(f"\nAttempting Spatial Analysis of Phi_field at t = {profile_time:.2f}...")
    try:
        if len(phi_profile) != len(x_coords) or len(phi_profile) == 0:
            print(f"  Spatial Analysis Error: phi_profile (len {len(phi_profile)}) and x_coords (len {len(x_coords)}) mismatch or empty.")
            return

        peak_value = np.max(phi_profile)
        peak_idx = np.argmax(phi_profile)
        peak_location_x = x_coords[peak_idx]
        print(f"  Peak Phi value: {peak_value:.4f} at x = {peak_location_x:.4f}")

        if peak_value <= 1e-9: # Avoid FWHM for zero or negative peaks
            print("  Peak value is near zero or negative, FWHM not calculated.")
            print(f"Completed Spatial Analysis for t={profile_time:.2f} (FWHM not applicable).")
            return

        half_max_val = peak_value / 2.0
        above_half_max_indices = np.where(phi_profile >= half_max_val)[0]

        if len(above_half_max_indices) == 0:
            print("  Profile does not exceed half maximum, FWHM not calculated.")
            print(f"Completed Spatial Analysis for t={profile_time:.2f} (FWHM not applicable).")
            return

        first_idx_above = above_half_max_indices[0]
        last_idx_above = above_half_max_indices[-1]

        # Interpolate left boundary of FWHM
        x_left_half = x_coords[first_idx_above]
        if first_idx_above > 0 and phi_profile[first_idx_above - 1] < half_max_val:
            # Ensure interpolation points are monotonic for phi_profile values
            phi_interp_left = [phi_profile[first_idx_above - 1], phi_profile[first_idx_above]]
            x_interp_left = [x_coords[first_idx_above - 1], x_coords[first_idx_above]]
            if phi_interp_left[0] < phi_interp_left[1]: # Normal case
                 x_left_half = np.interp(half_max_val, phi_interp_left, x_interp_left)
            elif phi_interp_left[0] > phi_interp_left[1]: # Flipped for interp if needed (should not happen if peak is clear)
                 x_left_half = np.interp(half_max_val, phi_interp_left[::-1], x_interp_left[::-1])
            # else: # phi_interp_left[0] == phi_interp_left[1], use original x_coords[first_idx_above]

        # Interpolate right boundary of FWHM
        x_right_half = x_coords[last_idx_above]
        if last_idx_above < len(x_coords) - 1 and phi_profile[last_idx_above + 1] < half_max_val:
            # Ensure interpolation points are monotonic for phi_profile values
            phi_interp_right = [phi_profile[last_idx_above + 1], phi_profile[last_idx_above]] # Note order for interpolation from outside to inside peak
            x_interp_right = [x_coords[last_idx_above + 1], x_coords[last_idx_above]]
            if phi_interp_right[0] < phi_interp_right[1]: # Normal case
                x_right_half = np.interp(half_max_val, phi_interp_right, x_interp_right)
            elif phi_interp_right[0] > phi_interp_right[1]: # Flipped
                x_right_half = np.interp(half_max_val, phi_interp_right[::-1], x_interp_right[::-1])
            # else: use original x_coords[last_idx_above]
        
        fwhm = x_right_half - x_left_half
        if fwhm < 0: # Should not happen with correct interpolation logic
            print(f"  Warning: Calculated FWHM is negative ({fwhm:.4f}). Defaulting to index-based width.")
            fwhm = x_coords[last_idx_above] - x_coords[first_idx_above]
        
        print(f"  Estimated FWHM of central peak: {fwhm:.4f}")
        print(f"Completed Spatial Analysis for t={profile_time:.2f}.")

    except IndexError as ie:
        print(f"  Could not reliably estimate FWHM for this profile at t={profile_time:.2f} (IndexError: {ie}).")
    except Exception as e:
        print(f"  Error during FWHM estimation at t={profile_time:.2f}: {e}")

print("DEBUG: Starting analysis section...", flush=True)
# --- Perform Analysis ---
# Ensure all necessary data arrays are available and correctly named from the simulation output

analyze_temporal_oscillations(t_values_out, Psi_ID14_values_out, "Psi_ID14", t_start_active_phase, dt_main_loop)
analyze_temporal_oscillations(t_values_out, I_values_out, "I_ID14", t_start_active_phase, dt_main_loop)
analyze_temporal_oscillations(t_values_out, R_values_out, "R_ID14", t_start_active_phase, dt_main_loop)

avg_Phi_field_out = np.mean(Phi_field_out, axis=1)
analyze_temporal_oscillations(t_values_out, avg_Phi_field_out, "avg(Phi_field)", t_start_active_phase, dt_main_loop)

analyze_phase_relationship(t_values_out, Psi_ID14_values_out, I_values_out, "Psi_ID14", "I_ID14", t_start_active_phase, dt_main_loop)
analyze_phase_relationship(t_values_out, avg_Phi_field_out, Psi_ID14_values_out, "avg(Phi_field)", "Psi_ID14", t_start_active_phase, dt_main_loop)

x_coordinates = np.linspace(x_min, x_max, N_x)
times_for_spatial_analysis = [90.0, 100.0, 110.0, 120.0]
for t_spatial in times_for_spatial_analysis:
    idx_spatial = np.argmin(np.abs(t_values_out - t_spatial))
    actual_time_spatial = t_values_out[idx_spatial]
    phi_profile_snapshot = Phi_field_out[idx_spatial, :]
    analyze_phi_spatial_profile(x_coordinates, phi_profile_snapshot, actual_time_spatial)

print("\n--- End of Quantitative Analysis ---")
print("DEBUG: Finished analysis section.", flush=True)


# --- Plotting ---
print("Plotting results...")

fig, axs = plt.subplots(4, 1, figsize=(12, 10), sharex=True)

# Plot I(t)
axs[0].plot(t_values_out, I_values_out, label='I(t)')
axs[0].set_ylabel('I (Current)')
axs[0].grid(True)
axs[0].legend()

# Plot R(t)
axs[1].plot(t_values_out, R_values_out, label='R(t)')
axs[1].set_ylabel('R (Radius)')
axs[1].grid(True)
axs[1].legend()

# Plot Psi_ID14(t)
axs[2].plot(t_values_out, Psi_ID14_values_out, label='Psi_ID14(t)')
axs[2].set_ylabel('Psi_ID14 (Consciousness)')
axs[2].grid(True)
axs[2].legend()

# Plot average Phi(x,t) over time
avg_Phi_t = np.mean(Phi_field_out, axis=1)
axs[3].plot(t_values_out, avg_Phi_t, label='avg(Phi(x,t)) over x')
axs[3].set_ylabel('Average Phi Field')
axs[3].set_xlabel('Time (model units)')
axs[3].grid(True)
axs[3].legend()

plt.tight_layout()
plot_filename = 'id14_id26_simulation_plot_summary.png'
plt.savefig(plot_filename)
print("DEBUG: Summary plot saved.", flush=True)
print(f"Summary plot saved to {plot_filename}")
plt.close(fig)

# Plot Phi(x,t) as a 2D heatmap
plt.figure(figsize=(10, 6))
plt.imshow(Phi_field_out.T, aspect='auto', origin='lower', 
           extent=[0, t_max, x_min, x_max], cmap='viridis') # Transpose for x on y-axis, t on x-axis
plt.colorbar(label='Phi(x,t) value')
plt.xlabel('Time (model units)')
plt.ylabel('Spatial Dimension x')
plt.title('ID26 Electric Potential Field Phi(x,t)')
heatmap_filename = 'id14_id26_simulation_plot_phi_heatmap.png'
plt.savefig(heatmap_filename)
print("DEBUG: Heatmap plot saved.", flush=True)
print(f"Phi field heatmap saved to {heatmap_filename}")
plt.close()

print("Done.")
print("--- Script End ---", flush=True)
