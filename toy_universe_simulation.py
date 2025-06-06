import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# --- Core Parameters from the Research Paper ---
U_E = 86.4  # Fundamental Energy Unit
omega_image = 0.5184  # Characteristic Angular Frequency
# R_char = (U_E**3 / (400 * np.pi))**(1/3) # Characteristic System Radius
R_char = 8.00651991 # Using the pre-calculated value from the paper for consistency
I_char = 52.20  # Characteristic Current (also I_char_param from simulation notes)

# --- Derived Constants for ID14 ---
alpha = np.pi**2 / 6
beta = U_E / R_char
gamma_Psi = 5.0  # Damping coefficient for Psi stabilization (Targeting avg(Psi) ~14)
phi_stillpoint = 0.5 # Assumption for the Φ term in dR/dt, consistent with Ψ(0)

# --- System of Differential Equations (ID14) ---
def model(t, y):
    I, R, Psi = y

    # REVISED dI/dt: Halved constant driving term U_E to target avg(I) ~ I_char/2
    # REVISED dI/dt: Further refined constant driving term to (U_E / 3.54) to target avg(I) ~ I_char/2
    # dI/dt = (U_E / 3.54) - (U_E / I_char) * I(t) * (1 + cos(ω*t))
    dIdt = (U_E / 3.54) - (U_E / I_char) * I * (1 + np.cos(omega_image * t))

    # REVISED dR/dt: Based on ID14 structure with oscillatory drive and ordinal factor O=2 for damping.
    # Φ_electric_effective_value = 2 * U_E (hypothesis for amplitude scaling)
    # k₂ = 1/I_char
    # γ_R_effective_damping = 2 * (U_E / R_char)
    # dR/dt = (Φ_electric_effective_value / I_char) * I * sin(ω*t) - γ_R_effective_damping * R
    phi_E_val = 2.0 * U_E
    dRdt = (phi_E_val / I_char) * I * np.sin(omega_image * t) - 2.0 * (U_E / R_char) * R

    # dΨ/dt with stabilization term
    dPsidt = alpha * I - beta * R - gamma_Psi * Psi

    return [dIdt, dRdt, dPsidt]

# --- Initial Conditions ---
I0 = 0.0
R0 = 0.0 # Starts at its mean for a stable oscillation around 0
Psi0 = phi_stillpoint # Start at the 'stillpoint' value from ID22
y0 = [I0, R0, Psi0]

# --- Time Span for Simulation ---
T_period = 1 / omega_image
t_max = 50 * T_period  # Simulate for 50 characteristic periods
t_eval_points = 1000 # Number of points for a smooth plot
t_eval = np.linspace(0, t_max, t_eval_points)

# --- Solve the ODE System ---
# Using 'RK45' which is a good general-purpose solver. Max_step can prevent too large steps.
sol = solve_ivp(model, [0, t_max], y0, t_eval=t_eval, dense_output=True, method='RK45', max_step=T_period/20)

# --- Extract Results ---
t = sol.t
I_t = sol.y[0]
R_t = sol.y[1]
Psi_t = sol.y[2]

# --- Plotting Results ---
fig, axs = plt.subplots(3, 1, figsize=(12, 15), sharex=True)

axs[0].plot(t, I_t, label='I(t) - Current')
axs[0].axhline(y=I_char/2, color='r', linestyle='--', label=f'Target Mean I(t) = {I_char/2:.2f}')
axs[0].axhline(y=I_char, color='g', linestyle=':', label=f'I_char = {I_char:.2f}')
axs[0].set_ylabel('I(t) [model current units]')
axs[0].set_title(f'Toy Universe: ID14 System Dynamics (U_E={U_E}, ω={omega_image:.4f})')
axs[0].grid(True)
axs[0].legend()

axs[1].plot(t, R_t, label='R(t) - Radius', color='orange')
axs[1].axhline(y=0, color='r', linestyle='--', label='Target Mean R(t) = 0')
axs[1].set_ylabel('R(t) [model length units]')
axs[1].grid(True)
axs[1].legend()

axs[2].plot(t, Psi_t, label='Ψ(t) - Consciousness', color='green')
axs[2].axhline(y=phi_stillpoint, color='r', linestyle='--', label=f'Ψ_stillpoint = {phi_stillpoint:.2f}')
axs[2].set_ylabel('Ψ(t) [model consciousness units]')
axs[2].set_xlabel('Time [model time units]')
axs[2].grid(True)
axs[2].legend()

plt.tight_layout()
plot_filename = 'id14_simulation_plot.png'
plt.savefig(plot_filename)
print(f"Simulation complete. Plot saved to {plot_filename}")

# --- Optional: Print some summary values ---
print(f"Core Parameters: U_E={U_E}, ω_image={omega_image:.4f}, R_char={R_char:.4f}, I_char={I_char:.2f}")
print(f"Derived Constants: α={alpha:.4f}, β={beta:.4f}, γ_Ψ={gamma_Psi:.4f}")
print(f"Initial Conditions: I0={I0}, R0={R0}, Psi0={Psi0}")
print(f"Simulation Time Span: 0 to {t_max:.2f} model time units ({t_max/T_period:.1f} periods)")
print(f"Number of time steps evaluated: {len(t)}")
print(f"Final values @ t={t[-1]:.2f}: I={I_t[-1]:.2f}, R={R_t[-1]:.2f}, Psi={Psi_t[-1]:.2f}")
print(f"Mean values over simulation: avg(I)={np.mean(I_t):.2f}, avg(R)={np.mean(R_t):.2f}, avg(Psi)={np.mean(Psi_t):.2f}")
