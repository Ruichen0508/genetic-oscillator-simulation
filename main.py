import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def Genetic_oscillator(t, y, params):
    D_A, D_R, D_A_prime, D_R_prime, M_A, A, M_R, R, C = y
    alpha_A, alpha_A_prime, alpha_R, alpha_R_prime, beta_A, beta_R, delta_MA, delta_MR, delta_A, delta_R, gamma_A, gamma_R, gamma_C, theta_A, theta_R = params

    dD_A_dt = theta_A * D_A_prime - gamma_A * D_A * A
    dD_R_dt = theta_R * D_R_prime - gamma_R * D_R * A
    dD_A_prime_dt = gamma_A * D_A * A - theta_A * D_A_prime
    dD_R_prime_dt = gamma_R * D_R * A - theta_R * D_R_prime
    dM_A_dt = alpha_A_prime * D_A_prime + alpha_A * D_A - delta_MA * M_A
    dA_dt = beta_A * M_A + theta_A * D_A_prime + theta_R * D_R_prime - A * (gamma_A * D_A + gamma_R * D_R + gamma_C * R + delta_A)
    dM_R_dt = alpha_R_prime * D_R_prime + alpha_R * D_R - delta_MR * M_R
    dR_dt = beta_R * M_R - gamma_C * A * R + delta_A * C - delta_R * R
    dC_dt = gamma_C * A * R - delta_A * C

    return [dD_A_dt, dD_R_dt, dD_A_prime_dt, dD_R_prime_dt, dM_A_dt, dA_dt, dM_R_dt, dR_dt, dC_dt]

initial_parameters = [
    50,   # alpha_A: 1/h
    500,  # alpha_A_prime: 1/h
    0.01, # alpha_R: 1/h
    50,   # alpha_R_prime: 1/h
    50,   # beta_A: 1/h
    5,    # beta_R: 1/h
    10,   # delta_MA: 1/h
    0.5,  # delta_MR: 1/h
    1,    # delta_A: 1/h
    0.2,  # delta_R: 1/h
    1,    # gamma_A: 1/(mol*h)
    1,    # gamma_R: 1/(mol*h)
    2,    # gamma_C: 1/(mol*h)
    50,   # theta_A: 1/h
    100   # theta_R: 1/h
]

initial_conditions = [1, 1, 0, 0, 0, 0, 0, 0, 0]

t_span = (0, 400)
t_eval = np.linspace(0, 400, 1000)

sol_RK45 = solve_ivp(Genetic_oscillator, t_span, initial_conditions, method='RK45', args=(initial_parameters,), t_eval=t_eval)
sol_BDF = solve_ivp(Genetic_oscillator, t_span, initial_conditions, method='BDF', args=(initial_parameters,), t_eval=t_eval)
sol_Radau = solve_ivp(Genetic_oscillator, t_span, initial_conditions, method='Radau', args=(initial_parameters,), t_eval=t_eval)

plt.figure(figsize=(12, 6))

plt.subplot(1, 3, 1)
plt.plot(sol_RK45.t, sol_RK45.y[5], label='Activator A')
plt.plot(sol_RK45.t, sol_RK45.y[7], label='Repressor R')
plt.title('RK45 Solver')
plt.xlabel('Time (hours)')
plt.ylabel('Concentration')
plt.legend()

plt.subplot(1, 3, 2)
plt.plot(sol_BDF.t, sol_BDF.y[5], label='Activator A')
plt.plot(sol_BDF.t, sol_BDF.y[7], label='Repressor R')
plt.title('BDF Solver')
plt.xlabel('Time (hours)')
plt.ylabel('Concentration')
plt.legend()

plt.subplot(1, 3, 3)
plt.plot(sol_Radau.t, sol_Radau.y[5], label='Activator A')
plt.plot(sol_Radau.t, sol_Radau.y[7], label='Repressor R')
plt.title('Radau Solver')
plt.xlabel('Time (hours)')
plt.ylabel('Concentration')
plt.legend()
plt.savefig('Different_Solvers.png')

plt.tight_layout()
plt.show()

#Task 2
# Parameters for the system
params = {
    'alpha_A': 50,           # 1/hr
    'alpha_A_prime': 500,    # 1/hr
    'alpha_R': 0.01,         # 1/hr
    'alpha_R_prime': 50,     # 1/hr
    'beta_A': 50,            # 1/hr
    'beta_R': 5,             # 1/hr
    'delta_MA': 10,          # 1/hr
    'delta_MR': 0.5,         # 1/hr
    'delta_A': 1,            # 1/hr
    'delta_R': 0.2,          # 1/hr
    'gamma_A': 1,            # 1/(mol * hr)
    'gamma_R': 1,            # 1/(mol * hr)
    'gamma_C': 2,            # 1/(mol * hr)
    'theta_A': 50,           # 1/hr
    'theta_R': 100           # 1/hr
}

initial_state = [1, 1, 0, 0, 0, 0, 0, 0, 0]

StateChangeMat = np.array([
    [0, 0, 0, 0, 0, -1, 0, -1, +1],  # A + R -> C
    [0, 0, 0, 0, 0, -1, 0, 0, 0],    # A -> ∅
    [0, 0, 0, 0, 0, 0, 0, +1, -1],   # C -> R
    [0, 0, 0, 0, 0, 0, 0, -1, 0],    # R -> ∅
    [-1, 0, +1, 0, 0, -1, 0, 0, 0],  # D_A + A -> D'_A
    [0, -1, 0, +1, 0, -1, 0, 0, 0],  # D_R + A -> D'_R
    [+1, 0, -1, 0, 0, +1, 0, 0, 0],  # D'_A -> D_A + A
    [0, 0, 0, 0, +1, 0, 0, 0, 0],    # D_A -> D_A + M_A
    [0, 0, 0, 0, +1, 0, 0, 0, 0],    # D'_A -> D'_A + M_A
    [0, 0, 0, 0, -1, 0, 0, 0, 0],    # M_A -> ∅
    [0, 0, 0, 0, 0, +1, 0, 0, 0],    # M_A -> A + M_A
    [0, +1, 0, -1, 0, +1, 0, 0, 0],  # D'_R -> A + D_R
    [0, 0, 0, 0, 0, 0, +1, 0, 0],    # D_R -> M_R + D_R
    [0, 0, 0, 0, 0, 0, +1, 0, 0],    # D'_R -> D'_R + M_R
    [0, 0, 0, 0, 0, 0, -1, 0, 0],    # M_R -> ∅
    [0, 0, 0, 0, 0, 0, 0, +1, 0],    # M_R -> M_R + R
])

def PropensityFunc(state, params):
    # Extract parameters for reactions
    gamma_C = params['gamma_C']
    delta_A = params['delta_A']
    delta_R = params['delta_R']
    gamma_A = params['gamma_A']
    gamma_R = params['gamma_R']
    theta_A = params['theta_A']
    alpha_A = params['alpha_A']
    alpha_A_prime = params['alpha_A_prime']
    delta_MA = params['delta_MA']
    beta_A = params['beta_A']
    theta_R = params['theta_R']
    alpha_R = params['alpha_R']
    alpha_R_prime = params['alpha_R_prime']
    delta_MR = params['delta_MR']
    beta_R = params['beta_R']

    D_A, D_R, D_A_prime, D_R_prime, M_A, A, M_R, R, C = state

    # Define propensity functions for each reaction
    return [
        gamma_C * A * R,          # r1: A + R -> C
        delta_A * A,              # r2: A -> ∅
        delta_A * C,              # r3: C -> R
        delta_R * R,              # r4: R -> ∅
        gamma_A * D_A * A,        # r5: D_A + A -> D'_A
        gamma_R * D_R * A,        # r6: D_R + A -> D'_R
        theta_A * D_A_prime,      # r7: D'_A -> A + D_A
        alpha_A * D_A,            # r8: D_A -> D_A + M_A
        alpha_A_prime * D_A_prime,# r9: D'_A -> D'_A + M_A
        delta_MA * M_A,           # r10: M_A -> ∅
        beta_A * M_A,             # r11: M_A -> A + M_A
        theta_R * D_R_prime,      # r12: D'_R -> A + D_R
        alpha_R * D_R,            # r13: D_R -> D_R + M_R
        alpha_R_prime * D_R_prime,# r14: D'_R -> D'_R + M_R
        delta_MR * M_R,           # r15: M_R -> ∅
        beta_R * M_R              # r16: M_R -> M_R + R
    ]

def RandExp(a):
    return -np.log(np.random.random()) / a

def RandDist(react, probs):
    # Pick up an event according to the ratio of the propensity
    return np.random.choice(react, p=probs)

# SSA Solver
def SSA_solver(initial_state, StateChangeMat, final_time):
    [m, n] = StateChangeMat.shape
    ReactNum = np.array(range(m))

    state = np.array(initial_state, dtype=float)
    times = [0]  # Initialize time
    states = [state.copy()]  # Initialize state history

    t = 0
    while t < final_time:
        w = PropensityFunc(state, params)
        a = np.sum(w)
        tau = RandExp(a)
        t += tau
        if t > final_time:
            break
        probs = w / a
        which = RandDist(ReactNum, probs)
        state += StateChangeMat[which, :]
        times.append(t)
        states.append(state.copy())

    return np.array(times), np.array(states)

# Run SSA Solver
final_time = 400  # 400 hours
times, states = SSA_solver(initial_state, StateChangeMat, final_time)

# Extract `A` and `R` values for plotting
A_values = states[:, 5]  # Column for A
R_values = states[:, 7]  # Column for R

# Plotting `A` and `R` values over time
plt.figure(figsize=(12, 6))
plt.plot(times, A_values, label='A (Activator)', color='b')
plt.xlabel('Time (hours)')
plt.ylabel('Concentration of A')
plt.title('Concentration of A Over Time')
plt.legend()
plt.grid(True)
plt.savefig('Concentration_of_A_Over_Time.png')
plt.show()

plt.figure(figsize=(12, 6))
plt.plot(times, R_values, label='R (Repressor)', color='r')
plt.xlabel('Time (hours)')
plt.ylabel('Concentration of R')
plt.title('Concentration of R Over Time')
plt.legend()
plt.grid(True)
plt.savefig('Concentration_of_R_Over_Time.png')
plt.show()

#Task 3
# Modify the parameters for Figure 5 (with delta_R = 0.05)
modified_params = [
    50,    # alpha_A
    500,   # alpha_A_prime
    0.01,  # alpha_R
    50,    # alpha_R_prime
    50,    # beta_A
    5,     # beta_R
    10,    # delta_MA
    0.5,   # delta_MR
    1,     # delta_A
    0.05,  # delta_R (changed from 0.2 to 0.05)
    1,     # gamma_A
    1,     # gamma_R
    2,     # gamma_C
    50,    # theta_A
    100    # theta_R
]

# Initial conditions
initial_conditions = [1, 1, 0, 0, 0, 0, 0, 0, 0]

# Define the ODE model as in part (a)
def Genetic_oscillator(t, y, params):
    D_A, D_R, D_A_prime, D_R_prime, M_A, A, M_R, R, C = y
    alpha_A, alpha_A_prime, alpha_R, alpha_R_prime, beta_A, beta_R, delta_MA, delta_MR, delta_A, delta_R, gamma_A, gamma_R, gamma_C, theta_A, theta_R = params

    dD_A_dt = theta_A * D_A_prime - gamma_A * D_A * A
    dD_R_dt = theta_R * D_R_prime - gamma_R * D_R * A
    dD_A_prime_dt = gamma_A * D_A * A - theta_A * D_A_prime
    dD_R_prime_dt = gamma_R * D_R * A - theta_R * D_R_prime
    dM_A_dt = alpha_A_prime * D_A_prime + alpha_A * D_A - delta_MA * M_A
    dA_dt = beta_A * M_A + theta_A * D_A_prime + theta_R * D_R_prime - A * (gamma_A * D_A + gamma_R * D_R + gamma_C * R + delta_A)
    dM_R_dt = alpha_R_prime * D_R_prime + alpha_R * D_R - delta_MR * M_R
    dR_dt = beta_R * M_R - gamma_C * A * R + delta_A * C - delta_R * R
    dC_dt = gamma_C * A * R - delta_A * C

    return [dD_A_dt, dD_R_dt, dD_A_prime_dt, dD_R_prime_dt, dM_A_dt, dA_dt, dM_R_dt, dR_dt, dC_dt]

# Time span (400 hours)
t_span = (0, 400)
t_eval = np.linspace(0, 400, 1000)

# Solve using solve_ivp
sol_deterministic = solve_ivp(Genetic_oscillator, t_span, initial_conditions, args=(modified_params,), t_eval=t_eval)

# Plot deterministic results
plt.plot(sol_deterministic.t, sol_deterministic.y[7], label='Repressor R (Deterministic)', color='blue')
plt.xlabel('Time (hours)')
plt.ylabel('Concentration of R')
plt.legend()
plt.title('Repressor R Over Time (Deterministic Model)')
plt.grid(True)
plt.savefig('Repressor_R_Over_Time_(Deterministic_Model).png')
plt.show()
# Reuse SSA solver from part (b) with modified parameters (delta_R = 0.05)

modified_params = {
    'alpha_A': 50, 'alpha_A_prime': 500, 'alpha_R': 0.01, 'alpha_R_prime': 50,
    'beta_A': 50, 'beta_R': 5, 'delta_MA': 10, 'delta_MR': 0.5, 'delta_A': 1,
    'delta_R': 0.05,  # Modified value
    'gamma_A': 1, 'gamma_R': 1, 'gamma_C': 2, 'theta_A': 50, 'theta_R': 100
}

# Run SSA solver (from part b)
final_time = 400  # 400 hours
times, states = SSA_solver(initial_state, StateChangeMat, final_time)

# Extract `R` values for plotting
R_values_stochastic = states[:, 7]  # Column for R

# Plot stochastic results
plt.plot(times, R_values_stochastic, label='Repressor R (Stochastic)', color='orange')
plt.xlabel('Time (hours)')
plt.ylabel('Concentration of R')
plt.legend()
plt.title('Repressor R Over Time (Stochastic Model)')
plt.grid(True)
plt.savefig('Repressor_R_Over_Time_(Stochastic_Model).png')
plt.show()
# Plot both deterministic and stochastic results together for comparison
plt.plot(sol_deterministic.t, sol_deterministic.y[7], label='Repressor R (Deterministic)', color='blue')
plt.plot(times, R_values_stochastic, label='Repressor R (Stochastic)', color='orange')
plt.xlabel('Time (hours)')
plt.ylabel('Concentration of R')
plt.legend()
plt.title('Comparison of Repressor R in Deterministic and Stochastic Models')
plt.grid(True)
plt.savefig('Comparison_of_Repressor_R_in_Deterministic_and_Stochastic_Models.png')
plt.show()