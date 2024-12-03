from sympy import symbols, Matrix, cos, sin, tan, Function

# Define time variable
t = symbols('t')

# Define state variables and their derivatives
u, alpha, q, theta, beta, p, r, phi = symbols('u alpha q theta beta p r phi')
u_dot, alpha_dot, q_dot, theta_dot, beta_dot, p_dot, r_dot, phi_dot = symbols('u_dot alpha_dot q_dot theta_dot beta_dot p_dot r_dot phi_dot')

# Define parameters (gravitational acceleration, mass, inertia matrix components)
g, m, I_xx, I_yy, I_zz, I_xy, I_xz, I_yz = symbols('g m I_xx I_yy I_zz I_xy I_xz I_yz')
phi_0, theta_0 = symbols('phi_0 theta_0')

# Define the delta terms for the external forces/moments
X_u, X_beta, X_alpha, X_p, X_q, X_r, X_alpha_dot = symbols("X'_u X'_beta X'_alpha X'_p X'_q X'_r X'_dot{alpha}")
Y_u, Y_beta, Y_alpha, Y_p, Y_q, Y_r, Y_alpha_dot = symbols("Y'_u Y'_beta Y'_alpha Y'_p Y'_q Y'_r Y'_dot{alpha}")
Z_u, Z_beta, Z_alpha, Z_p, Z_q, Z_r, Z_alpha_dot = symbols("Z'_u Z'_beta Z'_alpha Z'_p Z'_q Z'_r Z'_dot{alpha}")
L_u, L_beta, L_alpha, L_p, L_q, L_r, L_alpha_dot = symbols("L'_u L'_beta L'_alpha L'_p L'_q L'_r L'_dot{alpha}")
M_u, M_beta, M_alpha, M_p, M_q, M_r, M_alpha_dot = symbols("M'_u M'_beta M'_alpha M'_p M'_q M'_r M'_dot{alpha}")
N_u, N_beta, N_alpha, N_p, N_q, N_r, N_alpha_dot = symbols("N'_u N'_beta N'_alpha N'_p N'_q N'_r N'_dot{alpha}")

# Define the ΔX, ΔY, ΔZ, L, M, N terms based on inputs
Delta_X = X_u*u + X_beta*beta + X_alpha*alpha + X_p*p + X_q*q + X_r*r + X_alpha_dot*alpha_dot
Delta_Y = Y_u*u + Y_beta*beta + Y_alpha*alpha + Y_p*p + Y_q*q + Y_r*r + Y_alpha_dot*alpha_dot
Delta_Z = Z_u*u + Z_beta*beta + Z_alpha*alpha + Z_p*p + Z_q*q + Z_r*r + Z_alpha_dot*alpha_dot
L = L_u*u + L_beta*beta + L_alpha*alpha + L_p*p + L_q*q + L_r*r + L_alpha_dot*alpha_dot
M = M_u*u + M_beta*beta + M_alpha*alpha + M_p*p + M_q*q + M_r*r + M_alpha_dot*alpha_dot
N = N_u*u + N_beta*beta + N_alpha*alpha + N_p*p + N_q*q + N_r*r + N_alpha_dot*alpha_dot

# Equations of motion
u_dot_eq = (m*g*cos(phi_0)*cos(theta_0)*theta + Delta_X)/m
beta_dot_eq = (2*m*g*sin(phi_0) - m*g*cos(phi_0)*cos(theta_0) + Delta_Y)/m
alpha_dot_eq = (-2*m*g*cos(phi_0)*cos(theta_0) + Delta_Z)/m

# Inertia matrix
inertia_matrix = Matrix([[I_xx, -I_xy, -I_xz], [-I_xy, I_yy, -I_yz], [-I_xz, -I_yz, I_zz]])
torques = Matrix([L, M, N])
angular_accelerations = inertia_matrix.inv() * torques

# Rotation kinematics
phi_dot_eq = p + tan(theta_0)*(q*sin(phi_0) + r*cos(phi_0))
theta_dot_eq = q*cos(phi_0) - r*sin(phi_0)

# State space form setup: state variables [u, alpha, q, theta, beta, p, r, phi]
state_vars = Matrix([u, alpha, q, theta, beta, p, r, phi])
state_dot_eqs = Matrix([u_dot_eq, alpha_dot_eq, q_dot, theta_dot_eq, beta_dot_eq, p_dot, r_dot, phi_dot_eq])

# Display the result in terms of a state-space form equation
state_space_eq = state_dot_eqs.subs({p_dot: angular_accelerations[0], q_dot: angular_accelerations[1], r_dot: angular_accelerations[2]})

