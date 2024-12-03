import numpy as np
## PONTO DE ESTAB alpha = 3 deg
## DADOS DO XFLR 5
# Parâmetros para as derivadas de estabilidade aerodinâmica longitudinal
CD = 0.01279  # Coeficiente de arrasto
V0 = 12.8  # Velocidade de cruzeiro (m/s)
dCD_dV = -0.0035  # Derivada parcial de CD em relação à velocidade
rho = 1.225  # Densidade do ar (kg/m^3)
S = 0.56  # Área da asa (m^2)
dT_dV = 0  # Derivada parcial da tração em relação à velocidade
CL = 0.24481  # Coeficiente de sustentação
dCD_dalpha = 0.00177  # Derivada parcial de CD em relação ao ângulo de ataque (alpha, deg)
dCL_dV = 0  # Derivada parcial de CL em relação à velocidade
dCL_dalpha = 0.06722  # Derivada parcial de CL em relação ao ângulo de ataque (alpha, deg)
dCm_dV = 0.00001  # Derivada parcial do coeficiente de momento de arfagem em relação à velocidade
dCm_dalpha = -0.01377  # Derivada parcial do coeficiente de momento de arfagem em relação ao ângulo de ataque



# Derivadas de estabilidade aerodinâmica longitudinal
X_u = -2 * CD - V0 * dCD_dV + (2 / (rho * V0 * S)) * dT_dV  # Força axial em função da velocidade
X_alpha = CL - dCD_dalpha  # Força axial em função do ângulo de ataque
X_q = 0 # Força axial em função da taxa de arfagem --> nula, não há cauda no avião

Z_u = -2 * CL - V0 * dCL_dV  # Força normal em função da velocidade
Z_alpha = -CD - dCL_dalpha  # Força normal em função do ângulo de ataque
Z_q = 0 # Força normal em função da taxa de arfagem --> nula, não há cauda no avião

M_u = V0 * dCm_dV  # Momento de arfagem em função da velocidade
M_alpha = dCm_dalpha  # Momento de arfagem em função do ângulo de ataque
M_q = 0 # Momento de arfagem em função da taxa de arfagem --> nula, não há cauda no avião

# Parâmetros para as derivadas de estabilidade lateral-direcional
SB = 0.05582743  # Área projetada da fuselagem lateral (m^2) 
yB = 0.2 # Coeficiente de arrasto lateral da fuselagem, estimado da literatura
SF = 0.03  # Área do leme (m^2)
a1F = 00  # Inclinação da curva de sustentação do leme
VF = 00  # Razão de volume do estabilizador vertical
CL = 00  # Coeficiente de sustentação para os derivados laterais
Gamma = 2/180*np.pi  # Ângulo de enflechamento da asa (em radianos)
cy = [00, 00, 00]  # Comprimentos das cordas dos elementos da asa (m)
y = [00, 00, 00]  # Posições ao longo da envergadura (m)
ay = [00, 00, 00] # dCl/dalpha local
Cdy = [00, 00, 00] # Cd local
Cly = [00, 00, 00] # Cl local
s = 1.8/2 # Metade da envergadura (m)
dCd_dalpha_y = [00,00,00] # dCd/dalpha local

# Derivadas de estabilidade lateral-direcional
Y_beta = (SB / S) * yB - (SF / S) * a1F  # Força lateral em função do ângulo de deslizamento lateral (beta)
Y_p = 0 # Força lateral em função da taxa de arfagem --> nula, dependente de altura do leme
Y_r = VF*a1F # Força lateral em função da taxa de guinada

L_beta = -1/(S*s) * sum(c * a * Gamma * yi for c, a, yi in zip(cy,ay, y))  # Momento de rolamento em função de beta
L_p = -1/(2*S*s**2) * sum( (a+Cd)*c*yi**2 for a,Cd, c, yi in zip(ay, Cdy, cy, y)) # Momento de rolamento em função da taxa de arfagem
L_r = 1/(S*s**2) * sum( Cl*c*yi**2 for Cl, c, yi in zip(Cly, cy, y)) # Momento de rolamento em função da taxa de guinada

N_beta = VF * a1F  # Momento de guinada em função de beta
N_p = -1/(2*S*s**2) * sum( (Cl-dCd)*c*yi**2 for Cl, dCd, c, yi in zip(Cly, dCd_dalpha_y, cy, y)) # Momento de guinada em função da taxa de arfagem
N_r = -1/(S*s**2) * sum( Cd*c*yi**2 for Cd, c, yi in zip(Cdy, cy, y)) # Momento de guinada em função da taxa de guinada
