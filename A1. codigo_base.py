import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

id_var = 1 # MUDAR PARA A VARIÁVEL DE INTERESSE NOS GRÁFICOS


## Constantes
g = 9.81 # m/s²
m = 2.4 # kg
u0 = 23.43245 # m/s
theta0 = 2.559*180/np.pi # rad, considerar theta
phi0 = 3*180/np.pi # rad

Ixx = 0.05598  # kgm²
Iyy = 0.04788  # kgm²
Izz = 0.10352  # kgm²
Ixz = -0.00211 # kgm²
Ixy = -0.0000075 # kgm²
Iyz = -0.0000002 # kgm²

## Funcoes
def dbeta_dp(Yp, Yr, Ybeta, Lp, Lr, Lbeta, Np, Nr, Nbeta): # Mudança de beta pela rolagem, provem de derivacao das eqs de lateral
    denominador = Ybeta + Yr * ((Lr * Nbeta - Lbeta * Nr) / (Lr * Np - Lp * Nr))
    numerador = Yp + Yr * ((Lp * Nbeta - Lbeta * Np) / (Lr * Np - Lp * Nr))
    resultado = -numerador / denominador
    return resultado

def dbeta_dr(Yp, Yr, Ybeta, Lp, Lr, Lbeta, Np, Nr, Nbeta): # Mudança de beta pela guinada, provem de derivacao das eqs de lateral
        numerador = (Np * Lr / Lp) - Nr
        denominador = Nbeta - (Np * Lbeta / Lp)      
        resultado = numerador / denominador
        return resultado

## Determinados do XFLR5
Xu=	-0.17604
Xalpha=	1.3633
Zu=	-1.9892
Zalpha=	-18.566
Zq=	-4.9505
Mu=	0.0010976
Malpha=	-1.1415
Mq=	-0.82869

	
Ybeta=	-0.32304
Yp=	-0.068839
Yr=	0.15623
Lbeta=	-0.3551
Lp=	-1.1237
Lr=	0.18423
Nbeta=	0.13096
Np=	-0.099993
Nr=	-0.045389


## Calculados por proxy
    
db_dp = dbeta_dp(Yp, Yr, Ybeta, Lp, Lr, Lbeta, Np, Nr, Nbeta)
db_dr = dbeta_dr(Yp, Yr, Ybeta, Lp, Lr, Lbeta, Np, Nr, Nbeta)
Xq = 0
Xbeta =  (2.220175-2.2207)/0.25 # dFx/dbeta, do XFLR5
Xp = Xbeta * db_dp # dCx/dbeta * dbeta_dp
Xr = Xbeta * db_dr

Zbeta = (23.30638-23.32382)/0.25
Zp = Zbeta*db_dp
Zr = Zbeta*db_dr

Mbeta = (-0.01431499 +0.000351469)/0.25
Mp = Mbeta * db_dp
Mr = Mbeta * db_dr

## Para o controle, estimados individualmente no XFLR5

# Subscrito δelevador
X_de=	-2.4249
Y_de=	0.0040723
Z_de=	-52.927
L_de=	-0.00034691
M_de=	-14.969
N_de=	-0.00087569


# Subscrito δrudder (leme)
X_dr =	 -96.58
Y_dr = 	-0.0022158
Z_dr =	 -384.27
L_dr = 	0.0012264
M_dr =	 -11.459
N_dr =	 -0.01483


# Subscrito δaileron
X_da=	-2.7749
Y_da=	0.00056597
Z_da=	-128.04
L_da=	-0.00005845
M_da=	-22.489
N_da=	0.000060727

# Matrizes intermediarias

F = np.array([
    [Xu / m, Xalpha / m, Xq / m, -g * np.cos(theta0) * np.cos(phi0), Xbeta / m, Xp / m, Xr / m, 0],
    [Zu / m, Zalpha / m, Zq/m + u0, g * np.sin(theta0) * np.cos(phi0),
     Zbeta / (m), Zp / (m), Zr / (m), g * np.sin(phi0)],
    [Mu, Malpha, Mq, 0, Mbeta, Mp, Mr, 0],
    [0, 0, np.cos(phi0), 0, 0, 0, -np.sin(phi0), 0],
    [0, 0, 0, 0, Ybeta / (m ), Yp / (m ), (Yr) / (m) - u0, g * np.cos(phi0) * np.cos(theta0)],
    [0, 0, 0, 0, Lbeta, Lp, Lr, 0],
    [0, 0, 0, 0, Nbeta, Np, Nr, 0],
    [0, 0, np.tan(theta0) * np.sin(phi0), 0, 0, 1, np.tan(theta0) * np.cos(phi0), 0]
])

E = np.array([
    [1, 0, 0, 0,  0,  0, 0, 0],
    [0, 1, 0, 0,  0,  0, 0, 0],
    [0, 0, Iyy, 0,  -Ixy,  -Iyz, 0, 0],
    [0, 0, 0, 1,  0,  0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, -Ixy, 0,  0, Ixx,  -Ixz, 0],
    [0, 0, Iyz, 0, 0, -Ixz,  Izz, 0],
    [0, 0, 0, 0,  0,   0,  0, 1]
])


G = np.array([
    [X_de/m, X_dr/m, X_da/m],
    [Z_de/(m*u0), Z_dr/(m*u0), Z_da/(m*u0)],
    [M_de, M_dr, M_da],
    [0, 0, 0],
    [Y_de/(m*u0), Y_dr/(m*u0), Y_da/(m*u0)],
    [L_de, L_dr, L_da],
    [N_de, N_dr, N_da],
    [0, 0, 0]
])

# Para o Espaço de Estados

A = np.linalg.inv(E) @ F
B = np.linalg.inv(E) @ G
C = np.eye(len(A))
D = np.zeros((A.shape[0], B.shape[1])) 


## Análise de estabilidade


eigenvalues = np.linalg.eigvals(A)

reais = np.real(eigenvalues)
imaginarias = np.imag(eigenvalues)

plt.figure()
plt.scatter(reais, imaginarias, color='red', marker='o')
plt.axhline(0, color='black', linewidth=1.5, linestyle='--')
plt.axvline(0, color='black', linewidth=1.5, linestyle='--')
plt.title('Autovalores do sistema')
plt.xlabel('Parte real')
plt.ylabel('Parte imaginária')
plt.grid()
plt.show()

## Análise de resposta de frequência

w = np.logspace(-2, 4, 1000)  # (rad/s)

resp_freq = np.zeros((8, 3, len(w)), dtype=complex)
for idx, omega in enumerate(w): # Computar H(jw) para cada frequêcia
    jwI = 1j * omega * np.eye(A.shape[0])  # jωI
    H = C @ np.linalg.inv(jwI - A) @ B                  # (jωI - A)^-1 , cada elemento representa como cada i_esima variavel de estado responde a j_esima entrada
    resp_freq[:, :, idx] = H         # Resposta para a dada frequencia

nomes_x = ['$u$', '$\\alpha$', '$q$', '$\\theta$', '$\\beta$', '$p$', '$r$', '$\phi$']
nomes_u = ["$\delta_e$", "$\delta_r$", "$\delta_a$"]

# id_var = 0  # MUDAR PARA VARIÁVEL DE INTERESSE NOS GRÁFICOS

for id_u in range(resp_freq.shape[1]):
    magnitude = np.abs(resp_freq[id_var, id_u, :])
    phase = np.angle(resp_freq[id_var, id_u, :], deg=True)

    plt.figure()
    plt.subplot(2, 1, 1)
    plt.loglog(w, magnitude)
    plt.title(f'Relação da entrada {nomes_u[id_u]} para a variável {nomes_x[id_var]}')
    plt.xlabel('Frequency (rad/s)')
    plt.ylabel('Magnitude (dB)')
    y_ticks = [10**-3,10**-2, 10**-1, 10**0, 10**1, 10**2]
    plt.yticks(y_ticks)
    plt.grid(True)
    
    plt.subplot(2, 1, 2)
    plt.semilogx(w, phase)
    plt.xlabel('Frequência (rad/s)')
    plt.ylabel('Fase (graus)')
    plt.grid(True)
    y_ticks = [-180,-90, 0, 90, 180]
    plt.yticks(y_ticks)
    plt.tight_layout()
    plt.show()
    
## Análise de resposta temporal
def dx(t, x):
    if t == 0: 
        u = np.zeros(B.shape[1])
        # u[id_u] = -10  # PARA ENTRADA EM IMPULSO
    else:
        u = np.zeros(B.shape[1])
    u[id_u] = 0.006 # PARA ENTRADA EM DEGRAU
    return A @ x + B @ u

u = np.zeros(B.shape[1])
x0 = np.zeros(A.shape[0])
x0[0] = 1
tempo = [0,5] # t0, tf

id_var = 7  # MUDAR PARA VARIÁVEL DE INTERESSE NOS GRÁFICOS

plt.figure()
temp_n = np.linspace(tempo[0], tempo[1])
x_natural =  np.array([scipy.linalg.expm(A * time).dot(x0) for time in temp_n])
plt.plot(temp_n, x_natural[:, id_var])
plt.title(f'Resposta natural da variável {nomes_x[id_var]}')
plt.xlabel('Time (s)')
plt.ylabel(f'{nomes_x[id_var]}')
plt.grid(True)
plt.show()

for id_u in range(len(nomes_u)):
    methods = ['RK45', 'RK23', 'BDF', 'Radau']  # Available methods
    plt.figure()
    for method in methods:
        sol = solve_ivp(dx, tempo, x0, method=method, t_eval=np.linspace(tempo[0], tempo[1], 500))
        t = sol.t
        x = sol.y
        plt.plot(t, x[id_var, :], label=f'Método {method}')
        print(x)
    plt.title(f'Resposta para degrau de {nomes_u[id_u]} para a variável {nomes_x[id_var]}')
    plt.xlabel('Time (s)')
    plt.ylabel(f'{nomes_x[id_var]}')
    plt.legend()
    plt.grid(True)
    plt.show()

