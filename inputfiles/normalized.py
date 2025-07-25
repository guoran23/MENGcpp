import numpy as np
import matplotlib.pyplot as plt

# input 
massbk = 1.0 # [mp] for background plasma
Te=1.0 #[kev]
mass = 2.0 # EP mass [mp]
freq = 4.1*100 /2.0/ np.pi #[kHz]
B = 3.0 #[T]
B_ref = 3.0 #[T]
zcharge = 1.0 # EP charge
Rmajor = 10.0 #[m]
q0= 1.7 # safety factor

Tf = 100 #[kev]
n0= 2.0e19 # [m^-3]

#
mproton = 1.6726219e-27 #kg
qe = 1.60217662e-19 #C
#Permeability of free space, mu0, [henry/m].
munort = 1.25663706143591729e-6# [H/m]

#
B_N = 1

v_N = np.sqrt(2.0 * Te *1000* qe / mproton) #m/s
Tf_norm = Tf/ Te
vts = np.sqrt(Tf *1000* qe / (mass*mproton)) #m/s
t_N = 1.0 / v_N
v_A=  (B/(np.sqrt(munort*n0*massbk*mproton)))
omega_A = v_A/Rmajor
print("v_N [m/s] = {:.3e}".format(v_N))
print("v_A [m/s] = {:.3e}".format(v_A))
omega_N = 1.0 /t_N
omega_gy = zcharge * qe * B / (mass * mproton) #rad/s

paux_T_transit = 2.0 * np.pi * q0 * Rmajor/ vts #s
t_A = Rmajor / v_A #s

#
omega_norm = freq * 1000 * 2.0 * np.pi / omega_N
t_mode_norm = 1.0 / freq / 1000 / t_N
omega_gy_norm = omega_gy / omega_N
paux_T_transit_norm = paux_T_transit / t_N
t_A_norm = t_A / t_N
omega_A_norm = omega_A/omega_N
print("omega_A_norm =", omega_A_norm)
print("omega_norm = ", omega_norm)
print("omega_gy_norm = ", omega_gy_norm)
print("t_mode_norm = ", t_mode_norm)
print("paux_T_transit_norm = ", paux_T_transit_norm)
print("t_A_norm = ", t_A_norm)
print("Tf_norm = ", Tf_norm)