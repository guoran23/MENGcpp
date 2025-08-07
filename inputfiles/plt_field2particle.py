import numpy as np
import matplotlib.pyplot as plt
from compare_derivative import compare_derivative

# === Parameter Setup ===
Nrad = 32
Nthe = 64
Nphi = 64
n = 6
m = [10, 11] # not used in plot

filename = 'data_field2particle.txt'
# filename = 'data_field2particle_apark.txt'

# === Read Data ===
with open(filename, 'r') as f:
    header = f.readline()  # skip header
    data = np.loadtxt(f)

# Assign columns
ptrad  = data[:, 0]
ptthe  = data[:, 1]
ptphi  = data[:, 2]
ptRR   = data[:, 3]
ptZZ   = data[:, 4]
ff     = data[:, 5]
dfdrad = data[:, 6]
dfdthe = data[:, 7]
dfdphi = data[:, 8]
dfdpar = data[:, 9]

# === Plot ff vs ptrad ===
plt.figure(figsize=(12, 4))
plt.subplot(1, 3, 1)
plt.plot(ptrad, ff)
plt.xlabel('ptrad')
plt.ylabel('ff')
plt.title('ff vs ptrad')
plt.grid(True)

# === Reshape to 3D ===
shape = (Nrad, Nthe, Nphi)
ptrad  = ptrad.reshape(shape,order='F')
ptthe  = ptthe.reshape(shape,order='F')
ptphi  = ptphi.reshape(shape,order='F')
ff     = ff.reshape(shape,order='F')
dfdrad = dfdrad.reshape(shape,order='F')
dfdthe = dfdthe.reshape(shape,order='F')
dfdphi = dfdphi.reshape(shape,order='F')

# === Visualize: ff(r, theta) at fixed phi ===
phi_index = Nphi // 2

R = ptrad[:, :, phi_index]
T = ptthe[:, :, phi_index]
FF = ff[:, :, phi_index]

plt.subplot(1, 3, 2)
contour = plt.contourf(R, T, FF, levels=50, cmap='viridis')
plt.colorbar(contour)
plt.xlabel('ptrad')
plt.ylabel('ptthe')
plt.title(f'ff at phi index = {phi_index}')

# === Check if ff ~ cos(n * phi) ===
i_rad = Nrad // 2
i_the = 0  # MATLAB i_the = 1 means index 0 in Python

phi_line = ptphi[i_rad, i_the, :]
ff_line = ff[i_rad, i_the, :]
ff_norm = ff_line / np.max(np.abs(ff_line))
model = np.cos(n * phi_line)

plt.subplot(1, 3, 3)
plt.plot(phi_line, ff_norm, 'o-', label='ff (normalized)')
plt.plot(phi_line, model, '-', label=f'cos({n}·phi)')
plt.xlabel('phi')
plt.ylabel('ff')
plt.title('ff vs phi at fixed (r, theta)')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()

# === Compare derivatives ===
compare_derivative(ff, ptrad, dfdrad, 'r', 'df/dr')
compare_derivative(ff, ptthe, dfdthe, 'theta', 'df/dtheta')
compare_derivative(ff, ptphi, dfdphi, 'phi', 'df/dphi')
