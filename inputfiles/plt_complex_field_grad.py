import numpy as np
import matplotlib.pyplot as plt
from compare_derivative import compare_derivative

# === Parameter Setup ===
Nrad = 64
Nthe = 128
Nphi = 16
n = 6
m = [-11] 
rc = 0.5  # 中心值，根据你 C++ 中的代码设置（你可以改）
sigma = 0.1

filename = 'data_complex_field_gradient.txt'

# === Read Data ===
with open(filename, 'r') as f:
    header = f.readline()  # skip header
    data = np.loadtxt(f)

# Assign columns
ptrad  = data[:, 0]
ptthe  = data[:, 1]
ptphi  = data[:, 2]
dfdrad_re   = data[:, 3]
dfdrad_im   = data[:, 4]
dfdthe_re    = data[:, 5]
dfdthe_im  = data[:, 6]
dfdphi_re = data[:, 7]
dfdphi_im = data[:, 8]

dfdrad = dfdrad_re + 1j * dfdrad_im
dfdthe = dfdthe_re + 1j * dfdthe_im
dfdphi = dfdphi_re + 1j * dfdphi_im


# %%
# === Plot dfdrad, dfdthe, dfdphi vs ptrad ===
plt.figure(figsize=(12, 9))
# Row 1: dfdrad
plt.subplot(3, 2, 1)  # 左侧: Re
plt.plot(ptrad, dfdrad_re, label='Re(dfdrad)')
plt.xlabel('ptrad')
plt.ylabel('Re(dfdrad)')
plt.title('Real Part of dfdrad')
plt.grid(True)

plt.subplot(3, 2, 2)  # 右侧: Im
plt.plot(ptrad, dfdrad_im, label='Im(dfdrad)', color='orange')
plt.xlabel('ptrad')
plt.ylabel('Im(dfdrad)')
plt.title('Imaginary Part of dfdrad')
plt.grid(True)

# Row 2: dfdthe
plt.subplot(3, 2, 3)
plt.plot(ptrad, dfdthe_re, label='Re(dfdthe)')
plt.xlabel('ptrad')
plt.ylabel('Re(dfdthe)')
plt.title('Real Part of dfdthe')
plt.grid(True)

plt.subplot(3, 2, 4)
plt.plot(ptrad, dfdthe_im, label='Im(dfdthe)', color='orange')
plt.xlabel('ptrad')
plt.ylabel('Im(dfdthe)')
plt.title('Imaginary Part of dfdthe')
plt.grid(True)

# Row 3: dfdphi
plt.subplot(3, 2, 5)
plt.plot(ptrad, dfdphi_re, label='Re(dfdphi)')
plt.xlabel('ptrad')
plt.ylabel('Re(dfdphi)')
plt.title('Real Part of dfdphi')
plt.grid(True)

plt.subplot(3, 2, 6)
plt.plot(ptrad, dfdphi_im, label='Im(dfdphi)', color='orange')
plt.xlabel('ptrad')
plt.ylabel('Im(dfdphi)')
plt.title('Imaginary Part of dfdphi')
plt.grid(True)

plt.tight_layout()
plt.show()


# === Reshape to 3D ===
shape = (Nrad, Nthe, Nphi)
ptrad  = ptrad.reshape(shape,order='F')
ptthe  = ptthe.reshape(shape,order='F')
ptphi  = ptphi.reshape(shape,order='F')
dfdrad = dfdrad.reshape(shape,order='F')
dfdthe = dfdthe.reshape(shape,order='F')
dfdphi = dfdphi.reshape(shape,order='F')

# === Visualize: ff(r, theta) at fixed phi ===
phi_index = Nphi // 2

R = ptrad[:, :, phi_index]
T = ptthe[:, :, phi_index]

# %%
# 创建解析表达式
f_analytic = np.exp(1j * n * ptphi) * np.exp(1j * m[0] * ptthe) * np.exp(-((ptrad - rc) / sigma) ** 2)
dfdrad_exact = f_analytic * (-2 * (ptrad - rc) / sigma**2)
dfdthe_exact = 1j * m[0] * f_analytic
dfdphi_exact = 1j * n * f_analytic
# === 复数误差 ===
error_dfdrad_c = dfdrad - dfdrad_exact
error_dfdthe_c = dfdthe - dfdthe_exact
error_dfdphi_c = dfdphi - dfdphi_exact

# === 绝对误差（模） ===
error_dfdrad = np.abs(error_dfdrad_c)
error_dfdthe = np.abs(error_dfdthe_c)
error_dfdphi = np.abs(error_dfdphi_c)

# === 截面索引 ===
i_theta = 0 #Nthe // 2
i_phi = 0 #Nphi // 2

# === 画模长误差 ===
plt.figure(figsize=(12, 4))
plt.plot(ptrad[:, i_theta, i_phi], error_dfdrad[:, i_theta, i_phi], label='|dfdrad - exact|')
plt.plot(ptrad[:, i_theta, i_phi], error_dfdthe[:, i_theta, i_phi], label='|dfdthe - exact|')
plt.plot(ptrad[:, i_theta, i_phi], error_dfdphi[:, i_theta, i_phi], label='|dfdphi - exact|')
plt.xlabel("r")
plt.ylabel("Abs Error")
plt.yscale("log")
plt.title("Gradient Absolute Error at θ=π, φ=π")
plt.legend()
plt.grid(True)
plt.show()

# === 画实部/虚部误差 ===
plt.figure(figsize=(12, 4))

plt.subplot(1, 2, 1)
plt.plot(ptrad[:, i_theta, i_phi], np.real(error_dfdrad_c[:, i_theta, i_phi]), label='Re(dfdrad - exact)')
plt.plot(ptrad[:, i_theta, i_phi], np.real(error_dfdthe_c[:, i_theta, i_phi]), label='Re(dfdthe - exact)')
plt.plot(ptrad[:, i_theta, i_phi], np.real(error_dfdphi_c[:, i_theta, i_phi]), label='Re(dfdphi - exact)')
plt.xlabel("r")
plt.ylabel("Error (Real part)")
plt.title("Real Part of Gradient Error")
plt.legend()
plt.grid(True)

plt.subplot(1, 2, 2)
plt.plot(ptrad[:, i_theta, i_phi], np.imag(error_dfdrad_c[:, i_theta, i_phi]), label='Im(dfdrad - exact)')
plt.plot(ptrad[:, i_theta, i_phi], np.imag(error_dfdthe_c[:, i_theta, i_phi]), label='Im(dfdthe - exact)')
plt.plot(ptrad[:, i_theta, i_phi], np.imag(error_dfdphi_c[:, i_theta, i_phi]), label='Im(dfdphi - exact)')
plt.xlabel("r")
plt.ylabel("Error (Imag part)")
plt.title("Imaginary Part of Gradient Error")
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()

# === 对比图：数值 vs 解析 ===
plt.figure(figsize=(15, 8))

# 实部对比
plt.subplot(2, 1, 1)
plt.plot(ptrad[:, i_theta, i_phi], np.real(dfdrad[:, i_theta, i_phi]), label='Re(dfdrad num)', linestyle='-')
plt.plot(ptrad[:, i_theta, i_phi], np.real(dfdrad_exact[:, i_theta, i_phi]), label='Re(dfdrad ana)', linestyle='--')

plt.plot(ptrad[:, i_theta, i_phi], np.real(dfdthe[:, i_theta, i_phi]), label='Re(dfdthe num)', linestyle='-')
plt.plot(ptrad[:, i_theta, i_phi], np.real(dfdthe_exact[:, i_theta, i_phi]), label='Re(dfdthe ana)', linestyle='--')

plt.plot(ptrad[:, i_theta, i_phi], np.real(dfdphi[:, i_theta, i_phi]), label='Re(dfdphi num)', linestyle='-')
plt.plot(ptrad[:, i_theta, i_phi], np.real(dfdphi_exact[:, i_theta, i_phi]), label='Re(dfdphi ana)', linestyle='--')

plt.xlabel("r")
plt.ylabel("Real Part")
plt.title("Real Part: Numerical vs Analytic Gradient")
plt.legend()
plt.grid(True)

# 虚部对比
plt.subplot(2, 1, 2)
plt.plot(ptrad[:, i_theta, i_phi], np.imag(dfdrad[:, i_theta, i_phi]), label='Im(dfdrad num)', linestyle='-')
plt.plot(ptrad[:, i_theta, i_phi], np.imag(dfdrad_exact[:, i_theta, i_phi]), label='Im(dfdrad ana)', linestyle='--')

plt.plot(ptrad[:, i_theta, i_phi], np.imag(dfdthe[:, i_theta, i_phi]), label='Im(dfdthe num)', linestyle='-')
plt.plot(ptrad[:, i_theta, i_phi], np.imag(dfdthe_exact[:, i_theta, i_phi]), label='Im(dfdthe ana)', linestyle='--')

plt.plot(ptrad[:, i_theta, i_phi], np.imag(dfdphi[:, i_theta, i_phi]), label='Im(dfdphi num)', linestyle='-')
plt.plot(ptrad[:, i_theta, i_phi], np.imag(dfdphi_exact[:, i_theta, i_phi]), label='Im(dfdphi ana)', linestyle='--')

plt.xlabel("r")
plt.ylabel("Imag Part")
plt.title("Imaginary Part: Numerical vs Analytic Gradient")
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()
