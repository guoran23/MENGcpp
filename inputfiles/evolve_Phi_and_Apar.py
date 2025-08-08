import numpy as np
import matplotlib.pyplot as plt

# 参数设置
omega1 = -0.1 + 0.4j       # omega1(t)
gamma_d = 0.0              # gamma_d
omega_0 = 1.0              # omega_0
dt = 0.01                  # 时间步长
T = 10                     # 总时间
N = int(T / dt)            # 步数

# 初始化数组（复数）
phi = np.zeros(N, dtype=complex)
A = np.zeros(N, dtype=complex)
t = np.arange(N) * dt

# 初始条件
phi[0] = 1.0e-5 + 0j
A[0] = 1.0e-5 + 0j

# 主循环 - Euler 方法
for n in range(N - 1):
    omega1_final = omega1 - 1j * gamma_d
    dphi_dt = -1j * omega_0 * (A[n] - phi[n]) - 2j * omega1_final * phi[n]
    dA_dt = -1j * omega_0 * (phi[n] - A[n])

    # Euler 更新
    phi[n + 1] = phi[n] + dt * dphi_dt
    A[n + 1] = A[n] + dt * dA_dt

# 绘图
plt.figure(figsize=(8, 10))

# 实部
plt.subplot(3, 1, 1)
plt.plot(t, phi.real, 'b', label='Re(phi)')
plt.plot(t, A.real, 'r--', label='Re(A)')
plt.xlabel('Time')
plt.ylabel('Real Part')
plt.legend()
plt.title('Real Parts')

# 虚部
plt.subplot(3, 1, 2)
plt.plot(t, phi.imag, 'b', label='Im(phi)')
plt.plot(t, A.imag, 'r--', label='Im(A)')
plt.xlabel('Time')
plt.ylabel('Imag Part')
plt.legend()
plt.title('Imaginary Parts')

# 绝对值
plt.subplot(3, 1, 3)
plt.plot(t, np.abs(phi), 'b', label='Abs(phi)')
plt.plot(t, np.abs(A), 'r--', label='Abs(A)')
plt.xlabel('Time')
plt.ylabel('Abs')
plt.legend()
plt.title('Absolute Values')

plt.tight_layout()
plt.show()
