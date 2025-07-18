import numpy as np
import matplotlib.pyplot as plt

# 参数设定
q0 = 1.71
q2 = 0.16
q4 = 0.5
B0 = 3.0

# ρ从0到1
rho = np.linspace(0, 1, 500)

# 第一种形式
q_bar_1 = q0 + q2 * rho**2
psi_1 = (B0 / (2 * q2)) * np.log(1 + (q2 / q0) * rho**2)
psi_1_0 = (B0 / (2 * q2)) * np.log(1  + (q2 / q0) * 0**2)  # ρ=0时的值
psi_1_edge = (B0 / (2 * q2)) * np.log(1 + (q2 / q0) * 1**2)  # ρ=1时的值
psi_1_norm = (psi_1 - psi_1_0) / (psi_1_edge - psi_1_0)  # 归一化
print("psi_1 at ρ=0:", psi_1_0)
print("psi_1 at ρ=1:", psi_1_edge)

# 第二种形式
q_bar_2 = q0 + q2 * rho**2 + q4 * rho**4
delta = np.sqrt(4 * q0 * q4 - q2**2)

# 避免除以0
if delta == 0:
    raise ValueError("delta 为零，无法计算 arctan 表达式，请检查 q0, q2, q4 的值")

term1 = (q2 + 2 * q4 * rho**2) / delta
term2 = q2 / delta
psi_2 = (B0 / delta) * (np.arctan(term1) - np.arctan(term2))

# 绘图
plt.figure(figsize=(12, 6))

# q_bar 对 rho
plt.subplot(1, 2, 1)
plt.plot(rho, q_bar_1, label=r'$\bar{q} = q_0 + q_2 \rho^2$', lw=2)
plt.plot(rho, q_bar_2, label=r'$\bar{q} = q_0 + q_2 \rho^2 + q_4 \rho^4$', lw=2, linestyle='--')
plt.xlabel(r'$\rho$')
plt.ylabel(r'$\bar{q}$')
plt.title(r'$\bar{q}$ vs $\rho$')
plt.legend()
plt.grid(True)

# psi 对 rho
plt.subplot(1, 2, 2)
plt.plot(rho, psi_1, label=r'$\psi_1$', lw=2)
plt.plot(rho, psi_2, label=r'$\psi_2$', lw=2, linestyle='--')
plt.xlabel(r'$\rho$')
plt.ylabel(r'$\psi$')
plt.title(r'$\psi$ vs $\rho$')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()

plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.plot(rho, psi_1_norm, label=r'$\psi_{1,norm}$', lw=2, linestyle='--')
plt.xlabel(r'$\rho$')
plt.ylabel(r'$\psi$ (normalized)')
plt.title(r'Normalized $\psi$ vs $\rho$')
plt.legend()
plt.grid(True)
plt.subplot(1, 2, 2)
plt.plot(rho, q_bar_1, label=r'$\bar{q} = q_0 + q_2 \rho^2$', lw=2)
plt.plot(rho, q0+psi_1_norm, label=r'$\bar{q} = q_0 + \psi_{1,norm}$', lw=2, linestyle='--')
plt.xlabel(r'$\rho$')
plt.ylabel(r'$\bar{q}$')
plt.title(r'$\bar{q}$ vs $q = q0 + \psi_norm$')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
