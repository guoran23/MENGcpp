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

# %% 
# # 第二种形式
# q_bar_2 = q0 + q2 * rho**2 + q4 * rho**4
# delta = np.sqrt(4 * q0 * q4 - q2**2)

# # 避免除以0
# if delta == 0:
#     raise ValueError("delta 为零，无法计算 arctan 表达式，请检查 q0, q2, q4 的值")

# term1 = (q2 + 2 * q4 * rho**2) / delta
# term2 = q2 / delta
# psi_2 = (B0 / delta) * (np.arctan(term1) - np.arctan(term2))


# 绘图
plt.figure(figsize=(12, 6))

# q_bar 对 rho
plt.subplot(1, 2, 1)
plt.plot(rho, q_bar_1, label=r'$\bar{q} = q_0 + q_2 \rho^2$', lw=2)
# plt.plot(rho, q_bar_2, label=r'$\bar{q} = q_0 + q_2 \rho^2 + q_4 \rho^4$', lw=2, linestyle='--')
plt.xlabel(r'$\rho$')
plt.ylabel(r'$\bar{q}$')
plt.title(r'$\bar{q}$ vs $\rho$')
plt.legend()
plt.grid(True)

# psi 对 rho
plt.subplot(1, 2, 2)
plt.plot(rho, psi_1, label=r'$\psi_1$', lw=2)
# plt.plot(rho, psi_2, label=r'$\psi_2$', lw=2, linestyle='--')
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

#%% Plot Alfven continuum
# Parameters
nc = -6
m1d = np.arange(10, 12)  # Equivalent to 3:8 in MATLAB
nr = len(rho)

# --- First Plot: q profiles ---
plt.figure()
plt.plot(rho, q_bar_1, label='q')
plt.legend()
plt.xlabel('rho')
plt.ylabel('q')
plt.title('Profile of q vs √ψ_p')
plt.grid(True)

# --- TAE Frequency Plots ---
fig, axs = plt.subplots(1, 2, figsize=(12, 5))
for fnc in range(2):
    ax = axs[fnc]

    if fnc == 0:
        betati = 0.0028295
    else:
        betati = 0.0

    beta_enhance = 1
    wbae2 = (11 / 4) * betati * beta_enhance

    w2d = np.zeros((len(m1d), nr))

    for i, m in enumerate(m1d):
        w2d[i, :] = np.sqrt(((nc * q_bar_1 + m) / q_bar_1)**2 + wbae2)

    for i, m in enumerate(m1d):
        ax.plot(rho, np.abs(w2d[i, :]), label=f'm={m}')

    ax.set_ylim(0, 1.0)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.legend(loc='best', frameon=False)
    ax.set_ylabel('n + m/q')
    ax.set_xlabel('√ψ_p')

    if fnc == 0:
        ax.set_title(f'n={nc} With β_ti')
    else:
        ax.set_title(f'n={nc} Without β_ti')

plt.tight_layout()
plt.show()


# %%
