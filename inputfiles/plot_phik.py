import numpy as np
import matplotlib.pyplot as plt

# === CONFIGURATION ===
filename = "data_phik_initial.txt"
filename_apark = "data_apark_initial.txt"

lenntor = 1          # number of toroidal modes (adjust as needed)
nradfem = 34         # number of radial FEM nodes (adjust as needed)
nthefem = 64        # number of poloidal FEM nodes (adjust as needed)
radmin = 0.1         # min radius
radmax = 0.9         # max radius
themin = 0.0         # min theta (usually 0)
themax = 2 * np.pi   # max theta

m = [10,11]
m_amp= [1.0, 0.9] # Amplitude for each m mode
rc = [0.48,0.5] # Center radius for each m mode
rwidth = [0.1,0.1] # Width of Gaussian for each m mode

dtheta = (themax - themin) / nthefem  # poloidal step size

# === LOAD DATA ===
data = np.loadtxt(filename)
apark_data = np.loadtxt(filename_apark)

assert data.size == apark_data.size == lenntor * nradfem * nthefem, "Data size mismatch!"
# === SELECT A TOROIDAL MODE ===
itor = 0
offset = itor * nradfem * nthefem
phik_slice = data[offset:offset + nradfem * nthefem]
phik_2d = phik_slice.reshape((nthefem, nradfem))
apark_2d = apark_data[offset:offset + nradfem * nthefem].reshape((nthefem, nradfem))

# === CREATE R, THETA GRIDS ===
r = np.linspace(radmin, radmax, nradfem)
theta = np.linspace(themin, themax-dtheta, nthefem)
R, Theta = np.meshgrid(r, theta)

#%% === PLOT PHIK & APARK SIDE-BY-SIDE ===
fig, axs = plt.subplots(1, 2, figsize=(7, 10), constrained_layout=True)
# --- Plot PHIK ---
pcm1 = axs[0].pcolormesh(R, Theta, phik_2d, shading='auto', cmap='viridis')
axs[0].set_title(f'Re[phik], itor = {itor}')
axs[0].set_xlabel("Radius r")
axs[0].set_ylabel("Poloidal angle θ (rad)")
fig.colorbar(pcm1, ax=axs[0], orientation='vertical', label='Re[phik]')

# --- Plot APARK ---
pcm2 = axs[1].pcolormesh(R, Theta, apark_2d, shading='auto', cmap='plasma')
axs[1].set_title(f'Re[apark], itor = {itor}')
axs[1].set_xlabel("Radius r")
axs[1].set_ylabel("Poloidal angle θ (rad)")
fig.colorbar(pcm2, ax=axs[1], orientation='vertical', label='Re[apark]')

plt.suptitle('Comparison of phik and apark in Poloidal Plane', fontsize=14)
plt.show()

#%%Plot phik vs phik_2d_analytical
phik_2d_analytical = np.zeros_like(phik_2d, dtype=np.complex_)
for im, m_val in enumerate(m):
    phik_2d_analytical += m_amp[im] * np.exp(1j * m_val * Theta) * np.exp(-((R - rc[im]) / rwidth[im])**2)

# === PHIK VS R + 高斯理论函数 ===
# 选项1：固定一个 theta（如 theta = 0）
ith = 0
phik_fixed_theta = phik_2d[ith, :]

# 选项2：对所有 theta 平均
phik_avg_theta = np.mean(phik_2d, axis=0)
# phik_analytical at ith
gauss_theory = np.real(phik_2d_analytical[ith,:])  #amp * np.exp(-((r - 0.5) / 0.1)**2)

# === 画图 ===
plt.figure(figsize=(8, 5))
plt.plot(r, phik_fixed_theta.real, 'o-', label=f'phik_FEM at θ = {theta[ith]:.2f} rad')
# plt.plot(r, phik_avg_theta.real, 's--', label='phik averaged over θ')
plt.plot(r, gauss_theory, 'k-*', label='Analytical')
plt.xlabel('Radius r')
plt.ylabel('Re[phik] / f(r)')
plt.title('Comparison: phik(r) vs Analytical Gaussian')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
# %%
# plot phik at r=0.5
ith = np.argmin(np.abs(r - 0.5))
phik_at_r05 = phik_2d[:, ith]
phik_at_r05_ana = np.real(phik_2d_analytical[:, ith])
# === PLOTTING PHIK AT R=0.5 ===
plt.figure(figsize=(8, 5))
plt.plot(theta, phik_at_r05.real, 'o-', label='phik at r=0.5')
plt.plot(theta, phik_at_r05_ana.real, 'k-', label='Theoretical phik at r=0.5')
plt.xlabel('Poloidal angle θ (rad)')
plt.ylabel('Re[phik]')
plt.title(f'phik at r = 0.5, m={m}')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# %%
# === R, THETA GRIDS ===
X = R * np.cos(Theta)
Y = R * np.sin(Theta)

# === SETUP SUBPLOTS ===
fig, axs = plt.subplots(lenntor, 2, figsize=(7, 5 * lenntor), constrained_layout=True)
if lenntor == 1:
    axs = np.reshape(axs, (1, 2))
# === PLOT EACH TOROIDAL MODE ===
for itor in range(lenntor):
    offset = itor * nradfem * nthefem
    phik_slice = data[offset:offset + nradfem * nthefem]
    phik_2d = phik_slice.reshape((nthefem, nradfem))
    apark_slice = apark_data[offset:offset + nradfem * nthefem]
    apark_2d = apark_slice.reshape((nthefem, nradfem))

    # --- Plot PHIK ---
    ax = axs[itor][0]
    pcm = ax.pcolormesh(X, Y, phik_2d, shading='auto', cmap='viridis')
    fig.colorbar(pcm, ax=ax, orientation='vertical', label='Re[phik]')
    ax.set_title(f'Phi k, itor = {itor}')
    ax.set_xlabel('x = r cos(θ)')
    ax.set_ylabel('y = r sin(θ)')
    ax.set_aspect('equal')

    # --- Plot APARK ---
    ax = axs[itor][1]
    pcm = ax.pcolormesh(X, Y, apark_2d, shading='auto', cmap='plasma')
    fig.colorbar(pcm, ax=ax, orientation='vertical', label='Re[A_|| k]')
    ax.set_title(f'A_|| k, itor = {itor}')
    ax.set_xlabel('x = r cos(θ)')
    ax.set_ylabel('y = r sin(θ)')
    ax.set_aspect('equal')

plt.suptitle('Re[phik] in Poloidal Plane for All Toroidal Modes', fontsize=14)
plt.show()