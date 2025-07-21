import numpy as np
import matplotlib.pyplot as plt

# === CONFIGURATION ===
filename = "data_phik_initial.txt"
lenntor = 1          # number of toroidal modes (adjust as needed)
nradfem = 22         # number of radial FEM nodes (adjust as needed)
nthefem = 64        # number of poloidal FEM nodes (adjust as needed)
radmin = 0.2         # min radius
radmax = 0.8         # max radius
themin = 0.0         # min theta (usually 0)
themax = 2 * np.pi   # max theta
amp = 0.002
m = [10,11]

dtheta = (themax - themin) / nthefem  # poloidal step size

# === LOAD DATA ===
data = np.loadtxt(filename)
assert data.size == lenntor * nradfem * nthefem, "Data size mismatch!"
# === SELECT A TOROIDAL MODE ===
itor = 0
offset = itor * nradfem * nthefem
phik_slice = data[offset:offset + nradfem * nthefem]
phik_2d = phik_slice.reshape((nthefem, nradfem))

# === CREATE R, THETA GRIDS ===
r = np.linspace(radmin, radmax, nradfem)
theta = np.linspace(themin, themax-dtheta, nthefem)
R, Theta = np.meshgrid(r, theta)

# === PLOTTING ===
plt.figure(figsize=(8, 6))
plt.pcolormesh(R, Theta, phik_2d, shading='auto', cmap='viridis')
plt.colorbar(label='Re[phik]')
plt.xlabel("Radius r")
plt.ylabel("Poloidal angle θ (rad)")
plt.title(f"phik (Re) for toroidal mode itor = {itor}")
plt.tight_layout()
plt.show()

# === PHIK VS R + 高斯理论函数 ===

# 选项1：固定一个 theta（如 theta = 0）
ith = 0
phik_fixed_theta = phik_2d[ith, :]

# 选项2：对所有 theta 平均
phik_avg_theta = np.mean(phik_2d, axis=0)

# === 理论高斯函数 ===
gauss_theory = amp * np.exp(-((r - 0.5) / 0.1)**2)

# === 画图 ===
plt.figure(figsize=(8, 5))
plt.plot(r, phik_fixed_theta.real, 'o-', label=f'phik_FEM at θ = {theta[ith]:.2f} rad')
# plt.plot(r, phik_avg_theta.real, 's--', label='phik averaged over θ')
plt.plot(r, gauss_theory, 'k-', label='Theoretical Gaussian')
plt.xlabel('Radius r')
plt.ylabel('Re[phik] / f(r)')
plt.title('Comparison: phik(r) vs Gaussian')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
# %%
# plot phik at r=0.5
ith = np.argmin(np.abs(r - 0.5))
phik_at_r05 = phik_2d[:, ith]
phik_at_r05_ana = np.zeros_like(theta)
for im in m:
    phik_at_r05_ana += amp * np.cos(im * theta)  # Add theoretical contribution for each mode
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
fig, axs = plt.subplots(1, lenntor, figsize=(5 * lenntor, 5), constrained_layout=True)

# === PLOT EACH TOROIDAL MODE ===
for itor in range(lenntor):
    offset = itor * nradfem * nthefem
    phik_slice = data[offset:offset + nradfem * nthefem]
    phik_2d = phik_slice.reshape((nthefem, nradfem))

    ax = axs[itor] if lenntor > 1 else axs
    pcm = ax.pcolormesh(X, Y, phik_2d, shading='auto', cmap='viridis')
    fig.colorbar(pcm, ax=ax, orientation='vertical', label='Re[phik]')
    ax.set_title(f'Toroidal Mode itor = {itor}')
    ax.set_xlabel('x = r cos(θ)')
    ax.set_ylabel('y = r sin(θ)')
    ax.set_aspect('equal')

plt.suptitle('Re[phik] in Poloidal Plane for All Toroidal Modes', fontsize=14)
plt.show()