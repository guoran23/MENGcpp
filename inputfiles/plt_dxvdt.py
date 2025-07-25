import numpy as np
import matplotlib.pyplot as plt

# ---- 读取数据 ----
xv0_data = np.loadtxt("data_xv0_SP_0.txt")     # shape: (N, 5)
dxvdt_data = np.loadtxt("data_dxvdt_SP_0.txt") # shape: (N, 5)
#%%
rad_min = 0.45
rad_max = 0.55

# ---- 拆分字段 ----
xv0_rad    = xv0_data[:, 0]
xv0_theta  = xv0_data[:, 1]
xv0_phi = xv0_data[:,2]

dxvdt_rad  = dxvdt_data[:, 0]
dxvdt_partw = dxvdt_data[:, 4]

mask = (xv0_rad >= rad_min) & (xv0_rad <= rad_max)

rad_selected = xv0_rad[mask]
theta_selected = xv0_theta[mask]
phitor_selected = xv0_phi[mask]
w_selected = dxvdt_partw[mask]

# ---- 画图 ----
plt.figure(figsize=(15, 4))

# 图 1: dxvdt.rad vs. xv0.rad
plt.subplot(1, 3, 1)
plt.plot(xv0_rad, dxvdt_rad, 'o', markersize=2)
plt.xlabel("xv0.rad")
plt.ylabel("dxvdt.rad")
plt.title("dxvdt.rad vs xv0.rad")
plt.grid(True)

# 图 2: dxvdt.partw vs. xv0.rad
plt.subplot(1, 3, 2)
plt.plot(xv0_rad, dxvdt_partw, 'o', markersize=2)
plt.xlabel("xv0.rad")
plt.ylabel("dxvdt.partw")
plt.title("dxvdt.partw vs xv0.rad")
plt.grid(True)

# 图 3: dxvdt.partw vs. xv0.theta
plt.subplot(1, 3, 3)
plt.plot(xv0_theta, dxvdt_partw, 'o', markersize=2)
plt.xlabel("xv0.theta")
plt.ylabel("dxvdt.partw")
plt.title("dxvdt.partw vs xv0.theta")
plt.grid(True)

plt.tight_layout()
plt.show()

# %%
# %matplotlib notebook
# ---- 拆分变量 ----
x = xv0_data[:, 0]  # partrad
y = xv0_data[:, 1]  # parttheta
z = xv0_data[:, 2]  # partphitor

c = dxvdt_data[:, 4]  # partw (颜色)

# ---- 绘制 3D 散点图 ----
def plot_xyzc(x,y,z,c):
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    sc = ax.scatter(x, y, z, c=c, cmap='viridis', marker='.', s=5)
    plt.colorbar(sc, label='dxvdt.partw')

    ax.set_xlabel("partrad")
    ax.set_ylabel("parttheta")
    ax.set_zlabel("partphitor")
    ax.set_title("3D scatter: dxvdt.partw color-coded")

    ax.view_init(elev=0, azim=0)

    plt.tight_layout()
    plt.show()

# %%
plot_xyzc(x,y,z,c)
plot_xyzc(rad_selected,theta_selected,phitor_selected,w_selected)

# %%
