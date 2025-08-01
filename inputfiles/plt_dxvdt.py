import numpy as np
import matplotlib.pyplot as plt

# ---- 读取数据 ----
xv0_data = np.loadtxt("data_xv0_SP_0.txt")     # shape: (Nptot*nrun, 5)
dxvdt_data = np.loadtxt("data_dxvdt_SP_0.txt") # shape: (Nptot*nrun, 5)
#%%
np_tot = 10000
nrun = xv0_data.shape[0] // np_tot
rad_min = 0.45
rad_max = 0.55
plt_run = 3

# ---- 拆分字段 ----
xv0_rad    = xv0_data[:, 0]
xv0_theta = np.mod(xv0_data[:, 1], 2 * np.pi)
xv0_phi   = np.mod(xv0_data[:, 2], 2 * np.pi)
xv0_vpar = xv0_data[:,3]
xv0_w =xv0_data[:,4]

dxvdt_rad  = dxvdt_data[:, 0]
dxvdt_the = dxvdt_data[:, 1]
dxvdt_phi = dxvdt_data[:,2]
dxcdt_vpar = dxvdt_data[:,3]
dxvdt_partw = dxvdt_data[:, 4]

mask = (xv0_rad >= rad_min) & (xv0_rad <= rad_max)

rad_selected = xv0_rad[mask]
theta_selected = xv0_theta[mask]
phitor_selected = xv0_phi[mask]
w_selected = dxvdt_partw[mask]
#
plt.figure(figsize=(15, 4))
rad2d= xv0_rad.reshape(np_tot,nrun)
the2d= xv0_theta.reshape(np_tot,nrun)
phi2d=xv0_phi.reshape(np_tot,nrun)
weight2d=xv0_w.reshape(np_tot,nrun)
draddt_2d=dxvdt_rad.reshape(np_tot,nrun)
dthedt_2d=dxvdt_the.reshape(np_tot,nrun)
plt.subplot(1,3,1)
plt.plot(rad2d[:,plt_run], weight2d[:,plt_run],'o', markersize=2)
plt.xlabel(" rad")
plt.ylabel("weight")
plt.title(f'plt_run={plt_run}')
plt.subplot(1,3,2)
plt.plot( draddt_2d[:,plt_run])
plt.ylabel(" d rad /dt")
plt.title(f'plt_run={plt_run}')
# 图 2: dxvdt.partw vs. xv0.rad
plt.subplot(1, 3, 3)
plt.plot( dthedt_2d[:,plt_run])
plt.ylabel(" d the /dt")
plt.title(f'plt_run={plt_run}')
plt.grid(True)
rad=rad2d[:,plt_run]
weight=weight2d[:,plt_run]
# 设置 rad 的 bin 边界，例如从 0 到 1 分成 20 个区间
bins = np.linspace(0, 1, 21)  # 20 bins between 0 and 1

# 用 np.histogram 对 weight 按 rad 分 bin 求和
weight_sums, bin_edges = np.histogram(rad, bins=bins, weights=weight)
# 计算每个 bin 的中心点用于横坐标
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

plt.figure(figsize=(8, 5))
plt.bar(bin_centers, weight_sums, width=np.diff(bins), align='center', edgecolor='black')
plt.xlabel('rad')
plt.ylabel('Sum of weights')
plt.title('Weight Sum in rad bins')
plt.grid(True)
plt.tight_layout()
plt.show()

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
