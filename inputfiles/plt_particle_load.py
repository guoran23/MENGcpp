import numpy as np
import matplotlib.pyplot as plt
import glob
import os

VAR_NAMES = [
    "partrad", "parttheta", "partphitor", "partvpar",
    "partmu", "partw", "partfog", "R", "Z"
]

def load_particle_file(filename):
    """加载一个 data_particle 文件，返回 Nx9 的数组"""
    data = np.loadtxt(filename)
    return data
def plot_all_variables(data, spid):
    """绘制9个变量随粒子编号的变化图"""
    fig, axs = plt.subplots(3, 3, figsize=(14, 10))
    fig.suptitle(f"All Particle Variables - Species {spid}", fontsize=16)

    for i in range(9):
        ax = axs[i // 3][i % 3]
        ax.plot(data[:, i], '.', markersize=2)
        ax.set_title(VAR_NAMES[i])
        ax.set_xlabel("Particle Index")
        ax.set_ylabel(VAR_NAMES[i])
        ax.grid(True)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(f"particle_vars_{spid}.png", dpi=300)
    plt.show()

def plot_phase_space_projections(data, spid):
    """绘制多个相空间子图"""
    # 数据列映射
    rad      = data[:, 0]
    theta    = data[:, 1]
    phitor   = data[:, 2]
    vpar     = data[:, 3]
    mu       = data[:, 4]
    w        = data[:, 5]
    fog      = data[:, 6]
    R        = data[:, 7]
    Z        = data[:, 8]

    # 设置画布
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle(f"Phase Space Projections - Species {spid}", fontsize=18)

    # 1. R vs Z
    axs[0, 0].scatter(R, Z, s=2, alpha=0.5)
    axs[0, 0].set_xlabel("R")
    axs[0, 0].set_ylabel("Z")
    axs[0, 0].set_title("R vs Z")
    axs[0, 0].grid(True)
    axs[0, 0].axis("equal")

    # 2. vpar vs mu
    axs[0, 1].scatter(vpar, mu, s=2, alpha=0.5)
    axs[0, 1].set_xlabel("vpar")
    axs[0, 1].set_ylabel("mu")
    axs[0, 1].set_title("vpar vs mu")
    axs[0, 1].grid(True)

    # 3. mu vs w
    axs[0, 2].scatter(mu, w, s=2, alpha=0.5)
    axs[0, 2].set_xlabel("mu")
    axs[0, 2].set_ylabel("w")
    axs[0, 2].set_title("mu vs w")
    axs[0, 2].grid(True)

    # 4. rad vs vpar
    axs[1, 0].scatter(rad, vpar, s=2, alpha=0.5)
    axs[1, 0].set_xlabel("rad")
    axs[1, 0].set_ylabel("vpar")
    axs[1, 0].set_title("rad vs vpar")
    axs[1, 0].grid(True)

    # 5. theta vs phitor
    axs[1, 1].scatter(theta, phitor, s=2, alpha=0.5)
    axs[1, 1].set_xlabel("theta")
    axs[1, 1].set_ylabel("phitor")
    axs[1, 1].set_title("theta vs phitor")
    axs[1, 1].grid(True)

    # 最后一张子图空着，清除刻度
    axs[1, 2].axis('off')

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(f"phase_space_species_{spid}.png", dpi=300)
    plt.show()

def main():
    files = sorted(glob.glob("data_particle???.txt"))
    if not files:
        print("❌ 没有找到 data_particle???.txt 文件")
        return

    for file in files:
        data = load_particle_file(file)
        spid = os.path.basename(file)[-7:-4]  
        plot_all_variables(data, spid)
        plot_phase_space_projections(data, spid)

if __name__ == "__main__":
    main()
