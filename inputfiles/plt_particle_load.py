import numpy as np  
import matplotlib.pyplot as plt
import glob
import os

# 粒子数据中的列名
VAR_NAMES = [
    "partrad", "parttheta", "partphitor", "partvpar",
    "partmu", "partw", "partfog", "R", "Z"
]

def load_particle_file(filename):
    """加载一个 data_particle 文件，返回 Nx9 的数组"""
    data = np.loadtxt(filename)
    return data

def plot_all_variables_multi(data_list, spid):
    """多个 rank 的变量随粒子编号变化图，用不同颜色区分"""
    fig, axs = plt.subplots(3, 3, figsize=(14, 10))
    fig.suptitle(f"All Particle Variables - Species {spid}", fontsize=16)

    for i in range(9):  # 9 个变量
        ax = axs[i // 3][i % 3]
        for data, color, label in data_list:
            ax.plot(data[:, i], '.', markersize=2, color=color, label=label if i == 0 else None)
        ax.set_title(VAR_NAMES[i])
        ax.set_xlabel("Particle Index")
        ax.set_ylabel(VAR_NAMES[i])
        ax.grid(True)
        if i == 0:
            ax.legend(title="MPI Rank")

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()

def plot_phase_space_projections_multi(data_list, spid):
    """多个 rank 的数据绘在一张图上，用不同颜色区分"""
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle(f"Phase Space Projections - Species {spid}", fontsize=18)

    for data, color, label in data_list:
        rad      = data[:, 0]
        theta    = data[:, 1]
        phitor   = data[:, 2]
        vpar     = data[:, 3]
        mu       = data[:, 4]
        w        = data[:, 5]
        fog      = data[:, 6]
        R        = data[:, 7]
        Z        = data[:, 8]

        axs[0, 0].scatter(R, Z, s=2, alpha=0.5, c=color, label=label)
        axs[0, 1].scatter(vpar, mu, s=2, alpha=0.5, c=color)
        axs[0, 2].scatter(mu, w, s=2, alpha=0.5, c=color)
        axs[1, 0].scatter(rad, vpar, s=2, alpha=0.5, c=color)
        axs[1, 1].scatter(theta, phitor, s=2, alpha=0.5, c=color)
        axs[1, 2].scatter(rad, fog, s=2, alpha=0.5, c=color)

    # 设置子图标签
    titles = [
        "R vs Z", "vpar vs mu", "mu vs w",
        "rad vs vpar", "theta vs phitor", "f/g vs r"
    ]
    xlabels = ["R", "vpar", "mu", "rad", "theta", "r"]
    ylabels = ["Z", "mu", "w", "vpar", "phitor", "f/g"]

    for ax, title, xlabel, ylabel in zip(axs.flat, titles, xlabels, ylabels):
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(True)

    axs[0, 0].axis("equal")
    axs[0, 0].legend(title="MPI Rank")

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()

def main():
    # 匹配所有 RANK 文件
    files = sorted(glob.glob("data_particle_SP000_RANK???.txt"))
    if not files:
        print("❌ 没有找到 data_particle_SP000_RANK???.txt 文件")
        return

    # 颜色表，支持超过 8 个 rank 可扩展
    colors = [
        'red', 'blue', 'green', 'orange', 'purple', 'brown',
        'pink', 'cyan', 'olive', 'gray', 'teal', 'magenta'
    ]

    data_list = []
    all_data = []

    for i, file in enumerate(files):
        data = load_particle_file(file)
        color = colors[i % len(colors)]
        rank_str = os.path.basename(file)[-7:-4]  # 提取 "000"、"001" 等
        label = f"RANK{rank_str}"

        print(f"Loaded {file}: shape = {data.shape}, color = {color}, label = {label}")
        data_list.append((data, color, label))
        all_data.append(data)

    # 合并所有数据用于 plot_all_variables
    data_combined = np.vstack(all_data)

    # 提取 species id（例如 SP000）
    spid = os.path.basename(files[0])[17:20]

    # 画图
    plot_all_variables_multi(data_list, spid)
    plot_phase_space_projections_multi(data_list, spid)

if __name__ == "__main__":
    main()
