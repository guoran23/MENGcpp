import os
import glob
import numpy as np
import matplotlib.pyplot as plt

def load_profile_data(filename):
    """
    读取单个 species 文件，返回 rad, dens, dlndensdr 三个 numpy 数组。
    """
    data = np.loadtxt(filename)
    rad = data[:, 0]
    dens = data[:, 1]
    dlndensdr = data[:, 2]
    tem = data[:,3]
    return rad, dens, dlndensdr, tem

def plot_profiles(directory="."):
    """
    在指定目录中查找所有 species_*_dens.txt 文件，逐个读取并画图。
    """
    file_pattern = os.path.join(directory, "data_species_*_dens.txt")
    files = sorted(glob.glob(file_pattern))

    if not files:
        print("❌ 没有找到匹配的 profile 文件。")
        return

    for file in files:
        rad, dens, dlndensdr, tem = load_profile_data(file)
        spid = os.path.basename(file).split("_")[1]  # 从文件名中提取 species ID

        plt.figure(figsize=(10, 5))

        # 子图 1: 密度
        plt.subplot(2, 2, 1)
        plt.plot(rad, dens, label=f"Species {spid}")
        plt.xlabel("Radius")
        plt.ylabel("Density")
        plt.title(f"Density Profile - Species {spid}")
        plt.yscale('log')
        plt.grid(True)

        # 子图 2: 对数导数
        plt.subplot(2, 2, 2)
        plt.plot(rad, dlndensdr, label=f"Species {spid}")
        plt.xlabel("Radius")
        plt.ylabel("dln(Density)/dr")
        plt.title(f"dln(Density)/dr - Species {spid}")
        plt.grid(True)
        # 子图 3: 温度
        plt.subplot(2, 2, 3)
        tem = np.loadtxt(file, usecols=(3,))
        plt.plot(rad, tem, label=f"Species {spid}")
        plt.xlabel("Radius")
        plt.ylabel("Temperature")
        plt.title(f"Temperature Profile - Species {spid}")
        plt.grid(True)

        plt.tight_layout()
        # plt.savefig(f"profile_species_{spid}.png")
        plt.show()

if __name__ == "__main__":
    plot_profiles()
