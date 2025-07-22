import matplotlib.pyplot as plt
import numpy as np

# Step 1: Read data from file
filename = "particle_tot_energy.txt"

try:
    data = np.loadtxt(filename)
    energy = data[:, 0]  # Assuming the first column contains total energy values
    pert_energy = data[:, 1]  # Assuming the second column contains perturbed energy values 

    # Step 2: Create time array assuming constant time step (e.g. 1 unit)
    time = np.arange(len(energy))  # or replace with actual time steps if known

    # Step 3: Plot
    # 创建图和双轴
    fig, ax1 = plt.subplots(figsize=(8, 5))

    # 左轴：总能量
    ax1.plot(time, energy, color='blue', label='Total Energy')
    ax1.set_xlabel("Time Step")
    ax1.set_ylabel("Total Energy", color='blue')
    ax1.tick_params(axis='y', labelcolor='blue')

    # 右轴：扰动能量
    ax2 = ax1.twinx()
    ax2.plot(time, pert_energy, color='red', linestyle='--', label='Perturbed Energy')
    ax2.set_ylabel("Perturbed Energy", color='red')
    ax2.tick_params(axis='y', labelcolor='red')

    # 网格和美化
    ax1.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.title("Particle Energy vs Time")
    fig.tight_layout()
    plt.show()

except OSError:
    print(f"Error: Could not open file '{filename}'")
except Exception as e:
    print("Error while reading or plotting data:", e)
