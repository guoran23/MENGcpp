import matplotlib.pyplot as plt
import numpy as np

# Step 1: Read data from file
filename = "particle_tot_energy.txt"

try:
    energy = np.loadtxt(filename)

    # Step 2: Create time array assuming constant time step (e.g. 1 unit)
    time = np.arange(len(energy))  # or replace with actual time steps if known

    # Step 3: Plot
    plt.figure(figsize=(8, 5))
    plt.plot(time, energy, label='Total Particle Energy', color='blue')
    plt.xlabel("Time Step")
    plt.ylabel("Total Energy")
    plt.title("Particle Total Energy vs Time")
    plt.grid(True, which='both', linestyle='--', alpha=0.6)
    plt.minorticks_on()
    plt.legend()
    plt.tight_layout()
    plt.show()

except OSError:
    print(f"Error: Could not open file '{filename}'")
except Exception as e:
    print("Error while reading or plotting data:", e)
