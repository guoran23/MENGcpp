import matplotlib.pyplot as plt
import numpy as np
import re

def parse_complex_line(line):
    """Parse one line of (real,imag) formatted complex numbers"""
    matches = re.findall(r"\(([-+eE0-9\.]+),([-+eE0-9\.]+)\)", line)
    return [complex(float(r), float(i)) for r, i in matches]

def load_amplitude_data(filename):
    all_steps = []
    with open(filename, 'r') as f:
        for line in f:
            amps = parse_complex_line(line)
            all_steps.append(amps)
    return np.array(all_steps)  # shape: (n_steps, n_modes)

def plot_amplitude_evolution(amplitude_arr):
    n_steps, n_modes = amplitude_arr.shape

    # Plot magnitude vs time
    plt.figure(figsize=(10, 6))
    for mode in range(n_modes):
        amp_mag = np.abs(amplitude_arr[:, mode])
        plt.plot(range(n_steps), amp_mag, label=f'Mode {mode}')
    plt.xlabel("Time Step")
    plt.ylabel("Amplitude Magnitude")
    plt.title("Evolution of Mode Amplitudes")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # Plot complex trajectories
    plt.figure(figsize=(6, 6))    
    for mode in range(n_modes):
        plt.plot(amplitude_arr[:, mode].real, amplitude_arr[:, mode].imag, label=f'Mode {mode}')
    plt.xlabel("Re")
    plt.ylabel("Im")
    plt.title("Amplitude Complex Trajectories")
    plt.axis("equal")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    filename = "amplitude.txt"  # 改成你的输出文件名
    amplitude_arr = load_amplitude_data(filename)
    plot_amplitude_evolution(amplitude_arr)
