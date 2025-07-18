import matplotlib.pyplot as plt
import numpy as np
import re

def parse_complex_line(line):
    """Parse one line of (real,imag) formatted complex numbers, with NaN filtering"""
    matches = re.findall(r"\(([-+eE0-9\.nanNAN]+),([-+eE0-9\.nanNAN]+)\)", line)
    result = []
    for r, i in matches:
        real = float(r)
        imag = float(i)
        if not np.isfinite(real) or not np.isfinite(imag):
            return None  # 如果有 NaN 或 Inf，整行丢弃
        result.append(complex(real, imag))
    return result

def load_omega1_data(filename):
    all_steps = []
    with open(filename, 'r') as f:
        for line in f:
            omg1s = parse_complex_line(line)
            all_steps.append(omg1s)
    return np.array(all_steps)  # shape: (n_steps, n_modes)

def plot_omega1_evolution(omega1_arr):
    n_steps, n_modes = omega1_arr.shape

    # Plot Imag(omega1) vs time
    plt.figure(figsize=(10, 6))
    for mode in range(n_modes):
        omg1_Im = np.imag(omega1_arr[:, mode])
        plt.plot(range(n_steps), omg1_Im, label=f'Mode {mode}')
    plt.xlabel("Time Step")
    plt.ylabel("Im(omega1)")
    plt.title("Evolution of Mode Im(omega1_n)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    filename = "omega1.txt"  # 改成你的输出文件名
    omega1_arr = load_omega1_data(filename)
    plot_omega1_evolution(omega1_arr)
