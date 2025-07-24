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

def load_amplitude_data(filename):
    all_steps = []
    skipped_steps = 0
    expected_length = None
    with open(filename, 'r') as f:
        for lineno, line in enumerate(f, 1):
            amps = parse_complex_line(line)
            if amps is None:
                print(f"⚠️  Skipping line {lineno}: contains NaN or Inf -> {line.strip()}")
                skipped_steps += 1
                continue
            if expected_length is None:
                expected_length = len(amps)
            elif len(amps) != expected_length:
                print(f"⚠️  Skipping line {lineno}: inconsistent mode count (expected {expected_length}, got {len(amps)})")
                skipped_steps += 1
                continue
            all_steps.append(amps)
    print(f"✅ Loaded {len(all_steps)} valid steps with {expected_length} modes, skipped {skipped_steps} bad steps.")
    return np.array(all_steps)

def plot_amplitude_evolution(amplitude_arr):
    n_steps, n_modes = amplitude_arr.shape

    # Plot magnitude over time
    plt.figure(figsize=(10, 6))
    for mode in range(n_modes):
        amp_mag = np.abs(amplitude_arr[:, mode])
        plt.semilogy(range(n_steps), (amp_mag), label=f'Mode {mode}')
    plt.xlabel("Time Step")
    plt.ylabel("Amplitude Magnitude")
    plt.title("Evolution of Mode Amplitudes")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # Plot complex trajectory
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
    filename = "data_Amplitude.txt"  # 修改为你的文件名
    amplitude_arr = load_amplitude_data(filename)
    if amplitude_arr.size > 0:
        plot_amplitude_evolution(amplitude_arr)
    else:
        print("❌ 没有有效的 amplitude 数据可供绘图。")
