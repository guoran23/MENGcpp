import numpy as np

def load_complex_data(filename):
    """
    从文本文件中读取复数数据，每个复数由连续两列（实部和虚部）组成。
    """
    data = np.loadtxt(filename)
    re_parts = data[:, ::2]
    im_parts = data[:, 1::2]
    complex_data = re_parts + 1j * im_parts
    return complex_data

def main():
    filename = "data_Amplitude.txt"
    Amp_complex = load_complex_data(filename)
    n_steps, n_modes = Amp_complex.shape
    
    print("n_steps, n_modes:", Amp_complex.shape)
    print("Initial Amp:", np.abs(Amp_complex[0]))
    omega1_complex = load_complex_data("data_omega1.txt")
    Apar_complex = load_complex_data("data_Amplitude_Apar.txt")

    # 可视化示例
    import matplotlib.pyplot as plt
    
    plt.figure(figsize=(12, 5))
    plt.subplot(1,2,1)
    plt.semilogy(np.abs(Amp_complex[:, 0]), label="|A|")
    plt.xlabel("Time step")
    plt.ylabel("Amplitude |A|")
    plt.legend()
    plt.grid(True)
    plt.subplot(1,2,2)
    plt.plot(np.real(omega1_complex[:, 0]), label="Re(ω_num)")
    plt.plot(np.imag(omega1_complex[:, 0]), label="Im(ω_num)")
    plt.legend()
    plt.grid(True)
    plt.minorticks_on()   
    plt.grid(True, which='minor', linestyle=':', linewidth=0.5, alpha=0.6)
    plt.tight_layout()
    plt.show()
    
    plt.figure(figsize=(12, 5))
    plt.subplot(1,2,1)
    plt.plot(np.real(Amp_complex[:, 0]), label="Re Phi")
    plt.plot(np.real(Apar_complex[:, 0]), label="Re Apar")
    plt.xlabel("Time step")
    plt.ylabel("Amplitude")
    plt.legend()
    plt.grid(True)
    plt.subplot(1,2,2)
    plt.plot(np.imag(Amp_complex[:, 0]), label="Im Phi")
    plt.plot(np.imag(Apar_complex[:, 0]), label="Im Apar")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    
    
    # filename = "data_amp_with_phase.txt"
    # Amp_with_phase = load_complex_data(filename)
    # plt.figure(figsize=(12, 5))
    # plt.subplot(1,2,1)
    # plt.semilogy(np.abs(Amp_with_phase[:, 0]), label="|A with phase e^{-i omega0 t}|")
    # plt.xlabel("Time step")
    # plt.ylabel("|A with phase e^{-i omega0 t}|")
    # plt.legend()
    # plt.grid(True)
    # plt.subplot(1,2,2)
    # plt.plot(np.real(Amp_with_phase[:, 0]), label="A with phase")
    # plt.plot(np.imag(Amp_with_phase[:, 0]), label="A with phase")
    # plt.legend()
    # plt.grid(True)
    # plt.tight_layout()
    # plt.show()
    
    if False:
        dt = 0.1
        A_t = Amp_complex[:, 0]       # shape: (n_steps,)
        omega1_ana = np.empty(len(A_t) - 1, dtype=complex)

        An = A_t[:-1]
        An1 = A_t[1:]

        with np.errstate(divide='ignore', invalid='ignore'):
            omega1_ana = np.log(An1 / An) / (-1j * dt)

        # 可选：把无效项强制设置为 nan（如果需要）
        omega1_ana[np.isinf(omega1_ana)] = np.nan
        # 
        A_rec = np.zeros(n_steps, dtype=complex)
        A_rec[0] = Amp_complex[0, 0] #1.0 + 0j  # A0

        # 用欧拉法迭代计算 A_n
        for n in range(n_steps - 1):
            A_rec[n+1] = A_rec[n] - 1j * omega1_complex[n] * A_rec[n] * dt
            
        plt.figure(figsize=(12, 5))  # 更宽一点
        # -------- 左图：omega1 比较 --------
        plt.subplot(1, 2, 1)
        plt.plot(np.real(Amp_complex[:, 0]), label="Re(A_sim)")
        plt.plot(np.imag(Amp_complex[:, 0]), label="Im(A_sim)")
        plt.plot(A_rec[1:].real, linestyle="--", label='Re(A_rec)')
        plt.plot(A_rec[1:].imag, linestyle="--", label='Im(A_rec)')
        plt.title("Amplitude Evolution")
        plt.xlabel("Time step")
        plt.ylabel("Amplitude A")
        plt.legend()
        plt.grid(True)   
        # -------- 右图：A 的演化比较 --------
        plt.subplot(1, 2, 2)
        plt.plot(np.real(omega1_complex[:, 0]), label="Re(ω_num)")
        plt.plot(np.imag(omega1_complex[:, 0]), label="Im(ω_num)")
        plt.plot(omega1_ana[1:].real, '--', label='Re(ω_ana)')
        plt.plot(omega1_ana[1:].imag, '--', label='Im(ω_ana)')
        plt.title("Instantaneous Frequency ω")
        plt.xlabel("Time step")
        plt.ylabel("ω (complex)")
        plt.legend()
        plt.grid(True)
        
        plt.tight_layout()
        plt.show()

# 只有当直接运行这个脚本时才执行 main()
if __name__ == "__main__":
    main()
