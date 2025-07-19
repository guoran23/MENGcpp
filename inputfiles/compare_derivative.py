import numpy as np
import matplotlib.pyplot as plt

def compare_derivative(f, coord, df_analytical, dim_name, title_str):
    """
    Compare numerical and analytical derivatives in one of the 3 dimensions.

    Parameters:
        f (ndarray): 3D scalar field
        coord (ndarray): 3D coordinate field (e.g., ptrad, ptthe, ptphi)
        df_analytical (ndarray): 3D analytical partial derivative
        dim_name (str): 'r', 'theta', or 'phi'
        title_str (str): Title string for the plot
    """

    Nrad, Nthe, Nphi = f.shape
    df_num = np.zeros_like(f)

    if dim_name == 'r':
        for i in range(1, Nrad - 1):
            df_num[i, :, :] = (f[i + 1, :, :] - f[i - 1, :, :]) / (coord[i + 1, :, :] - coord[i - 1, :, :])
        df_num[0, :, :] = (f[1, :, :] - f[0, :, :]) / (coord[1, :, :] - coord[0, :, :])
        df_num[-1, :, :] = (f[-1, :, :] - f[-2, :, :]) / (coord[-1, :, :] - coord[-2, :, :])

        i_the = Nthe // 2
        i_phi = Nphi // 2
        x_line = coord[:, i_the, i_phi]
        df_a = df_analytical[:, i_the, i_phi]
        df_n = df_num[:, i_the, i_phi]

    elif dim_name == 'theta':
        for j in range(1, Nthe - 1):
            df_num[:, j, :] = (f[:, j + 1, :] - f[:, j - 1, :]) / (coord[:, j + 1, :] - coord[:, j - 1, :])
        df_num[:, 0, :] = (f[:, 1, :] - f[:, 0, :]) / (coord[:, 1, :] - coord[:, 0, :])
        df_num[:, -1, :] = (f[:, -1, :] - f[:, -2, :]) / (coord[:, -1, :] - coord[:, -2, :])

        i_rad = Nrad // 2
        i_phi = Nphi // 2
        x_line = coord[i_rad, :, i_phi]
        df_a = df_analytical[i_rad, :, i_phi]
        df_n = df_num[i_rad, :, i_phi]

    elif dim_name == 'phi':
        for k in range(1, Nphi - 1):
            df_num[:, :, k] = (f[:, :, k + 1] - f[:, :, k - 1]) / (coord[:, :, k + 1] - coord[:, :, k - 1])
        df_num[:, :, 0] = (f[:, :, 1] - f[:, :, 0]) / (coord[:, :, 1] - coord[:, :, 0])
        df_num[:, :, -1] = (f[:, :, -1] - f[:, :, -2]) / (coord[:, :, -1] - coord[:, :, -2])

        i_rad = Nrad // 2
        i_the = Nthe // 2
        x_line = coord[i_rad, i_the, :]
        df_a = df_analytical[i_rad, i_the, :]
        df_n = df_num[i_rad, i_the, :]

    else:
        raise ValueError("Invalid dimension name: must be 'r', 'theta', or 'phi'")

    # === Plot ===
    plt.figure()
    plt.plot(x_line, df_a, 'o-', label='df/dx (analytical)')
    plt.plot(x_line, df_n, '--', label='df/dx (numerical)')
    plt.xlabel(dim_name)
    plt.ylabel(title_str)
    plt.legend()
    plt.title(f'{title_str} comparison at center slice')
    plt.grid(True)
    plt.show()
