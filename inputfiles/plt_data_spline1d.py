import numpy as np
import matplotlib.pyplot as plt

# Load data
original_data = np.loadtxt('data_fval.txt')         # Coarse grid data
spline_data = np.loadtxt('data_fspl.txt')  # Spline interpolated data

interpolated_data = np.loadtxt('data_fval_intp_dense.txt')  # Interpolated to finer grid

# Define the x-axis based on grid sizes
nx_coarse = 12
nx_fine = 40

dx_coarse = 1 / nx_coarse  # Coarse grid spacing
x_coarse = np.linspace(0, 1-1/nx_coarse, nx_coarse)
x_fine = np.linspace(0 - dx_coarse, 1 + dx_coarse, nx_fine)
y_ana = np.cos(2 * np.pi * x_fine)  # Analytical function for comparison
# Plotting
plt.figure(figsize=(10, 6))
plt.plot(x_coarse, original_data, 'o-', label='Original Data (nx=12)', color='blue', markersize=6)
plt.plot(x_coarse, spline_data, 's-', label='Spline Interpolated Data (nx=12)', color='orange', markersize=4)
plt.plot(x_fine, interpolated_data, 'x-', label='Interpolated Data (nx=40)', color='red', markersize=5)
plt.plot(x_fine, y_ana, '--', label='Analytical Function', color='green')
# Add plot decorations
plt.xlabel('Normalized Position')
plt.ylabel('Function Value')
plt.title('Comparison of Original and Interpolated Data')
plt.grid(True)
plt.legend()
plt.tight_layout()

# Show the plot
plt.show()
