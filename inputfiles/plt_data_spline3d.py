import numpy as np
import matplotlib.pyplot as plt

# Load data from the text file
data = np.loadtxt('data_fspl.txt')
nx=12
ny=16
nz=18
#data = np.loadtxt('data_fval_intp_dense.txt')
#nx=40
#ny=41
#nz=42

data = data.reshape(nz,ny,nx)

# Create a figure with subplots
plt.figure(figsize=(5, 8))  # Adjust figure size if needed

# Loop through the first 4 slices (fic = 1 to 4)
for fic in range(8):
    plt.subplot(4, 2, fic + 1)  # Create a 2x2 grid of subplots
    plt.pcolormesh(data[:, :, fic], shading='auto')  # Pseudocolor plot with interpolation
    plt.colorbar()  # Add a colorbar
    plt.title(f'Slice {fic + 1}')  # Add a title to each subplot

# Adjust layout and display the plot
plt.tight_layout()
plt.show()