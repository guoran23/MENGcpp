import numpy as np
import matplotlib.pyplot as plt

# Load data from the text file
data = np.loadtxt('data_fval.txt')
nx=12
ny=16
#data = np.loadtxt('data_fval_intp_dense.txt')
#nx=40
#ny=41

data = data.reshape(ny,nx)

# Create a pseudocolor plot
plt.pcolormesh(data, shading='auto')  # 'auto' is equivalent to 'interp' in MATLAB
plt.colorbar()  # Add a colorbar

# Show the plot
plt.show()