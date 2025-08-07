# import numpy as np 
import matplotlib.pyplot as plt

# Set the file name
spid = 0  # Assuming spid is an integer
sfile = f"data_track{spid}.txt"
print(f"====loading {sfile}====")

# Load data
var = np.loadtxt(sfile)

# Set parameters
ntrack = 10  # Total number of particles
ntime = var.shape[0] // ntrack  # Number of time steps = total rows / 10
nvar = var.shape[1]  # Number of variables (10)

# Check data consistency
assert var.shape[0] % ntrack == 0, "Number of rows is not divisible by ntrack, data format may be incorrect!"

# Parse data and reshape to make time steps explicit
data = {
    "partrad": var[:, 0].reshape(ntime, ntrack),
    "partthe": var[:, 1].reshape(ntime, ntrack),
    "partphi": var[:, 2].reshape(ntime, ntrack),
    "partvpar": var[:, 3].reshape(ntime, ntrack),
    "partmu": var[:, 4].reshape(ntime, ntrack),
    "partw": var[:, 5].reshape(ntime, ntrack),
    "partR": var[:, 6].reshape(ntime, ntrack),
    "partZ": var[:, 7].reshape(ntime, ntrack),
    "partE": var[:, 8].reshape(ntime, ntrack),
    "partPcan": var[:, 9].reshape(ntime, ntrack),
}

# Load boundary data 'rad=1.0'
datafile = "rzbbbs.txt"
dataedge = np.loadtxt(datafile)
var1 = {
    "Rbbbs": dataedge[:, 0],
    "Zbbbs": dataedge[:, 1]
}

# Plotting
ncol, nrow = 3, 3
fig, axes = plt.subplots(nrow, ncol, figsize=(3 * ncol, 3 * nrow))
axes = axes.flatten()  # Convert 2D array to 1D for easier indexing

# Trajectory plot: R-Z
for fic in range(ntrack):
    axes[0].plot(data["partR"][:, fic], data["partZ"][:, fic])
axes[0].plot(var1["Rbbbs"], var1["Zbbbs"], 'r--')
axes[0].set_aspect('equal')
axes[0].grid(True)               # 显示主网格线
axes[0].minorticks_on()   
axes[0].grid(True, which='minor', linestyle=':', linewidth=0.5, alpha=0.6)
axes[0].set_xlabel("R")
axes[0].set_ylabel("Z")

# Trajectory plot: φ-Z
for fic in range(ntrack):
    axes[1].plot(data["partphi"][:, fic], data["partZ"][:, fic])
axes[1].set_xlabel("φ")
axes[1].set_ylabel("Z")

# Trajectory plot: φ-θ
for fic in range(ntrack):
    axes[2].plot(data["partphi"][:, fic], data["partthe"][:, fic])
axes[2].set_xlabel("φ")
axes[2].set_ylabel("θ")

# Velocity component v_par
axes[3].plot(data["partvpar"])
axes[3].set_title(r"$v_{||}$")

# Variable w
axes[4].plot(data["partw"])
axes[4].set_title("w")

# Variable P_phi
pltPcan = data["partPcan"]
axes[5].plot(pltPcan[:,:]-pltPcan[0,:])
axes[5].set_title(r"$P_{\phi}-Pcan_0$")

# Variable E
axes[6].plot(data["partE"])
axes[6].set_title("E")
# Variable rad
axes[7].plot(data["partrad"])
axes[7].set_title("rad")
#
pltpar=2
axes[8].plot(data["partPcan"][:,pltpar])
axes[8].set_title("Pcan")
plt.tight_layout()
plt.show()
