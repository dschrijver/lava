from pathlib import Path
from natsort import natsorted
import numpy as np
import h5py
import matplotlib.pyplot as plt


data_path = Path(".")
data_files = data_path.glob("data_*.h5")
data_files = natsorted(data_files, key=lambda x: x.name)

with h5py.File(data_files[-1], "r") as data_file:
    u = np.array(data_file["u"])

NX = 32
NY = 64
d = NY/2.0
y = np.arange(-d + 1/2, d, 1)
tau = 1.0
nu = 1.0/3.0*(tau-0.5)
F = 1e-5
u_exact = [F/(2*nu)*(d**2 - y_i**2) for y_i in y]

plt.plot(u[0,:], y)
plt.plot(u_exact, y, ls="--", color="red")
plt.show()