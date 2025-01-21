import numpy as np
from pathlib import Path
from natsort import natsorted
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import scienceplots


data_path = Path(".")
data_files = data_path.glob("data_*.h5")
data_files = natsorted(data_files, key=lambda x: x.name)
n_files = len(data_files)

with h5py.File(data_files[-1], 'r') as data_file:
    rho = np.array(data_file["rho"])
    u = np.array(data_file["u"])

print(rho[0,0])

NX = 200
NY = 200
d = NY/2.0
y = np.arange(-d + 1/2, d, 1)
nu = 1.0/3.0*(1.0-0.5)
F = 1e-5/NX
u_exact = [F/(2*nu)*(d**2 - y_i**2) for y_i in y]
plt.plot(y, u[NX//2,:])
plt.plot(y, u_exact, ls="--")
plt.show()

