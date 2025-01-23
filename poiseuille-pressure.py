import numpy as np
from pathlib import Path
from natsort import natsorted
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import scienceplots

plt.style.use('science')
mpl.rcParams['lines.linewidth'] = 3
plt.rcParams.update({'font.size': 30})


data_path = Path(".")
data_files = data_path.glob("data_*.h5")
data_files = natsorted(data_files, key=lambda x: x.name)
n_files = len(data_files)

with h5py.File(data_files[-1], 'r') as data_file:
    rho = np.array(data_file["rho"])
    u = np.array(data_file["u"])

print(rho[0,0])

NX = 32
NY = 128
d = NY/2.0
y = np.arange(-d + 1/2, d, 1)
nu = 1.0/3.0*(1.0-0.5)
c_s = np.sqrt(1.0/3.0)
F = 1e-5
u_exact = np.array([F/(2*nu)*(d**2 - y_i**2) for y_i in y])

fig, ax = plt.subplots(figsize=(7, 6.0))

ax.plot(u[NX//2,:]/c_s, y, label="LBE")
ax.plot(u_exact/c_s, y, ls="--", label="Analytical", color="red")
ax.tick_params(axis='y', which='minor', left=False, right=False)
ax.tick_params(axis='x', which='minor', bottom=False, top=False)
ax.set_xlabel("$u/c_s$")
ax.set_ylabel("$y$")
ax.set_ylim(-64, 64)
ax.set_yticks([-64, 0, 64])
ax.legend()

plt.savefig("poiseuille-pressure.pdf")

