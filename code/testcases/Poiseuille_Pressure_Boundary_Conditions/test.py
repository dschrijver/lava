from pathlib import Path
from natsort import natsorted
import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scienceplots

plt.style.use('science')
mpl.rcParams['lines.linewidth'] = 3
plt.rcParams.update({'font.size': 30})


data_path = Path("../../")
data_files = data_path.glob("data_*.h5")
data_files = natsorted(data_files, key=lambda x: x.name)

with h5py.File(data_files[-1], 'r') as data_file:
    u = np.array(data_file["u"])

NX = 32
NY = 128
d = NY/2.0
y = np.arange(-d + 1/2, d, 1)
nu = 1.0/3.0*(1.0-0.5)
c_s = np.sqrt(1.0/3.0)
F = 1e-5
u_exact = np.array([F/(2*nu)*(d**2 - y_i**2) for y_i in y])

Q = np.sum(np.sqrt((u_exact-u[NX//2,:])**2))

print("\n\nTotal deviation of numerical result from analytical result:", Q, "\n\n")


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

plt.savefig("result.pdf")
