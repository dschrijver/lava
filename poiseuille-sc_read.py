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
data_files = data_path.glob("data*.h5")
data_files = natsorted(data_files, key=lambda x: x.name)

n_files = len(data_files)

with h5py.File(data_files[-1], 'r') as data_file:
    rho = np.array(data_file["rho"])
    rho_lava = np.array(data_file["rho_lava"])
    rho_air = np.array(data_file["rho_air"])
    u = np.array(data_file["u"])


NY = u[0,:].shape[0]
u_analytical = np.zeros(NY, dtype=np.float64)
y = np.zeros(NY, dtype=np.float64)
F = 1e-7
a = 32.0
L_y = float(NY)/2.0

tau_lava = 1.0
tau_air = 0.625
nu_lava = 1.0/3.0*(tau_lava - 0.5)
nu_air = 1.0/3.0*(tau_air - 0.5)

for j in range(NY):
    y_i = float(j) + 0.5 - L_y
    y[j] = y_i
    if (np.abs(y_i) > a):
        u_analytical[j] = F/(2.0*nu_lava)*(L_y**2 - y_i**2)
    else:
        u_analytical[j] = F/(2.0*nu_lava)*(L_y**2 - a**2) + F/(2.0*nu_air)*(a**2 - y_i**2)

fig, ax = plt.subplots(nrows = 1, ncols=2, figsize=(14, 6.0), sharey=True)
fig.tight_layout()
ax[0].plot(u[0,:]/(np.sqrt(1.0/3.0)), y, lw=3, label="LBE")
ax[0].plot(u_analytical/(np.sqrt(1.0/3.0)), y, lw=3,ls="--", color="red", label="Analytical")
ax[0].legend()
ax[0].set_xlabel("$u / c_s$")
ax[0].set_ylabel("y")
ax[0].set_ylim(-NY//2, NY//2)
ax[0].set_yticks([-NY//2, 0, NY//2])
ax[0].tick_params(axis='y', which='minor', left=False, right=False)
ax[0].tick_params(axis='x', which='minor', bottom=False, top=False)


ax[1].plot(rho_lava[0,:], y, label="$\\rho_{\\text{lava}}$")
ax[1].plot(rho_air[0,:], y, label="$\\rho_{\\text{air}}$")
ax[1].plot(rho[0,:], y, label="$\\rho$")
ax[1].legend()
ax[1].tick_params(axis='y', which='minor', left=False, right=False)
ax[1].tick_params(axis='x', which='minor', bottom=False, top=False)
ax[1].set_xlabel("$\\rho$")
# plt.savefig("poiseuille_mc.pdf")
plt.show()
