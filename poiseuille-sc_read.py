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

with h5py.File(data_files[-1], 'r') as data_file:
    rho = np.array(data_file["rho"])
    rho_lava = np.array(data_file["rho_lava"])
    rho_air = np.array(data_file["rho_air"])
    u = np.array(data_file["u"])


NY = u[0,:].shape[0]
u_analytical = np.zeros(NY)
u_analytical_2 = np.zeros(NY)
F = 1e-5
a = NY/4
d = NY/2
nu_lava = 1.0/3.0*(1.0-0.5)
nu_air = 1.0/3.0*(0.625-0.5)
print(nu_lava, nu_air)
y = np.arange(-d + 1/2, d, 1)

for j in range(NY):
    if (j < NY/4) or (j >= NY*3/4):
        u_analytical[j] = F/(2*nu_lava)*(d**2 - y[j]**2)
    else:
        u_analytical[j] = F/(2*nu_lava)*(d**2 - a**2) + F/(2*nu_air)*(a**2 - y[j]**2)
    u_analytical_2[j] = F/(2*nu_lava)*(d**2 - y[j]**2)
        

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
