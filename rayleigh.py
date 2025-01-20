import numpy as np
from pathlib import Path
from natsort import natsorted
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import scienceplots
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.style.use('science')
mpl.rcParams['lines.linewidth'] = 3
plt.rcParams.update({'font.size': 30})


CS2 = 1.0/3.0

# Lattice quantities
NX = 100
NY = 50
NT = 100_000
U_lattice = 0.01 # U_lattice < sqrt(1/3) = 0.577

# Dimensional quantities
T_top = 293                 # K
T_bottom = 294              # K
Delta_T = T_bottom - T_top  # K

LX = 2          # m
LY = 1          # m
L = LY          # m
g = 9.81        # m/s^2
alpha = 207e-6  # 1/K

# Dimensionless numbers
Ra = 1000000
Pr = 0.71

# Computed dimensional quantities
dx = LX / NX                            # m
U = np.sqrt(alpha * Delta_T * g * L)    # m/s
nu = np.sqrt(U**2 * L**2 * Pr / Ra)     # m^2/s
kappa = nu / Pr                         # m^2/s

# Conversion factors
C_l = dx
C_u = U / U_lattice
C_t = C_l / C_u
C_nu = C_l * C_u
C_g = C_u**2 / C_l

# Computed dimensionless quantities
nu_lattice = nu / C_nu
kappa_lattice = kappa / C_nu
tau = 1.0/CS2*nu_lattice + 0.5
tau_g = 1.0/CS2*kappa_lattice + 0.5
g_lattice = g / C_g

# Compute total time passed
t = NT * C_t # s
dt = C_t

# Printing computed quantities
print("dt =", dt)
print("tau =", tau, ", tau_g =", tau_g, ", g_lattice =", g_lattice)


data_path = Path(".")
data_files = data_path.glob("data*.h5")
data_files = natsorted(data_files, key=lambda x: x.name)

with h5py.File(data_files[-1], 'r') as data_file:
    rho = np.array(data_file["rho"]).T
    u = np.array(data_file["u"]).T
    v = np.array(data_file["v"]).T
    T = np.array(data_file["T"]).T

x = np.arange(0.5, NX, 1)
y = np.arange(0.5, NY, 1)

fig, ax = plt.subplots(nrows = 1, ncols=2, figsize=(14, 6.0))
fig.tight_layout()

contour = ax[1].contour(x, y, T, levels=25)
ax[1].tick_params(axis='y', which='minor', left=False, right=False)
ax[1].tick_params(axis='x', which='minor', bottom=False, top=False)
ax[1].set_xticks([0, 50, 100])
ax[1].set_yticks([50])
ax[1].set_xlabel("x")
ax[1].set_ylabel("y")
divider = make_axes_locatable(ax[1])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(contour, cax=cax, orientation='vertical', ticks=[293, 294], label="T")


ax[0].streamplot(x, y, u, v)
ax[0].tick_params(axis='y', which='minor', left=False, right=False)
ax[0].tick_params(axis='x', which='minor', bottom=False, top=False)
ax[0].set_xticks([0, 50, 100])
ax[0].set_yticks([50])
ax[0].set_xlabel("x")
ax[0].set_ylabel("y")
ax[0].set_xlim(0, 100)


plt.savefig("Ra5000.pdf")
# plt.show()
