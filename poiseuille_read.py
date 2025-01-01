import numpy as np
from pathlib import Path
from natsort import natsorted
import matplotlib as mpl
import matplotlib.pyplot as plt
import scienceplots

plt.style.use('science')
mpl.rcParams['lines.linewidth'] = 3
plt.rcParams.update({'font.size': 30})

time = np.loadtxt("time.dat").T

data_path = Path(".")
u_files = data_path.glob("u_*.dat")
u_files = natsorted(u_files, key=lambda x: x.name)

d = 128.0/2.0
y = np.arange(-d + 1/2, d, 1)
nu = 1.0/3.0*(1.0-0.5)
F = 1e-5
u = np.loadtxt(u_files[-1])
u_exact = [F/(2*nu)*(d**2 - y_i**2) for y_i in y]

fig, ax = plt.subplots(nrows = 1, ncols=2, figsize=(14, 6.0))
fig.tight_layout()
ax[0].plot(u[0,:]/(np.sqrt(1.0/3.0)), y, lw=3, label="LBE")
ax[0].plot(u_exact/(np.sqrt(1.0/3.0)), y, lw=3,ls="--", color="red", label="Analytical")
ax[0].legend()
ax[0].set_xlabel("$u / c_s$")
ax[0].set_ylabel("y")
ax[0].set_ylim(-64, 64)
ax[0].set_yticks([-64, 0, 64])
ax[0].tick_params(axis='y', which='minor', left=False, right=False)
ax[0].tick_params(axis='x', which='minor', bottom=False, top=False)
plt.savefig("poiseuille_mc.pdf")