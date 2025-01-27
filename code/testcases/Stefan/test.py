import numpy as np
from pathlib import Path
from natsort import natsorted
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import scienceplots
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit, root
from scipy.special import erf

plt.style.use('science')
mpl.rcParams['lines.linewidth'] = 3
plt.rcParams.update({'font.size': 30})


def transcendental(x, St):
    return x*np.exp(x*x)*erf(x)-St/np.sqrt(np.pi)

CS2 = 1.0/3.0
NX = 4
NY = 2048


data_path = Path("../../")
data_files = data_path.glob("data_*.h5")
data_files = natsorted(data_files, key=lambda x: x.name)

phi_limit = 0.9
time = []
h = []

for t in range(len(data_files)):
    with h5py.File(data_files[t], 'r') as data_file:
        t_i = data_file["time"][0]
        time.append(t_i)
        T = np.array(data_file["T"])
        phi = np.array(data_file["phi"])

    NX, NY = phi.shape
    for j in reversed(range(NY)):
        if phi[0,j] < phi_limit:
            h.append(0.5+j)
            break
    else:
        h.append(0)

time = np.array(time)
h = np.array(h)

St = 1.0
tau_g = 0.50498
kappa = 1.0/3.0*(tau_g-0.5)
params, _ = curve_fit(lambda x, a: a*np.sqrt(x), time, h)

sol = root(lambda x: transcendental(x, St), 0.6)
prefactor = 2.0*sol["x"][0]*np.sqrt(kappa)


fig, ax = plt.subplots(figsize=(7, 6.0))

ax.scatter(time, h)
ax.loglog(time, prefactor*np.sqrt(time), ls="--", label="$y_i(t) = 2 \\lambda \\sqrt{\\kappa t} $", color="red")
# ax.tick_params(axis='y', which='minor', left=False, right=False)
ax.tick_params(axis='x', which='minor', bottom=False, top=False)
ax.set_xlabel("$t$")
ax.set_ylabel("$y_i$")
ax.legend()


plt.savefig("result.pdf")