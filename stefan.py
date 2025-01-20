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

# plt.style.use('science')
# mpl.rcParams['lines.linewidth'] = 3
# plt.rcParams.update({'font.size': 30})


def transcendental(x, St):
    return x*np.exp(x*x)*erf(x)-St/np.sqrt(np.pi)

CS2 = 1.0/3.0
NX = 50
NY = 400


data_path = Path(".")
data_files = data_path.glob("data*.h5")
data_files = natsorted(data_files, key=lambda x: x.name)

phi_limit = 0.5
time = []
h = []

for t in range(len(data_files)):
    with h5py.File(data_files[t], 'r') as data_file:
        t_i = data_file["time"][0]
        time.append(t_i)
        rho = np.array(data_file["rho"])
        u = np.array(data_file["u"])
        v = np.array(data_file["v"])
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
tau_g = 0.52
kappa = 1.0/3.0*(tau_g-0.5)
print(kappa)
H = 2048
params, _ = curve_fit(lambda x, a: a*np.sqrt(x), time, h)

sol = root(lambda x: transcendental(x, St), 0.6)
prefactor = 2.0*sol["x"][0]*np.sqrt(kappa)
print(prefactor, params[0])


plt.loglog(time, h)
plt.loglog(time, prefactor*np.sqrt(time), ls="--")


plt.show()