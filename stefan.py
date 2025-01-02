import numpy as np
from pathlib import Path
from natsort import natsorted
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import scienceplots
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit

# plt.style.use('science')
# mpl.rcParams['lines.linewidth'] = 3
# plt.rcParams.update({'font.size': 30})


CS2 = 1.0/3.0
NX = 50
NY = 400


data_path = Path(".")
data_files = data_path.glob("data*.h5")
data_files = natsorted(data_files, key=lambda x: x.name)

phi_limit = 0.001
t_min = 0
t_max = 25000
time = []
h = []

for t in range(len(data_files)):
    with h5py.File(data_files[t], 'r') as data_file:
        t_i = data_file["time"][0]
        if (t_i > t_min) and (t_i < t_max):
            time.append(t_i)
        rho = np.array(data_file["rho"])
        u = np.array(data_file["u"])
        v = np.array(data_file["v"])
        T = np.array(data_file["T"])
        phi = np.array(data_file["phi"])

    if (t_i > t_min) and (t_i < t_max):
        NX, NY = phi.shape
        for j in range(NY):
            if phi[0,j] > phi_limit:
                h.append(0.5+j)
                break


params, _ = curve_fit(lambda x, a: a*np.sqrt(x), time, h)
print(params)

plt.loglog(time, h)
plt.loglog(time, params[0]*np.sqrt(time), ls="--")
# plt.yscale("log")
# plt.xscale("log")


plt.show()