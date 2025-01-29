import numpy as np
import h5py
import matplotlib.pyplot as plt


with h5py.File("data_100.h5", "r") as data_file:
    rho = np.array(data_file["rho"])
    u = np.array(data_file["u"])

plt.plot(u[0,:])
# plt.plot(rho[0,:])
plt.show()