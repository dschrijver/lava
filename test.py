import numpy as np
import h5py
import matplotlib.pyplot as plt


with h5py.File("data_2.h5", "r") as data_file:
    rho = np.array(data_file["rho"])

plt.imshow(rho.T, origin="lower")
plt.colorbar()
plt.show()