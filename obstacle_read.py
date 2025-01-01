import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from natsort import natsorted


time = np.loadtxt("time.dat").T

data_path = Path(".")
u_files = data_path.glob("u_*.dat")
u_files = natsorted(u_files, key=lambda x: x.name)

u = np.loadtxt(u_files[-1]).T

plt.imshow(u, origin="lower")
plt.colorbar()
plt.show()
