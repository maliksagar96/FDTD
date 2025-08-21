import numpy as np
import matplotlib.pyplot as plt

import matplotlib
matplotlib.use("TkAgg")   # or "Qt5Agg"

vmin = -0.1
vmax = 0.1

# Load field data
for t in np.arange(0, 2000, 50):

	Hz = np.loadtxt(f"data/Hz/Hz{t}.txt")
	plt.figure(figsize=(6,6))
	plt.imshow(Hz, cmap="rainbow",origin="lower")
	plt.colorbar(label="Hz field")
	plt.title("Hz field snapshot")
	plt.xlabel("x index")
	plt.ylabel("y index")
	plt.savefig(f"data/Hz_frames/Hz{t}.jpg")
# plt.show()
