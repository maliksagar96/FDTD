import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Needed for 3D plotting
from matplotlib import cm  # Colormap
from tqdm import tqdm


dataDir = "data/"
framesDir = "frames/"

for i in range(0, 1000):
  filename = "Hy_1D_" + str(i)
  data = np.loadtxt(dataDir+filename+".txt")  
  plt.plot(data)
  plt.ylim([-1.2, 1.2])
  plt.savefig(framesDir + filename + ".png")
  plt.clf()

# plt.show()
