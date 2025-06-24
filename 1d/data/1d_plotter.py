import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Needed for 3D plotting
from matplotlib import cm  # Colormap
from tqdm import tqdm


plotEx = 1
plotHy = 0


if(plotEx):
  dataDir = "Ex/"
  framesDir = "Ex_frames/"

  for i in tqdm(range(0, 500), desc = "ElectricField"):
    filename = "ElectricField_" + str(i)
    data = np.loadtxt(dataDir+filename+".txt")  
    plt.plot(data)
    plt.ylim([-1.2, 1.2])
    plt.savefig(framesDir + filename + ".png")
    plt.clf()

# plt.show()


if(plotHy):
  dataDir = "Hy/"
  framesDir = "Hy_frames/"

  for i in tqdm(range(0, 500), desc = "MagneticField"):
    filename = "MagneticField_" + str(i)
    data = np.loadtxt(dataDir+filename+".txt")  
    plt.plot(data)
    plt.ylim([-1.2, 1.2])
    plt.savefig(framesDir + filename + ".png")
    plt.clf()

# plt.show()
