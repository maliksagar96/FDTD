import matplotlib.pyplot as plt  
from tqdm import tqdm
import numpy as np

filedir = "Ex_frames/"
dataDir = "Ex/"

# for i in tqdm(range(0, 1000), desc="Electric Field"):
#   filename = f"ElectricField_{i}.txt"
#   data = np.loadtxt(dataDir + filename)
#   plt.plot(data)
#   plt.ylim([-1.2, 1.2])
#   plt.savefig(filedir + filename.replace(".txt", ".png"))  # Save as .png
#   plt.clf()

filedir = "Hy_frames/"
dataDir = "Hy/"

for i in tqdm(range(0, 1000), desc="Magnetic Field"):
  filename = f"MagneticField_{i}.txt"
  data = np.loadtxt(dataDir + filename)
  plt.plot(data)
  plt.ylim([-1.2, 1.2])
  plt.savefig(filedir + filename.replace(".txt", ".png"))  # Save as .png
  plt.clf()
