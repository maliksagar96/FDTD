import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
from tqdm import tqdm

# Optional: Create output directory if it doesn't exist
os.makedirs("frames/Ez", exist_ok=True)

for i in tqdm(np.arange(0, 1000,2), desc = "ElectricFieldZ"):
  filename = f"Ez/Ez_{i}.txt"
  saveFile = f"frames/Ez/Ez_{i:04d}.png"  # zero-padded

  Ez = np.loadtxt(filename)
  ny, nx = Ez.shape
  x = np.arange(nx)
  y = np.arange(ny)
  X, Y = np.meshgrid(x, y)

  fig = plt.figure(figsize=(8, 6))
  ax = fig.add_subplot(111, projection='3d')
  # ax.plot_surface(X, Y, Ez, cmap='viridis', linewidth=0, antialiased=False)

  ax.plot_surface(X, Y, Ez, color='lightgray', linewidth=0, antialiased=False)
  ax.plot_wireframe(X, Y, Ez, color='black', linewidth=0.3)


  ax.set_xlabel('X')
  ax.set_ylabel('Y')
  ax.set_zlabel('Ez')
  ax.set_zlim(-1, 1)

  plt.tight_layout()
  plt.savefig(saveFile, dpi=150)
  plt.close(fig)

  # print(f"Saved frame {i}")
