import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from tqdm import tqdm
import json
import cv2

matplotlib.use("Agg")  # headless mode

# --- Check argument ---
if len(sys.argv) < 2:
  print("Usage: python moviemaker.py <field_index>")
  print(" 0=Ex, 1=Ey, 2=Ez, 3=Hx, 4=Hy, 5=Hz")
  sys.exit(1)

field_arg = int(sys.argv[1])
fields = {
  0: "Ex",
  1: "Ey",
  2: "Ez",
  3: "Hx",
  4: "Hy",
  5: "Hz"
}

if field_arg not in fields:
  print("Invalid argument. Use 0=Ex, 1=Ey, 2=Ez, 3=Hx, 4=Hy, 5=Hz.")
  sys.exit(1)

field_name = fields[field_arg]
folder = f"data/{field_name}"
video_name = f"{field_name}_animation.mp4"

# --- Load control file ---
with open("control_file.json", "r") as file:
  data = json.load(file)

c0 = 2.99792458e8
Iterations = data["Iterations"]
data_capture_interval = data["data_capture_interval"]
frequency = data["frequency"]
cells_per_wavelength = data["cells_per_wavelength"]
wavelength = c0 / frequency
dx = wavelength / cells_per_wavelength
pml_size = data["pml_size"]

vmap = np.array(data["cmap"])  # [vmin, vmax]

# --- Video writer ---
fps = 10
frame_width, frame_height = 600, 600
out = cv2.VideoWriter(video_name, cv2.VideoWriter_fourcc(*'mp4v'), fps, (frame_width, frame_height))

# --- Generate frames ---
for t in tqdm(range(0, Iterations, data_capture_interval), desc=f"{field_name}_movie"):
  filepath = os.path.join(folder, f"{field_name}{t}.txt")
  if not os.path.exists(filepath):
    print(f"Warning: {filepath} not found, skipping.")
    continue

  arr = np.loadtxt(filepath)
  Ny, Nx = arr.shape  # imshow expects arr[y, x]

  fig, ax = plt.subplots(figsize=(6, 6))
  im = ax.imshow(arr, cmap="rainbow", vmin=vmap[0], vmax=vmap[1], origin="lower")

  # PML strips
  ax.axhspan(0, pml_size, color="black", alpha=0.1)                # bottom
  ax.axhspan(Ny - pml_size, Ny, color="black", alpha=0.1)          # top
  ax.axvspan(0, pml_size, color="black", alpha=0.1)                # left
  ax.axvspan(Nx - pml_size, Nx, color="black", alpha=0.1)          # right

  ax.set_title(f"{field_name} snapshot (t={t})")
  ax.set_xlabel("x index")
  ax.set_ylabel("y index")
  fig.colorbar(im, ax=ax, label=f"{field_name} field")

  # Convert Matplotlib figure → OpenCV frame
  fig.canvas.draw()
  img = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
  img = img.reshape(fig.canvas.get_width_height()[::-1] + (3,))
  img = cv2.resize(img, (frame_width, frame_height))

  out.write(cv2.cvtColor(img, cv2.COLOR_RGB2BGR))
  plt.close(fig)

out.release()
print(f"Video saved as {video_name}")
