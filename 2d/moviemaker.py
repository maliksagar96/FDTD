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
    print(" 0 = Hz, 1 = Ex, 2 = Ey")
    sys.exit(1)

field_arg = int(sys.argv[1])
fields = {0: "Hz", 1: "Ex", 2: "Ey"}

if field_arg not in fields:
    print("Invalid argument. Use 0 (Hz), 1 (Ex), or 2 (Ey).")
    sys.exit(1)

field_name = fields[field_arg]
folder = f"data/{field_name}"
video_name = f"{field_name}_animation.mp4"

# --- Parameters ---
vmin, vmax = -1e-3, 1e-3
fps = 10
frame_width, frame_height = 600, 600

# --- Video writer ---
out = cv2.VideoWriter(video_name, cv2.VideoWriter_fourcc(*'mp4v'), fps, (frame_width, frame_height))

# --- Load control file ---
with open("control_file.scf", "r") as file:
    data = json.load(file)

Iterations = data["Iterations"]
data_capture_interval = data["data_capture_interval"]

# --- Generate frames ---
for t in tqdm(range(0, Iterations, data_capture_interval), desc=f"{field_name}_movie"):
    filepath = os.path.join(folder, f"{field_name}{t}.txt")
    if not os.path.exists(filepath):
        print(f"Warning: {filepath} not found, skipping.")
        continue

    arr = np.loadtxt(filepath)

    fig, ax = plt.subplots(figsize=(6,6))
    # im = ax.imshow(arr, cmap="rainbow", origin="lower")
    im = ax.imshow(arr, cmap="rainbow", vmin=vmin, vmax=vmax, origin="lower")
    ax.set_title(f"{field_name} field snapshot (t={t})")
    ax.set_xlabel("x index")
    ax.set_ylabel("y index")
    fig.colorbar(im, ax=ax, label=f"{field_name} field")

    # Convert Matplotlib figure → NumPy array
    fig.canvas.draw()
    img = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
    img = img.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    img = cv2.resize(img, (frame_width, frame_height))

    out.write(cv2.cvtColor(img, cv2.COLOR_RGB2BGR))  # Matplotlib = RGB, OpenCV = BGR
    plt.close(fig)

out.release()
print(f"Video saved as {video_name}")
