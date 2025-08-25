import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import matplotlib
import sys
import os
import json

matplotlib.use("Agg")   # headless, no GUI

# --- Field selection ---
if len(sys.argv) < 2:
    print("Usage: python fieldplotter.py <0:Hz | 1:Ex | 2:Ey>")
    sys.exit(1)

field_map = {0: "Hz", 1: "Ex", 2: "Ey"}
field_choice = int(sys.argv[1])
if field_choice not in field_map:
    print("Invalid choice. Use 0 for Hz, 1 for Ex, 2 for Ey.")
    sys.exit(1)

field_name = field_map[field_choice]
print(f"Generating images for {field_name}...")

# --- Load control file ---
with open("control_file.scf", "r") as file:
    data = json.load(file)

Iterations = data["Iterations"]
data_capture_interval = data["data_capture_interval"]

# --- Value range for plots ---
vmin, vmax = -1e-3, 1e-3

# --- Ensure output directory exists ---
frame_dir = f"data/{field_name}_frames"
os.makedirs(frame_dir, exist_ok=True)

# --- Generate and save images ---
for t in tqdm(range(0, Iterations, data_capture_interval), desc="Fieldplotter"):
    arr = np.loadtxt(f"data/{field_name}/{field_name}{t}.txt")

    plt.figure(figsize=(6, 6))
    # plt.imshow(arr, cmap="rainbow", origin="lower")
    plt.imshow(arr, cmap="rainbow", vmin=vmin, vmax=vmax, origin="lower")
    plt.colorbar(label=f"{field_name} field")
    plt.title(f"{field_name} field snapshot (t={t})")
    plt.xlabel("x index")
    plt.ylabel("y index")
    plt.savefig(f"{frame_dir}/{field_name}{t}.jpg")
    plt.close()

print(f"Saved images to {frame_dir}/")
