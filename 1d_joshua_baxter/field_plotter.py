import matplotlib
matplotlib.use('Agg')  # Faster non-GUI backend

import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np
import json
import imageio.v2 as imageio
import os
from multiprocessing import Pool

# Load config
with open("control_file.scf", "r") as file:
    json_data = json.load(file)

iterations = json_data["iterations"]
interval = json_data["data_capture_interval"]

filedir = "data/Ex_frames/"
dataDir = "data/Ex/"


os.makedirs(filedir, exist_ok=True)

indices = list(range(0, iterations, interval))

def process_frame(i):
    filename_txt = f"ElectricField_{i}.txt"

    filepath_txt = os.path.join(dataDir, filename_txt)
    data = np.loadtxt(filepath_txt)

    plt.figure()
    # plt.plot(x_axis, data)

    plt.plot(data)
    plt.ylim([-1, 1])
    
    filename_png = f"ElectricField_{i:04d}.png"
    filepath_png = os.path.join(filedir, filename_png)
    plt.savefig(filepath_png)
    plt.close()

    return filepath_png

# Parallel frame generation
with Pool() as pool:
    image_paths = list(tqdm(pool.imap(process_frame, indices), total=len(indices), desc="Electric Field"))

# Create movie
images = [imageio.imread(img_path) for img_path in image_paths]
imageio.mimsave("Ex_field_movie.mp4", images, fps=20)

# Cleanup
for img_file in image_paths:
    os.remove(img_file)
