import numpy as np
import matplotlib.pyplot as plt
import imageio.v2 as imageio
import os
from multiprocessing import Pool
from tqdm import tqdm
import json

interval = 10
frame_dir = 'frames/Ez'
os.makedirs(frame_dir, exist_ok=True)

# Function to generate a single frame
def generate_frame(i):
    data = np.loadtxt(f'Ez/Ez_{i}.txt')
    fig, ax = plt.subplots()
    im = ax.imshow(data, cmap='jet', origin='lower', aspect='auto', vmin=-0.5, vmax=0.5)
    # im = ax.imshow(data, cmap='jet', origin='lower', aspect='auto')
    plt.colorbar(im, ax=ax, label='Ez value')
    plt.title(f'Ez Field t={i}')
    plt.xlabel('X')
    plt.ylabel('Y')
    fname = f'{frame_dir}/Ez_{i}.png'
    plt.savefig(fname)
    plt.close()
    return fname

with open('../test/TEz/control_file.scf', 'r') as f:    
    data = json.load(f)
interval = data['data_capture_interval']
iterations = data['iterations']

# Create frames in parallel
indices = list(range(0, iterations, interval))

with Pool() as pool:
    filenames = pool.map(generate_frame, indices)

# Write to video
with imageio.get_writer('Ez_movie.mp4', fps=10) as writer:
    for fname in tqdm(filenames):
        image = imageio.imread(fname)
        writer.append_data(image)

# Optional: Remove temporary images
for fname in filenames:
    os.remove(fname)
