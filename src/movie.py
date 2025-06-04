import numpy as np
import matplotlib.pyplot as plt

Nt = 101  # number of frames (0 to 100)

for t in range(Nt):
    filename = f"Ey_{t}.txt"
    data = np.loadtxt(filename)

    plt.imshow(data, cmap='viridis', origin='lower', aspect='auto')
    plt.colorbar(label='Ey Value')
    plt.title(f"2D FDTD Ey Field at t={t}")
    plt.xlabel("Y index")
    plt.ylabel("X index")
    plt.tight_layout()

    plt.savefig(f"frame_{t}.png")
    plt.clf()  # Clear figure for next iteration
