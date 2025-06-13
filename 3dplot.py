import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Needed for 3D plotting
from matplotlib import cm  # Colormap
from tqdm import tqdm


Nt = 1001  # number of frames (0 to 100)



for t in tqdm(range(Nt), desc="Processing frames"):

    filename = f"data/Ey_{t}.txt"
    data = np.loadtxt(filename)

    x = np.arange(data.shape[0])
    y = np.arange(data.shape[1])
    X, Y = np.meshgrid(y, x)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X, Y, data, cmap=cm.viridis, vmin=-1, vmax=1)

    ax.set_title(f"3D FDTD Ey Field at t={t}")
    ax.set_xlabel("Y index")
    ax.set_ylabel("X index")
    ax.set_zlabel("Ey Value")
    fig.colorbar(surf, shrink=0.5, aspect=10)

    plt.tight_layout()
    plt.savefig(f"frames/frame_{t}.png")
    plt.clf()  # Clear figure for next iteration
    plt.close(fig)  # Close figure to release memory
