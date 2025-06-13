import numpy as np
import matplotlib.pyplot as plt

# Load the matrix from the text file
data = np.loadtxt("Ey_output.txt")

# Plot the heatmap
plt.imshow(data, cmap='viridis', origin='lower', aspect='auto')
plt.colorbar(label='Hz Value')
plt.title("2D FDTD Hz Field")
plt.xlabel("Y index")
plt.ylabel("X index")
plt.tight_layout()
plt.show()
