import numpy as np
import matplotlib.pyplot as plt

# Example: load or create your 2D sigma matrix
sigma = np.loadtxt("sigma.txt")  # or however you generate it

# Plot using imshow (similar to imagesc)
plt.imshow(sigma.T, origin='lower', cmap='hot', aspect='auto')
plt.colorbar(label='Sigma')
plt.title('PML Conductivity Profile')
plt.xlabel('i (x-direction)')
plt.ylabel('j (y-direction)')
plt.show()
