import numpy as np
import matplotlib.pyplot as plt

# Example: load or create your 2D sigma matrix
sigma = np.loadtxt("sigma.txt")  # or however you generate it
kappa_x = np.loadtxt("kappa_x.txt")
kappa_y = np.loadtxt("kappa_y.txt")
pml_a_x = np.loadtxt("pml_a_x.txt")
pml_a_y = np.loadtxt("pml_a_y.txt")
pml_b_x = np.loadtxt("pml_b_x.txt")
pml_b_y = np.loadtxt("pml_b_y.txt") 
pml_c_x = np.loadtxt("pml_c_x.txt")
pml_c_y = np.loadtxt("pml_c_y.txt")

# Plot using imshow (similar to imagesc)
plt.imshow(sigma.T, origin='lower', cmap='hot', aspect='auto')
plt.colorbar(label='Sigma')
plt.title('PML Conductivity Profile')
plt.xlabel('i (x-direction)')
plt.ylabel('j (y-direction)')

plt.figure()
plt.plot(kappa_x)
plt.title('PML Kappa_x Profile')
plt.xlabel('i (x-direction)')
plt.ylabel('Kappa_x')

plt.figure()
plt.plot(kappa_y)
plt.title('PML Kappa_y Profile')
plt.xlabel('i (y-direction)')
plt.ylabel('Kappa_y')

plt.figure()
plt.plot(pml_a_x)
plt.title('PML PML_a_x Profile')
plt.xlabel('i (x-direction)')
plt.ylabel('PML_a_x')

plt.figure()
plt.plot(pml_a_y)
plt.title('PML PML_a_y Profile')
plt.xlabel('i (y-direction)')
plt.ylabel('PML_a_y')

plt.figure()
plt.plot(pml_b_x)
plt.title('PML PML_b_x Profile')
plt.xlabel('i (x-direction)')
plt.ylabel('PML_b_x')

plt.figure()
plt.plot(pml_b_y)
plt.title('PML PML_b_y Profile')
plt.xlabel('i (y-direction)')
plt.ylabel('PML_b_y')

plt.figure()
plt.plot(pml_c_x)
plt.title('PML PML_c_x Profile')
plt.xlabel('i (x-direction)')
plt.ylabel('PML_c_x') 

plt.figure()
plt.plot(pml_c_y)
plt.title('PML PML_c_y Profile')
plt.xlabel('i (y-direction)')
plt.ylabel('PML_c_y')

plt.show()
