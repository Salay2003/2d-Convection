import numpy as np
import matplotlib.pyplot as plt

# defining constants
length = 10
nodes = 11
dx = length / (nodes - 1)
dt = 0.01
nu = 0.1
alpha = 0.1
g = 9.81
Tfinal = 10

""" In order to find temerature of each point in the  grid, we have to also update stream function and vorticities."""

Psi = np.zeros((nodes, nodes))
Omega = np.zeros((nodes, nodes))
Temperature = np.zeros((nodes, nodes))

# initial conditions
Temperature[-1, :] = 1

for time in range(0, int(Tfinal/dt)):
    for i in range(1, nodes - 1):
        for j in range(1, nodes - 1):
            
            # updating vorticity
            
            Omega[i, j] = Omega[i, j] + (dt / dx**2 * ((Psi[i - 1, j] + Psi[i + 1, j]) * (Omega[i, j - 1] - Omega[i, j + 1]) / 4 
                            + (Psi[i, j - 1] - Psi[i, j + 1]) * (Omega[i - 1, j] - Omega[i + 1, j]) / 4
                            + nu/2 * (Omega[i, j - 1] + Omega[i, j + 1] + Omega[i - 1, j] + Omega[i + 1, j] - 4 * Omega[i, j]))
                            - g * dt * (Temperature[i-1, j] - Temperature[i+1, j]) /2*dx)
            
            
            #updating stream function
            
            Psi[i, j] = (Omega[i ,j] * dx**2 + Psi[i, j-1] + Psi[i, j+1] + Psi[i-1, j] + Psi[i+1, j]) / 4
            
            #Updating temperature taking u and v with central difference
            
            u = (Psi[i, j + 1] - Psi[i, j - 1]) / (2 * dx)
            v = (Psi[i + 1, j] - Psi[i - 1, j]) / (2 * dx)
            
            Temperature[i,j] = Temperature[i,j] + dt* ( alpha * (Temperature[i, j - 1] + Temperature[i, j + 1] + Temperature[i - 1, j] + Temperature[i + 1, j] - 4 * Temperature[i, j]) / 2*dx- v*(Temperature[i, j + 1] - Temperature[i, j - 1]) / (2 * dx) - u*(Temperature[i + 1, j] - Temperature[i - 1, j]) / (2 * dx))
            
    if time % 100 == 0:
        print("time elapsed " + str(time*dt))
# plotting the results
plt.figure(figsize=(10, 10))
plt.imshow(Temperature, cmap='coolwarm', interpolation='bicubic')
plt.colorbar()
plt.show()
