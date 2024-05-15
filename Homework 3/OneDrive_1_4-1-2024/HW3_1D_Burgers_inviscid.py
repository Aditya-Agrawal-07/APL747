# %%
import numpy as np
import matplotlib.pyplot as plt
import HW3_utils_data as utils

# %%

T = 0.5               # Total time of integration
dt = 0.001          # Time step
N = int(T/dt)       # No of time points
a = 1               # Length of space
J = 256             # Space discretization points

x = np.arange(0,a,a/J)           # Discretization points
u0 = 0.25*np.sin(2*np.pi*x)                # Initial condition

ut, t = utils.pde_oned_Galerkin(u0, T, a, N, J, epsilon=0)

# %%
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 12
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

fig1 = plt.figure(figsize=(8,4), dpi=100)
plt.imshow(ut, aspect='auto', cmap='jet')
plt.xlabel('Time (T)')
plt.ylabel('Spatial dimension (x)')
plt.title('1D Burgers Equation', fontweight='bold')

print("Helo")
# %%
