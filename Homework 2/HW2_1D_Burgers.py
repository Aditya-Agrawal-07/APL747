# %%
import numpy as np
import matplotlib.pyplot as plt
import HW2_utils_data as utils
from scipy.stats import qmc

# %%

T = 1               # Total time of integration
dt = 0.001          # Time step
N = int(T/dt)       # No of time points
a = 1               # Length of space
J = 256             # Space discretization points
epsilon = 0.01/(np.random.normal(np.pi, 2*np.pi/4))                # viscosity

x = np.arange(-a,a,2*a/J)           # Discretization points
# u0 = np.sin(np.pi*x)                # Initial condition
u0 = np.sin(np.random.normal(np.pi, 2*0.005*np.pi)*x)
# u0 = utils.grf(x)                   # For random field

ut, t = utils.pde_oned_Galerkin(u0, T, a, N, J, epsilon)
print(ut.shape)
# %%
fig1 = plt.figure(figsize=(14,8), dpi=100)
plt.imshow(ut, aspect='auto', cmap='jet')
plt.xlabel('Time (T)')
plt.ylabel('Spatial dimension (x)')
plt.title('1D Burgers Equation')

# %%
