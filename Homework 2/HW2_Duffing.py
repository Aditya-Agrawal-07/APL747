import numpy as np
import matplotlib.pyplot as plt
import HW2_utils_data as utils

# %%
""" Generating system response """

# The time parameters:
x0 = np.array([0.1, 0])         # Initial condition
dt, t0, T = 0.0001, 0, 5        # Time discretization
tparam = [dt, t0, T]

k = np.random.normal(1000, 100) 
c = np.random.normal(2, 0.5) 
alpha = np.random.normal(100000, 2500)    # Oscillator parameters
sysparam = [k, c, alpha]

xt, t_eval = utils.duffing(x0, tparam, sysparam)

# %%
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams['font.size'] = 22
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

fig1 = plt.figure(figsize=(14,8), dpi=100)
plt.subplot(2,1,1); 
plt.plot(t_eval, xt[0, :], label='Displacement'); plt.ylabel('$x$(t)');
plt.grid(True); plt.legend(loc=1); plt.margins(0)
plt.subplot(2,1,2); 
plt.plot(t_eval, xt[1, :], label='Velocity'); plt.ylabel('$\dot{x}$(t)');
plt.legend(loc=1); plt.xlabel('Time (s)'); plt.grid(True); plt.margins(0)
plt.suptitle('Duffing Oscillator', y=0.95)
