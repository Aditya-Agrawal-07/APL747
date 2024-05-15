# -*- coding: utf-8 -*-
"""
This code is distributed for the course:
    Course code: APL 747
    Course name: Uncertainty Quantification and Propagation 
    Affiliation: Indian Institute of Technology Delhi 
    Semester: Spring 2024

@author: APL747
"""

# Load the libraries
import numpy as np
import HW5_utils_data as utils
import matplotlib.pyplot as plt
 
# %% Response simulation

# Parameters: Model
a = 0.5
b = 0.025
c = 0.25

d = 0.005
x0= [60, 50]

# Temporal parameters:
dt = 0.1
T = 100
sysparam = [a, b, c, d]
tparam = [dt, 0, T]
xt, t = utils.lotka_volterra(x0, tparam, sysparam) 

noise_std = 0.02
xt = xt + noise_std * np.std(xt, axis=1, keepdims=True) * np.random.randn(*xt.shape)

# %% Plot for training data and response
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams['font.size'] = 20

fig1, ax = plt.subplots(nrows=1, ncols=2, figsize =(12, 6), dpi=100)
plt.subplots_adjust(wspace = 0.25, hspace = 0.3)
ax = ax.flatten()

ax[0].plot(t, xt[0, :])
ax[0].plot(t, xt[1, :])
ax[0].set_xlabel('Time (s)'); 
ax[0].set_ylabel('Population Size')
ax[0].legend(['Prey','Predator'], loc=1)
ax[0].set_title('(a) States')
ax[0].margins(0)
ax[0].grid(alpha=0.3)

ax[1].plot(xt[0, :], xt[1, :])
ax[1].set_xlabel('Prey'); 
ax[1].set_ylabel('Predator')
ax[1].set_title('(b) Phase-portrait')
ax[1].margins(0)
ax[1].grid(alpha=0.3)
