#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This code is distributed for the course:
    Course code: APL 747
    Course name: Uncertainty Quantification and Propagation 
    Affiliation: Indian Institute of Technology Delhi 
    Semester: Spring 2024

@author: APL747
"""

import numpy as np
import matplotlib.pyplot as plt
import HW5_utils_data as utils

# %%
""" Generating system response """

# The time parameters:
x0 = np.array([1, 0])          # Initial condition
dt, t0, T = 0.001, 0, 1        # Time discretization
tparam = [dt, t0, T]

k = 1000
c = 1
alpha = 1000
sysparam = [k, c, alpha]

xt, t_eval = utils.duffing(x0, tparam, sysparam)

noise_std = 0.05
xt = xt + noise_std * np.std(xt, axis=1, keepdims=True) * np.random.randn(*xt.shape)

# %%
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams['font.size'] = 22
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

fig1, ax = plt.subplots(nrows=2, ncols=2, figsize=(14,8), dpi=100)
plt.subplots_adjust(wspace = 0.3, hspace = 0.3)
ax = ax.flatten()

ax[0].plot(t_eval, xt[0, :]); 
ax[0].set_ylabel('$x$(t)');
ax[0].grid(True, alpha=0.35); 
ax[0].set_xlabel('Time (s)');

ax[1].plot(t_eval, xt[1, :]); 
ax[1].set_ylabel('$\dot{x}$(t)');
ax[1].grid(True, alpha=0.35); 
ax[1].set_xlabel('Time (s)');

ax[2].plot(xt[0, :], xt[1, :]); 
ax[2].set_xlabel('${x}$(t)');
ax[2].set_ylabel('$\dot{x}$(t)');
ax[2].grid(True, alpha=0.35); 

# ax[3].hist(xt[0,:], bins=100, density=True, histtype='step', color='r')
ax[3].hist(xt[1,:], bins=100, density=True, histtype='step', color='b')
ax[3].set_ylabel('Density')
ax[3].set_xlabel('$\dot{x}$(t)')
ax[3].grid(True, alpha=0.35); 

plt.suptitle('Duffing Oscillator', y=0.95)
