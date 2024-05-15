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
import HW4_utils_data as utils

# %%
""" Generating system response """

# The time parameters:
x0 = np.array([1, 0])         # Initial condition
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

fig1 = plt.figure(figsize=(14,8), dpi=100)
plt.subplot(2,1,1); 
plt.plot(t_eval, xt[0, :], label='Displacement'); plt.ylabel('$x$(t)');
plt.grid(True); plt.legend(loc=1); plt.margins(0)
plt.subplot(2,1,2); 
plt.plot(t_eval, xt[1, :], label='Velocity'); plt.ylabel('$\dot{x}$(t)');
plt.legend(loc=1); plt.xlabel('Time (s)'); plt.grid(True); plt.margins(0)
plt.suptitle('Duffing Oscillator', y=0.95)
