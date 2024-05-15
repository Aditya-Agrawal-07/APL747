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

np.random.seed(2)

nx = 32                            # number of points in the x direction
ny = 32                            # number of points in the y direction
xmin, xmax = -1, 1                  # limits in the x direction
ymin, ymax = -1, 1                  # limits in the y direction
dx = (xmax - xmin) / (nx-1)         # grid spacing in the x direction
dy = (ymax - ymin) / (ny-1)         # grid spacing in the y direction

# Create the gridline locations and the mesh grid;
x = np.linspace(xmin, xmax, nx)
y = np.linspace(ymin, ymax, ny)
X, Y = np.meshgrid(x, y, indexing='ij')

# Compute the rhs
b = np.sin(2*np.pi*X)+np.sin(2*np.pi*Y)               # Source field 
b[0, :] = b[-1, :] = b[:, 0] = b[:, -1] = 0     

u0 = np.empty((nx, ny))
u0[0, :] = u0[-1, :] = u0[:, 0] = u0[:, -1] = 0     

u_plus_1 = u0 

max_it = 20000
alpha = 1

for it in range(max_it):
    # print('Iteration-{}'.format(it))
    u = u_plus_1 
    u_plus_1 = utils.jacobi(alpha, u, b, dx=dx, dy=dy)

# %%
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams['font.size'] = 22
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

fig1 = plt.figure(figsize=(8,6), dpi=100)

plt.contourf(X, Y, u_plus_1, cmap='jet', levels=32)
plt.xlabel('x')
plt.ylabel('y')
plt.colorbar(fraction=0.046)
