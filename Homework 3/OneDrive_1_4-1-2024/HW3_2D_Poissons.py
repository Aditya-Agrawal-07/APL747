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

np.random.seed(2)

def gauss_seidel(alpha, u, b, dx, dy):
    nx, ny = [*u.shape]
    u_plus_1 = u.copy()
    # apply evolution operator
    for j in range(1, ny-1):
        for i in range(1, nx-1):
            u_plus_1[i, j] = (0.25 * (alpha * (u_plus_1[i-1, j] + u[i+1, j] + u_plus_1[i, j-1]
                                               + u[i, j+1]) - b[i, j]*(dx*dy)))    
    return u_plus_1

def jacobi(alpha, u, b, dx, dy):
    u_plus_1 = np.zeros(u.shape)
    # apply evolution operator
    u_plus_1[1:-1, 1:-1] = (0.25 * (alpha * (u[:-2, 1:-1] + u[2:, 1:-1] + u[1:-1, :-2]
                                             + u[1:-1, 2:]) - b[1:-1, 1:-1]*(dx*dy))) 
    # copy boundary conditions
    u_plus_1[0,:]  = u[0,:]
    u_plus_1[-1,:] = u[-1,:]
    u_plus_1[:,0]  = u[:,0]
    u_plus_1[:,-1] = u[:,-1]
    return u_plus_1

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
b = np.ones((nx, ny))
b[0, :] = b[-1, :] = b[:, 0] = b[:, -1] = 0     # Boundary conditions

u0 = np.zeros((nx, ny))
u_plus_1 = u0 

mcs = 2000
max_it = 10000
sol = []
for samples in range(mcs):
    alpha = np.random.uniform(-1,1)
    print('Sample-{}, Alpha-{:0.4f}'.format(samples, alpha))
    
    for it in range(max_it):
        u = u_plus_1.copy()
        u_plus_1 = jacobi(alpha, u, b, dx=dx, dy=dy)
    sol.append(u_plus_1)

sol = np.stack(sol)
mean = np.nanmean(sol, axis=0)
std = np.nanstd(sol, axis=0)

# %%
fig1 = plt.figure(figsize=(10,4), dpi=100)
plt.subplots_adjust(wspace=0.5)

plt.subplot(1,2,1)
plt.contourf(X, Y, mean, cmap='jet', levels=32)
plt.colorbar(fraction=0.046)

plt.subplot(1,2,2)
plt.contourf(X, Y, std, cmap='jet', levels=32)
plt.colorbar(fraction=0.046)
