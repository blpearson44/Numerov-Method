#!/usr/bin/env python3
"""Iteratively find the wave funciton and energy for a particle."""

# Imports
import numpy as np
import matplotlib.pyplot as plt

# Constants
V0 = 83  # MeV
A = 2  # fm
mesh = np.linspace(
    -4, 4, 500
)  # grid is from -4fm to 4fm with 100 equally spaced pieces
C = 0.04829  # MeV^-1 fm^-2
# Only used prior to implementation of Energy finding function
E = -74.876956  # MeV


def numerov_method_negative(energy, grid):
    """Implement Numerov Method at index n to find n+1. From R -> L."""
    global A
    dx = grid[1] - grid[0]

    k = np.zeros(grid.size)

    # this will be the wave function when it is run
    val = np.zeros(grid.size)
    val[0] = 0
    val[1] = dx

    h_12 = dx ** 2 / 12

    # to simplify this, I am absorbing h**2/12 into this term
    k = []
    for i in range(grid.size):
        if abs(grid[i]) > A:
            k.append(h_12 * C * energy)
        else:
            k.append(h_12 * C * (energy + V0))

    print(k[0])

    # iterate through the array
    i = 2
    while grid[i] < - A:
        val[i] = 2 * (1 - 5 * k[i-1]) * val[i-1]\
            - (1 + k[i-2]) * val[i-2]
        val[i] = val[i]/(1 + k[i+1])
        i += 1

    val[val==0] = np.nan
    return val


def numerov_method_positive(energy, grid):
    """Implement Numerov Method at index n to find n+1. From R -> L."""
    dx = grid[1] - grid[0]

    k = np.zeros(grid.size)

    # this will be the wave function when it is run
    val = np.zeros(grid.size)
    val[-1] = 0
    val[-2] = dx

    h_12 = dx ** 2 / 12

    # to simplify this, I am absorbing h**2/12 into this term
    k = []
    for i in range(grid.size):
        if abs(grid[i]) > A:
            k.append(h_12 * C * energy)
        else:
            k.append(h_12 * C * (energy + V0))

    print(k[0])

    # iterate through the array
    i = grid.size - 3
    while grid[i] > -A:
        val[i] = 2 * (1 - 5 * k[i+1]) * val[i+1]\
            - (1 + k[i+2]) * val[i+2]
        val[i] = val[i]/(1 + k[i-1])
        i -= 1

    val[val==0] = np.nan
    return val


# data processing

# graphs
plt.style.use('ggplot')

fig = plt.figure(
    figsize=(9, 6),
    dpi=300
)

test_energies = [-80, -70, -74, -75]
axes = [plt.subplot(1, 4, i+1) for i in range(4)]
for i in range(4):
    ax = axes[i]
    psi = numerov_method_positive(test_energies[i], mesh)
    psi2 = numerov_method_negative(test_energies[i], mesh)
    ax.plot(
        mesh, psi,
        ls='none',
        ms=4,
        marker='.',
        mew=1.0,
        zorder=1,
        label=r'Backward integration'
    )
    ax.plot(
        mesh, psi2,
        ls='none',
        ms=4,
        marker='.',
        mew=1.0,
        zorder=1,
        label=r'Backward integration'
    )

ax.set_xlabel('x', labelpad=5)
ax.set_ylabel('y', labelpad=5)
plt.show()
plt.savefig('mismatch-energies.png')
