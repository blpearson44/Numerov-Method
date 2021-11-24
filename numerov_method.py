#!/usr/bin/env python3
"""Iteratively find the wave funciton and energy for a particle."""

# Imports
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Constants
V0 = 83  # MeV
A = 2  # fm
C = 0.04829  # MeV^-1 fm^-2


def create_karr(energy: float, grid: np.linspace) -> list:
    """Create array for K**2 at each point in grid."""
    dx = grid[0] - grid[1]
    h_12 = dx ** 2 / 12
    arr = []
    for i in range(grid.size):
        if abs(grid[i]) > A:
            arr.append(h_12 * C * energy)
        else:
            arr.append(h_12 * C * (energy + V0))
    return arr


def numerov_method_forward(energy: float, grid: np.linspace) -> np.linspace:
    """Implement Numerov Method at index n to find n+1. From R -> L."""
    global A
    k = create_karr(energy, grid)
    # this will be the wave function when it is run
    val = np.zeros(grid.size)
    val[0] = 0
    val[1] = 0.1

    # iterate through the array
    i = 2
    while grid[i] < A:
        val[i] = 2 * (1 - 5 * k[i-1]) * val[i-1]\
            - (1 + k[i-2]) * val[i-2]
        val[i] = val[i]/(1 + k[i+1])
        i += 1

    val[val == 0] = np.nan
    return val


def numerov_method_backward(energy: float, grid: np.linspace) -> np.linspace:
    """Implement Numerov Method at index n to find n+1. From L -> R."""
    global A
    k = create_karr(energy, grid)
    # this will be the wave function when it is run
    val = np.zeros(grid.size)
    val[-1] = 0
    val[-2] = 0.1

    # iterate through the array
    i = grid.size - 3
    while grid[i] > A:
        val[i] = 2 * (1 - 5 * k[i+1]) * val[i+1]\
            - (1 + k[i+2]) * val[i+2]
        val[i] = val[i]/(1 + k[i-1])
        i -= 1

    val[val == 0] = np.nan
    return val


# data processing
# graphs
plt.style.use('ggplot')


test_energies = [-74.876956, -51.282727, -16.0313114]
n = [10, 100, 500]
for energy in test_energies:
    fig = plt.figure(
        figsize=(9, 6),
        dpi=300
    )
    plot = [plt.subplot(1, 3, i+1) for i in range(3)]
    for subplot in range(len(plot)):
        grid = np.linspace(-4, 4, n[subplot])
        plot[subplot].plot(
            grid, numerov_method_forward(energy, grid),
            ls='none',
            ms=2, marker='o', mfc='white', mew=1.0,
            zorder=1,
            label=r'Forward integration'
        )
        plot[subplot].plot(
            grid, numerov_method_backward(energy, grid),
            ls='none',
            ms=2, marker='o', mfc='white', mew=1.0,
            zorder=1,
            label=r'Backward integration'
        )
        plot[subplot].set_xlabel(r'$x$', labelpad=5)
        plot[subplot].set_ylabel(r'$\psi$', labelpad=5)
    fig.tight_layout(h_pad=2)
    fig.suptitle(f"E = {energy}MeV")
    plt.subplots_adjust(top=0.85)
    forward_integration = mpatches.Patch(color='orange',
                                         label='Forward integration')
    backward_integration = mpatches.Patch(color='blue',
                                          label='Backward integration')
    fig.legend(handles=[forward_integration, backward_integration],
               loc='upper left', fontsize='small')
    plt.savefig(f"{energy}.png")
