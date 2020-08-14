#!/usr/bin/ipython3

import h5py
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 22})

f = h5py.File("../data/dump_50.h5", "r")
grid_pos   = f["grid_pos"][:]
rho        = f["rho"][:]
grid_pos_hd   = f["grid_pos_hd"][:]
rho_hd        = f["rho_hd"][:]

dx = grid_pos[2] - grid_pos[1]
dx_hd = grid_pos_hd[2] - grid_pos_hd[1]

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14,9))
fig.set_tight_layout(True)

ax1.plot(grid_pos, rho)
ax1.plot(grid_pos_hd, rho_hd)

ax2.axhline(0, color="black", linewidth=0.5)
ax2.plot(np.fft.fftfreq(len(rho), dx), np.fft.fft(rho), "x", label="rho", zorder=10)
ax2.plot(np.fft.fftfreq(len(rho_hd), dx_hd), np.fft.fft(rho_hd), "o", label="rho hd")
ax2.legend()

for i in range(-3,3):
    j = 16 * (2*i + 1)
    ax2.axvline(j, color="black", linewidth=1)

for i in range(-2,3):
    k = 16 * (2*i)
    ax2.axvline(k+2, color="black", linewidth=0.5, dashes=[10,10])
    ax2.axvline(k-2, color="black", linewidth=0.5, dashes=[10,10])

fig.savefig("rho_ft.pdf")
