#!/usr/bin/ipython3

import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

plt.rcParams.update({'font.size': 22})

files = os.listdir("../data")
num_dumps = len(files)

cur_dump_file = h5py.File("../data/dump_0.h5", "r")
positions  = cur_dump_file["positions"][:]
velocities = cur_dump_file["velocities"][:]
grid_pos   = cur_dump_file["grid_pos"][:]

delta_x = grid_pos[2] - grid_pos[1]

ims = []

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(14,18), sharex=True)
fig.set_tight_layout(True)

line1, = ax1.plot([0.5], [0], "bo", label="Electrons  ")
ax1.set_xlim(0,1)
ax1.set_ylim(-2e5,2e5)

label_line, = ax1.plot([], [], " ", label="Hello")
label_leg = ax1.legend(loc=1)

ax1.set_xlabel("Position")
ax1.set_ylabel("Velocity")

line4, = ax2.plot([0],[0], label="rho")
line5, = ax2.plot([0],[0], label="rho hd")
ax2.legend(loc=1)
ax2.set_ylim(-2e-9,2e-9)

ax2.set_ylabel("Charge density")

line2, = ax3.plot([0.5], [0], label="phi")
line3, = ax3.plot([0.5], [0], label="$E\Delta x$")
ax3.legend(loc=1)
ax3.set_xlim(0,1)
ax3.set_ylim(-1,1)

ax3.set_ylabel("Voltage")

for g in grid_pos:
    ax1.axvline(g, color="black", linewidth=0.5, dashes=[10, 10])
    ax2.axvline(g, color="black", linewidth=0.5, dashes=[10, 10])
    ax3.axvline(g, color="black", linewidth=0.5, dashes=[10, 10])

def animate(i):
    cur_dump_file = h5py.File("../data/dump_{0}.h5".format(i), "r")
    positions  = cur_dump_file["positions"][:]
    velocities = cur_dump_file["velocities"][:]

    grid_pos = cur_dump_file["grid_pos"][:]
    rho = cur_dump_file["rho"][:]
    phi =  cur_dump_file["phi"][:]
    elec = delta_x * cur_dump_file["elec"][:]

    grid_pos_hd = cur_dump_file["grid_pos_hd"][:]
    rho_hd = cur_dump_file["rho_hd"][:]

    line1.set_xdata(positions)
    line1.set_ydata(velocities)

    line2.set_xdata(grid_pos)
    line2.set_ydata(phi)

    line3.set_xdata(grid_pos)
    line3.set_ydata(elec)

    line4.set_xdata(grid_pos)
    line4.set_ydata(rho)

    line5.set_xdata(grid_pos_hd)
    line5.set_ydata(rho_hd)

    label_line.set_label('Dump {0}'.format(i))

    global label_leg
    label_leg.remove()
    label_leg = ax1.legend(loc=1)

    return line1, line2, line3, line4, line5, label_line, label_leg

ani = animation.FuncAnimation(fig, animate, frames=np.arange(0, num_dumps), interval=200)

ani.save("dist_func.gif", dpi=80, writer="imagemagick")
