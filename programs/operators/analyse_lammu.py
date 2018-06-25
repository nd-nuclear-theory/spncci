# Prototype
# Script comaring the rmes of truncated and original dat files.
# Robert Power
# 06/20/2018

import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

InFile = np.loadtxt('u3st_N0_20_double_1.0e-04.dat')

rme_matrix = np.zeros((21, 21))

for i in InFile:
	#coord = rme_matrix[int(i[1])][int(i[2])]
	#if (np.fabs(i[-1]) > coord):
	rme_matrix[int(i[1])][int(i[2])] = rme_matrix[int(i[1])][int(i[2])] + i[-1]

print(rme_matrix[0][5])

plt.imshow(rme_matrix, norm=LogNorm(vmin=0.01, vmax=1))
plt.colorbar()

plt.xlabel("mu")
plt.ylabel("lamda")
plt.title("RME untruncated magnitude")

# Major ticks
plt.set_xticks(np.arange(0, 21, 1));
plt.set_yticks(np.arange(0, 21, 1));

# Labels for major ticks
plt.set_xticklabels(np.arange(0, 21, 1));
plt.set_yticklabels(np.arange(0, 21, 1));

# Minor ticks
plt.set_xticks(np.arange(-.5, 21, 1), minor=True);
plt.set_yticks(np.arange(-.5, 21, 1), minor=True);

# Gridlines based on minor ticks
plt.grid(which='minor', linestyle='-', linewidth=2)


plt.savefig("heat_rmemag_lammu_sum_log.png")

plt.show(block=True)
