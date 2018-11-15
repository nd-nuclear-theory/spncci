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
sum = 0

for i in InFile:
	#if (i[3] == 0):
	rme_matrix[int(i[1])][int(i[2])] = rme_matrix[int(i[1])][int(i[2])] + i[-1]*i[-1]
	sum = sum + i[-1]*i[-1]

rme_matrix = rme_matrix/sum

print(np.sum(rme_matrix))

plt.imshow(rme_matrix, norm=LogNorm())
plt.colorbar()

plt.xlabel("mu")
plt.ylabel("lamda")
plt.title("RME untruncated magnitude")

ax = plt.gca()

# Major ticks
ax.set_xticks(np.arange(0, 21, 1));
ax.set_yticks(np.arange(0, 21, 1));

# Labels for major ticks
ax.set_xticklabels(np.arange(0, 21, 1));
ax.set_yticklabels(np.arange(0, 21, 1));

# Minor ticks
ax.set_xticks(np.arange(-.5, 21, 1), minor=True);
ax.set_yticks(np.arange(-.5, 21, 1), minor=True);

# Gridlines based on minor ticks
plt.grid(which='minor', linestyle='-', linewidth=2)


plt.savefig("heat_rmemag_lammu_log_sq.png")

plt.show(block=True)
