# Prototype
# Script comaring the rmes of truncated and original dat files.
# Robert Power
# 06/20/2018

import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

InFile = np.loadtxt('../../data/relative_interactions/jisp16_Nmax20_hw20.0_rel.dat', skiprows=8)
diff = []

np.savetxt('tempdiff.dat', diff, delimiter='\t')

unique = []
new = True

for f in InFile:
	new = True
	labels = (int(f[1]), int(f[2]), int(f[3]), int(f[4]), int(f[5]))
	for u in unique:
		if labels == u:
			new = False
			break
	if new is True:
		unique.append(labels)

for i in range(len(unique)):
	print(unique[i])

print(len(unique))

label_map = {}

for i in range(len(unique)):
	label_map[unique[i]] = i

diff_matrix = np.zeros((len(unique), len(unique)))

#for j in range(len(InFile)):
#	i = InFile[j]
#	if (i[3] == 1 and i[5] == 0):
#		diff_matrix[label_map[(i[1], i[2], i[3], i[4], i[5])]][label_map[(i[6], i[7], i[8], i[9], i[10])]] = np.fabs(i[-1])
#	else:
#		diff_matrix[label_map[(i[1], i[2], i[3], i[4], i[5])]][label_map[(i[6], i[7], i[8], i[9], i[10])]] = 0

#print(diff_matrix[0][5])

#plt.imshow(diff_matrix, norm=LogNorm(vmin=0.01, vmax=1))
#plt.colorbar()

#plt.title("RME untruncated magnitude")

#plt.savefig("heat_rmemag_s1t0.png")

#plt.show(block=True)