# Prototype
# Script comaring the rmes of truncated and original dat files.
# Robert Power
# 06/20/2018

import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

InFile = np.loadtxt('../../data/relative_interactions/jisp16_Nmax20_hw20.0_rel.dat', skiprows=8)
OutFile = np.loadtxt(sys.argv[1], skiprows=8)
diff = []


def modeOrder(l):
	order = 0
	orderList = []
	for i in l:
		temp = i
		if (i != 0):
			while (temp < 1):
				temp = temp*10
				order = order - 1
			orderList.append(order)
			order = 0
	return(max(set(orderList), key=orderList.count))


for i in range(len(InFile)):
	diff.append(np.fabs(InFile[i][-1] - OutFile[i][-1]))

print("Max difference = " + str(max(diff)))

print("Mean difference = " + str(sum(diff)/len(diff)))

print("Median difference = " + str(np.median(diff)))

print("Order of magnitude mode of difference = " + str(modeOrder(diff)))

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

print(label_map[(12,0,0,0,1)])

diff_matrix = np.zeros((len(unique), len(unique)))

for j in range(len(OutFile)):
	i = InFile[j]
	o = OutFile[j]
	diff_matrix[label_map[(i[1], i[2], i[3], i[4], i[5])]][label_map[(i[6], i[7], i[8], i[9], i[10])]] = np.fabs(i[-1] - o[-1])

print(diff_matrix[0][5])

plt.imshow(diff_matrix, norm=LogNorm())
plt.colorbar()

plt.title("RME " + sys.argv[1][32:-4] + " truncation difference magnitude")

plt.savefig("heat_" + sys.argv[1][32:-4] + ".png")

plt.show(block=True)