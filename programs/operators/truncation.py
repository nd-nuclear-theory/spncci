# Prototype
# Script comaring the rmes of truncated and original dat files.
# Robert Power
# 06/20/2018

import numpy as np
import sys

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

for i in unique:
	print(i)

