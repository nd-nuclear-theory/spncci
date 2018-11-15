# Prototype
# Code comparing the rmes of the same mu value
# Robert Power
# 06/20/2018

import numpy as np

InFile = np.loadtxt('u3st_N0_20_double_1.0e-04.dat')

mugroup = [[] for i in range(21)]


for i in InFile:
	index = int(i[2])
	mugroup[index].append(i[-1])

max_n = max([len(x) for x in mugroup])

for i in mugroup:
	length = len(i)
	for j in range(max_n - length):
		i.append(0)

for i in mugroup:
	print(len(i))

output = np.column_stack(mugroup)

#np.savetxt('mu.dat', output, delimiter='\t', fmt=' '.join(['%.2e']*21), 
#	header='mu = 0 \t 1 \t 2 \t 3 \t 4 \t 5 \t 6 \t 7 \t 8 \t' +
#			'9 \t 10 \t 11 \t 12 \t 13 \t 14 \t 15 \t 16 \t 17 \t ' +
#			'18 \t 19 \t 20 \t')
for i in range(len(mugroup)):
	print("mu = " + str(i))

	print("Max difference = " + str(max(mugroup[i])))

	print("Mean difference = " + str(sum(mugroup[i])/len(mugroup[i])))

	print("------------------------------------")
