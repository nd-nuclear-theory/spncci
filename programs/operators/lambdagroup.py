# Prototype
# Code comparing the rmes of the same lambda value
# Robert Power
# 06/20/2018

import numpy as np
import matplotlib.pyplot as plt

InFile = np.loadtxt('u3st_N0_04_double_1.0e-04.dat')

lamgroup = [[] for i in range(21)]


for i in InFile:
	index = int(i[1])
	lamgroup[index].append(i[-1])

for i in lamgroup:
	print(len(i))

#output = np.column_stack(lamgroup)

#np.savetxt('templambda.dat', output, delimiter='\t', fmt=' '.join(['%.2e']*21), 
#	header='lambda = 0 \t 1 \t 2 \t 3 \t 4 \t 5 \t 6 \t 7 \t 8 \t' +
#			'9 \t 10 \t 11 \t 12 \t 13 \t 14 \t 15 \t 16 \t 17 \t ' +
#			'18 \t 19 \t 20 \t')
for i in range(len(lamgroup)):
	print("lambda = " + str(i))

	print("Max difference = " + str(max(lamgroup[i])))

	print("Mean difference = " + str(sum(lamgroup[i])/len(lamgroup[i])))

	print("------------------------------------")

for i in InFile:
	if (i[1] - i[2] != i[0]):
		print("True")


# the histogram of the data
for i in range(len(lamgroup)):
	plt.scatter(np.arange(len(lamgroup[i])), np.log10(np.array(lamgroup[i])))
	plt.savefig("scatter" + str(i) + ".png")
	plt.show()
	plt.close()