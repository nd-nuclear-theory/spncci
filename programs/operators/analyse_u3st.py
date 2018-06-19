import numpy as np
import sys

#InFile = np.loadtxt('u3st_N0_10_double_1.0e-04.dat')
InFile = np.loadtxt(sys.argv[1])

unique = []
new = True

for f in InFile:
	new = True
	prefix = (int(f[0]), int(f[1]), int(f[2]), int(f[3]), int(f[4]))
	for u in unique:
		if prefix == u:
			new = False
			break
	if new is True:
		unique.append(prefix)
	print(prefix)

print(len(unique))
