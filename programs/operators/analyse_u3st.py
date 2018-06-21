# Prototype
# Script to find the number of rmes u3st truncated and
# unique N(lambda, mu) S0 T0 rows after truncation
# Robert Power
# 06/20/2018

import numpy as np
import sys
import matplotlib.pyplot as plt

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

print("Number of rmes u3st truncated = " + str(len(InFile)))

print("Unique N(lambda, mu) S0 T0 rows = " + str(len(unique)))
