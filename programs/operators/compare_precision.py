# Script for comparing double and quad precision
# Prototype
# Robert Power
# 06/05/18

import numpy as np

# Read .dat files
DataDouble = np.loadtxt('relative_double.dat')
DataQuad = np.loadtxt('relative_quad.dat')

double_rme = []
quad_rme = []

# Double precision rmes
for i in DataDouble:
	double_rme.append(i[-1])

# Quad precision rmes
for i in DataQuad:
	quad_rme.append(i[-1])

rme_difference = []

# Absolute Difference
for i in range(len(DataDouble)):
	rme_difference.append(np.fabs(quad_rme[i] - double_rme[i]))

output = DataDouble.tolist()

for i in range(len(DataDouble)):
	output[i].append(quad_rme[i])
	output[i].append(rme_difference[i])

# Write .dat file
np.savetxt('compare_precision.dat', output, delimiter='\t',
			fmt=' '.join(['%d']*12) + ", %.16e, %.16e, %.16e",
			header='Relative U3ST - Comparing Doubel and Quad precision\nN\' S\' T\' N S T lamda mu S0 T0 kappa0 L0 double quad difference\n')
