# Script for comparing double and quad precision for a range of zero thresholds
# Prototype
# Robert Power
# 06/09/18

import numpy as np

# Read .dat files
double_data_files = [""]*28
quad_data_files = [""]*28

double_data_files[0] = "double_1.0e-10.dat"

for i in range(9):
	double_data_files[i+1] = "double_" + str(i+1) + ".0e-03.dat"
	#quad_data_files[i] = "quad_" + str(i+1) + ".0e-01.dat"

for i in range(9):
	double_data_files[i+10] = "double_" + str(i+1) + ".0e-02.dat"
	#quad_data_files[i] = "quad_" + str(i+1) + ".0e-01.dat"

for i in range(9):
	double_data_files[i+19] = "double_" + str(i+1) + ".0e-01.dat"
	#quad_data_files[i] = "quad_" + str(i+1) + ".0e-01.dat"

for i in range(28):
	DataDouble = np.loadtxt(double_data_files[i], dtype=str)
	print(str(double_data_files[i]) + " - " + str(DataDouble[-1][-1]))
	#print(quad_data_files[i])

double_group = []
double_rme = []
DataDouble = []

#for i in double_data_files:
#	DataDouble = np.loadtxt(i, dtype=str)
#	print(DataDouble[-1][-1])
#	for j in DataDouble:
#		double_rme.append(j[-1])
#	double_group.append(double_rme)
#
#quad_group = []
#quad_rme = []
#DataQuad = []

#for i in quad_data_files:
#	DataQuad = np.loadtxt(i)
#	for j in DataQuad:
#		quad_rme.append(j[-1])
#	quad_group.append(quad_rme)
#
#threshold_compare = []
#compare = []

#for i in range(11):
#	for j in range(len(double_group[i])):
#		threshold_compare.append(np.fabs(double_group[i][j] - quad_group[i][j]))
#	compare.append(threshold_compare)
#	threshold_compare = []

#for i in range(9):
#	print(compare[i][4])
