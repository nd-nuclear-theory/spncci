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

#ordered = sorted(InFile,key=lambda InFile:np.fabs(InFile[-1]), reverse=False)

#np.savetxt('u3_ordered.dat', ordered, delimiter='\t',
#			fmt=' '.join(['%d']*8) + ", %.16e")
#			
sum = 0.0
sumbyN = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

for i in InFile:
	sum = sum + np.fabs(i[-1])
	factor = (20.0 + i[0])/2.0
	sumbyN[int(factor)] = sumbyN[int(factor)] + np.fabs(i[-1])

print(sumbyN)
print(np.sum(sumbyN))
print(sum)

weight = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
for s in range(len(sumbyN)):
	weight[s] = sumbyN[s]/sum

print(weight)

print(np.sum(weight))

for s in range(len(weight)):
	print("N = " + str(-20.0 + s*2) + " - " + str(weight[s]))

weight = np.array(weight)
vals = np.arange(0, 21, 1)*2 - 20
plt.pie(weight, labels=vals)
plt.axis('equal')
plt.title('N0 values')
plt.savefig('pie_N0.png')
plt.show()
plt.close()

###########################################################################

sum = 0.0
sumbylam = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

for i in InFile:
	sum = sum + np.fabs(i[-1])
	sumbylam[int(i[1])] = sumbylam[int(i[1])] + np.fabs(i[-1])

print(sumbylam)
print(np.sum(sumbylam))
print(sum)

weight = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
for s in range(len(sumbylam)):
	weight[s] = sumbylam[s]/sum

print(weight)

print(np.sum(weight))

for s in range(len(weight)):
	print("Lambda = " + str(s) + " - " + str(weight[s]))

weight = np.array(weight)
vals = np.arange(0, 21, 1)
plt.pie(weight, labels=vals)
plt.axis('equal')
plt.title('Lambda values')
plt.savefig('pie_lam.png')
plt.show()
plt.close()


###########################################################################

sum = 0.0
sumbymu = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

for i in InFile:
	sum = sum + np.fabs(i[-1])
	sumbymu[int(i[2])] = sumbymu[int(i[2])] + np.fabs(i[-1])

print(sumbymu)
print(np.sum(sumbymu))
print(sum)

weight = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
for s in range(len(sumbymu)):
	weight[s] = sumbymu[s]/sum

print(weight)

print(np.sum(weight))

for s in range(len(weight)):
	print("Mu = " + str(s) + " - " + str(weight[s]))

weight = np.array(weight)
vals = np.arange(0, 21, 1)
plt.pie(weight, labels=vals)
plt.axis('equal')
plt.title('Mu values')
plt.savefig('pie_mu.png')
plt.show()
plt.close()

###########################################################################

sum = 0.0
sumbyS = [0, 0, 0]

for i in InFile:
	sum = sum + np.fabs(i[-1])
	sumbyS[int(i[3])] = sumbyS[int(i[3])] + np.fabs(i[-1])

print(sumbyS)
print(np.sum(sumbyS))
print(sum)

weight = [0, 0, 0]
for s in range(len(sumbyS)):
	weight[s] = sumbyS[s]/sum

print(weight)

print(np.sum(weight))

for s in range(len(weight)):
	print("S = " + str(s) + " - " + str(weight[s]))

weight = np.array(weight)
vals = np.arange(0, 3, 1)
plt.pie(weight, labels=vals)
plt.axis('equal')
plt.title('Mu values')
plt.savefig('pie_S.png')
plt.show()
plt.close()