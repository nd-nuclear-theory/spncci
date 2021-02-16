##################################################################################
# Utilities
#
#   Anna E. McCoy
#   University of Notre Dame
#
#   SPDX-License-Identifier: MIT
##################################################################################
import math
import so3

##################################################################################################
def parity(n):
    """
    parity(n) determines the value (-)^n
    input
        n: integer
    output
        p: sign
    """
    if n%2==0:
        p=1
    else:
        p=-1
    return p

########################################################################################################################################################################################
def choose(n,r):
  	"""Computes n! / (r! (n-r)!) exactly. Returns a python long int."""
  	c=0
  	if (n>=0) and (0 <= r <= n):
	  	c = 1L
	  	denom = 1
	  	for (num,denom) in zip(xrange(n,n-r,-1), xrange(1,r+1,1)):
			c = (c * num) // denom
	return c

if __name__ == '__main__':
 	print choose(40,12)
############################################################################################

def d(J,Mp,M):
	"""
	Little d function defined according to mathematica documention
	Argument is pi/2
	Can use the scipy.misc.comb which can be give the option exact but may be slower. 

	uses convention for little d function as defined in Edmonds
	"""
	dd=0.0
	if so3.multJ(J,M)==1 and so3.multJ(J,Mp)==1:
		## defining the limits of k
		kmax=max(J+M,J-M,J-Mp)
		#kmin=max(Mp-M,0)
		##suming over k
		for k in range(int(kmax)+1):
			dd=dd+parity(k)*choose(int(J+M),int(J-Mp-k))*choose(int(J-M),k)
		dd=(parity(J-Mp)
			*math.sqrt(
				math.factorial(int(J+Mp))
				*math.factorial(int(J-Mp))
				/(
					2.**(2*J)
					*math.factorial(int(J+M))
					*math.factorial(int(J-M))
					)
				)
			)*dd
	#print "dd",dd
	return dd

##################################################################################################################

def delta(a,b):
	if a==b:
		d=1
	else:
		d=0
	return d

##################################################################################################################


if __name__ == '__main__':

	# print d(.5,.5,.5), 1/numpy.sqrt(2)
	# print d(.5,.5,-.5), -1/numpy.sqrt(2)
	# print d(.5,-.5,.5),1/numpy.sqrt(2)
	print [1,1],  d(1,1,1), (1./2)
	print [1,0],  d(1,1,0),   1./math.sqrt(2)
	print [1,-1], d(1,1,-1),  1./2
	print [0,1],  d(1,0,1),  -1./math.sqrt(2)
	print [0,0],  d(1,0,0),   0
	print [0,-1], d(1,0,-1),  1./math.sqrt(2)
	print [-1,1], d(1,-1,1),  1./2
	print [-1,0], d(1,-1,0), -1./math.sqrt(2)
	print [-1,-1],d(1,-1,-1), 1./2


if __name__ == '__main__':
    print "parity"
    print "2", parity(2)
    print "5", parity(5)



