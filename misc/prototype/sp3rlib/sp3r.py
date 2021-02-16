##################################################################################
#	Anna E. McCoy
#	University of Notre Dame
#
#	SPDX-License-Identifier: MIT
##################################################################################

import numpy
import numpy.linalg
from u3states import *
from spstates import *
import so3
import u3
import u3boson
import vcs
import time
import sys
numpy.set_printoptions(suppress=True,threshold=10000, linewidth=120)
#numpy.set_printoptions(threshold=10000, linewidth=120,formatter={"all": lambda x:format(x,"10.4f")})

A=Sp3RSU3Operator(U3State([2,0,0]),vcs.A,False)
B=Sp3RSU3Operator(U3State([0,0,-2]),vcs.B,False)
C=Sp3RSU3Operator(U3State([1,0,-1]),vcs.C,True)
Cn=Sp3RSU3Operator(U3State([1,0,-1]),vcs.C_n,False)
Cs=Sp3RSU3Operator(U3State([1,0,-1]),vcs.C_s1,False)
Cs=Sp3RSU3Operator(U3State([1,0,-1]),u3boson.C_s1,False)

#Csw=Sp3RSU3Operator(U3State([1,0,-1]),u3boson.Cs_weylboson2,False)

def CoupledOperators_Sp3R(operator1,operator2,w0,r0,psi2,psi1,r0p):
	s=psi1.s
	w1=psi1.w
	n1=psi1.n
	r1=psi1.r
	w2=psi2.w
	n2=psi2.n
	r2=psi2.r
	matrix=0.0
	m1=0
	m2=0
	m3=0
   	if (w1.N+operator1.u3.N+operator2.u3.N)!=w2.N:
   		return matrix
    #coupling w to first tensor
    #Checking if operators are C operators and restricting w12 accordingly
	w12set=u3.couple(w1,operator2.u3)
	r1max=u3.mult(s,n1,w1)
	r2max=u3.mult(s,n2,w2)
	for lambdamu12 in w12set:
		mult1=u3.mult(lambdamu12,operator1.u3,w2)
		#mult1=min(mult1max,r2max)
		if mult1==0:
			continue
		mult2=w12set[lambdamu12]
		#mult2=min(mult2max, r1max)
		w12=U3State(w1.N+operator2.u3.N,lambdamu12)
		n12set=u3boson.nlist(s,w12)
		if n12set=={}:
			continue
		for n12 in n12set:
			r12max=n12set[n12]
			#mult2=min(mult2max,r12max)
			for r12 in range(1,r12max+1):				
				psi12=SpStateReduced(s,n12,r12,w12)
				for r1p in range(1,mult1+1):
					for r2p in range(1,mult2+1):
						matrix=matrix+(
							u3.U(w1,operator2.u3, w2, operator1.u3, w12,r2p,r1p,w0,r0,r0p)
							*operator1.SU3RME(psi2,psi12,r1p)
							*operator2.SU3RME(psi12,psi1,r2p)
							)	
						
	return matrix
# if __name__ == '__main__':
# 	s=U3State([8,4,4])
# 	n=U3State([6,4,2])
# 	w=U3State([15,9,6])
# 	wset=u3.couple(s,n)
# 	psired=SpStateReduced(s,n,1,w)
# 	for lm in wset:
# 		w=U3State(s.N+n.N,lm)
# 		psired=SpStateReduced(s,n,1,w)
# 		Cs=Sp3RSU3Operator(U3State([1,0,-1]),vcs.C_s2,False)
# 	# print Cs.SU3RME(psired,psired,1)
# 	# print Cs.SU3RME(psired,psired,2)
# 		print CoupledOperators_Sp3R(Cs,Cs,SU3State(0,0),1,psired,psired,1)
# 	print u3.casimir2_su3(s)/numpy.sqrt(SU3State(1,1).dim())
# 	print 

def CoupledOperators_Sp3R_SO3(operator1,k1,L1,operator2,k2,L2,k0,L0,psip,psi,c1=0,c2=0):
    matrix=0.0
    w=psi.w
    wp=psip.w
    w1=operator1.u3
    w2=operator2.u3
    lm0set=u3.couple(w1,w2)
    for lm0 in lm0set:
    	k0max=so3.multk(lm0,L0)
    	if k0max==0:
    		continue
    	r0pmax=u3.mult(w,lm0,wp)
    	r0max=lm0set[lm0]
    	for r0 in range(1,r0max+1):
    		for r0p in range(1,r0pmax+1):
    			for k0 in range(1,k0max+1):
    				matrix=matrix+(
						u3.W(operator2.u3,k2,L2,operator1.u3,k1,L1,lm0,k0,L0,r0,c1=c2,c2=c1)
						*u3.W(w,psi.k,psi.L,lm0,k0,L0,wp,psip.k,psip.L,r0p)
						*CoupledOperators_Sp3R(operator1,operator2,lm0,r0,psip.red(),psi.red(),r0p)
						)
				
    #matrix=numpy.sqrt(2*Lp+1)*matrix
    if abs(matrix)<(10**(-10)):
        matrix=0.0
    return matrix

# if __name__ == '__main__':
# 	# s=U3State([8,4,4])
# 	# r=1
# 	# n=U3State([4,0,0])
# 	# w=U3State([10,6,4])
# 	# k=1
# 	# L=0
# 	# np=U3State([4,0,0])
# 	# rp=1
# 	# wp=U3State([10,6,4])
# 	# kp=1
# 	# Lp=L
# 	# L0=0
# 	s=U3State([326,250,250])
# 	r=1
# 	n=U3State([8,0,0])
# 	w=U3State([330,254,250])
# 	k=1
# 	L=0
# 	np=U3State([8,0,0])
# 	rp=1
# 	wp=U3State([332,252,250])
# 	kp=1
# 	Lp=L
# 	L0=0
# 	psired=SpStateReduced(s,n,r,w)
# 	psipred=SpStateReduced(s,np,rp,wp)
# 	psi0=SpState(s,U3State([0,0,0]),1,s,1,L)
# 	psi1=SpState(s,n,r,w,k,L)
# 	psi2=SpState(s,np,rp,wp,kp,Lp)
# 	print CoupledOperators_Sp3R_SO3(C,1,1,C,1,1,1,0,psi1,psi2)
# 	print -L*(L+1)*numpy.sqrt(3)
# 	print "Cs"
# 	print CoupledOperators_Sp3R_SO3(Cs,1,1,Cs,1,1,1,0,psi0,psi0)
# 	print CoupledOperators_Sp3R_SO3(Cs,1,1,Cs,1,1,1,0,psi1,psi1)
# 	print CoupledOperators_Sp3R_SO3(Cs,1,1,Cs,1,1,1,0,psi2,psi2)
# 	print CoupledOperators_Sp3R_SO3(Cs,1,1,Cs,1,1,1,0,psi1,psi2)

# 	a1=CoupledOperators_Sp3R_SO3(B,1,0,A,1,0,1,L0,psi1,psi2,c1=1)
# 	a2=CoupledOperators_Sp3R_SO3(A,1,0,B,1,0,1,L0,psi1,psi2,c2=1)
# 	print "testing commutators [B,A]_0", a1-a2, psi1.w.N*2./3
# 	print CoupledOperators_Sp3R(B,A,SU3State(2,2),1,psipred,psired,1)
# 	print CoupledOperators_Sp3R(A,B,SU3State(2,2),1,psipred,psired,1)
# 	print "a^\dagger a"
# 	a=Sp3RSU3Operator(U3State([2,0,0]),vcs.A_u3boson,False)
# 	b=Sp3RSU3Operator(U3State([0,0,-2]),vcs.B_u3boson,False)
# 	aa=Sp3RSU3Operator(U3State([2,0,0]),u3boson.A,False)
# 	bb=Sp3RSU3Operator(U3State([0,0,-2]),u3boson.B,False)
# 	print SU3State(2,0).dim()
# 	print CoupledOperators_Sp3R_SO3(a,1,2,b,1,2,1,L0,psi1,psi2,c2=1)
# 	aabb2=CoupledOperators_Sp3R_SO3(aa,1,2,bb,1,2,1,L0,psi1,psi2,c2=1)
# 	aabb0=CoupledOperators_Sp3R_SO3(aa,1,0,bb,1,0,1,L0,psi1,psi2,c2=1)
# 	coef2=u3.W(aa.u3,1,2,bb.u3,1,2,SU3State(2,2),1,0,1,c2=1)
# 	coef0=u3.W(aa.u3,1,0,bb.u3,1,0,SU3State(2,2),1,0,1,c2=1)
# 	print coef2,coef0
# 	print coef2*CoupledOperators_Sp3R(aa,bb,SU3State(0,0),1,psi1.red(),psi2.red(),1),aabb2
# 	print coef0*CoupledOperators_Sp3R(aa,bb,SU3State(0,0),1,psi1.red(),psi2.red(),1),aabb0
# 	print numpy.sqrt(6)*CoupledOperators_Sp3R(aa,bb,SU3State(0,0),1,psi1.red(),psi2.red(),1)
# 	print psi1.n.N
# 	print "Cs.A Cs.B"
# 	s=U3State([36,10,10])
# 	n=U3State([8,0,0])
# 	w=U3State([40,14,10])
# 	r=1
# 	psiA=SpState(s,n,r,w,k,L)
# 	s=U3State([36,10,10])
# 	np=U3State([8,2,0])
# 	wp=U3State([42,14,10])
# 	rp=1
# 	psiAp=SpState(s,np,rp,wp,k,L)
# 	Cs=Sp3RSU3Operator(U3State([1,0,-1]),u3boson.C_s1,False)
# 	A=Sp3RSU3Operator(U3State([2,0,0]),vcs.A_u3boson,False)
# 	B=Sp3RSU3Operator(U3State([0,0,-2]),vcs.B_u3boson,False)
# 	print CoupledOperators_Sp3R_SO3(Cs,1,2,A,1,2,1,0,psiAp, psiA)
# 	print CoupledOperators_Sp3R_SO3(B,1,2,Cs,1,2,1,0,psiA, psiAp)
# 	print CoupledOperators_Sp3R_SO3(A,1,2,Cs,1,2,1,0,psiAp, psiA)
# 	print "testing A_2B_2"
# 	s=U3State([326,250,250])
# 	r=1
# 	n=U3State([2,0,0])
# 	w=U3State([328,250,250])
# 	k=1
# 	L=0
# 	np=U3State([2,0,0])
# 	L0=0
# 	psi1=SpState(s,n,r,w,k,L)
# 	print CoupledOperators_Sp3R_SO3(A,1,2,B,1,2,1,L0,psi1,psi1,c2=1)
# bb=CoupledOperators_Sp3R_SO3(Cs,1,2,Cn,1,2,1,0,psi2,psi1)
# aa=CoupledOperators_Sp3R_SO3(Cn,1,2,Cs,1,2,1,0,psi1,psi2,c1=1)
# print aa, bb
# # print A.SU3RME(psi2,psi1,1)*W_SU3(w.SU3(),k,L,(2,0),1,2,wp.SU3(),kp,Lp,1), B.SU3RME(psi1,psi2,1)*W_SU3(wp.SU3(),kp,Lp,B.SU3(),1,2,w.SU3(),k,L,1)
# print CoupledOperators_Sp3R(C,C,0,0,1,psipred,psired,1)
# print CoupledOperators_Sp3R(C,C,2,2,1,psipred,psired,1)

#print C_SU3(w.lbda(),w.mu())





# if __name__ == '__main__':

# 	s=U3State([6,4,4])
# 	n=U3State([4,0,0])
# 	r=1
# 	w=U3State([8,6,4])
# 	wp=U3State([8,6,4])
	
	# #print
	# print cas, cscs#, (cas/cscs)**2

###############################################################
def Sp_bandhead_check(s):
	"checks that s is a valid unitary U(3) irrep label and thus a valid Sp bandhead"
	val=False
	if s.check() and (s.three()>0):
		t=s.three()
		b1=s.one()-t
		b2=s.two()-t
		tb1=0
		tb2=0
		if b2>1:
			tb1=2
			tb2=2
		else:
			if b2==1:
				tb1=2
				if b1>1:
					tb2=1
		if tb1+tb2<=(2*t):
			val=True
	return val
###############################################################
def Sp_irrep_red(s,Nmax):
	"""
	constructs a dictioanry of symplectic reduced states sorted according to excitation Nn
	currently just for s.three()>3 but will generalize soon 
	returns dictionary of the form {N:[state1, state2...]}
	"""
	StateDict={}
	#checking that s is valid unitary U(3) label to serve as bandhead
	if Sp_bandhead_check(s):
		PolyDict=u3boson.poly(Nmax)
		for N in range(0,Nmax+1,2):
			polylist=PolyDict[N]
			states=[]
			for n in polylist:
				CoupledDict=u3.couple(s,n)
				for lambdamu in CoupledDict:
					w=U3State(s.N+N,lambdamu)
					rmax=CoupledDict[lambdamu]
					for r in range(1,rmax+1):
						states.append(SpStateReduced(s,n,r,w))
			StateDict[N]=states
	else:
		print "invalid bandhead"
	return StateDict	





# if __name__ == '__main__':
# 	s=U3State(836,SU3State(78,10))
# 	w=U3State(844,SU3State(78,8))
# 	print u3boson.nlist(s,w)






# if __name__ == '__main__':
# 	## Calculating the dimension of the largest K matrix at each Nex
# 	s=U3State(49,SU3State(5,4))
# 	Nmax=20
# 	weight_states={}
# 	polydict=u3boson.poly(Nmax)
# 	for N in range(0,Nmax+1,2):
# 		weight_states[N]={}
# 		polyN=polydict[N]
# 		for n in polyN:
# 			lmset=u3.couple(s,n)
# 			for lm in lmset:
# 				rmax=lmset[lm]
# 				w=U3State(s.N+n.N,lm)
# 				if not w.check():
# 					continue
# 				if w in weight_states[N]:
# 					weight_states[N][w]=weight_states[N][w]+rmax
# 				else:
# 					weight_states[N][w]=rmax
# 	total=0
# 	for N in range(0,Nmax+1,2):
# 		totsize=0
# 		for w in weight_states[N]:
# 			imax=weight_states[N][w]
# 			#print "{:10} {:2}".format(w,imax)
# 			size=0
# 			K,nitr=vcs.kmatrix(s,w)
# 			for i in range(0,imax):
# 				 for j in range(0,imax):
# 					size=size+sys.getsizeof(K[i,j])
# 			totsize=totsize+size
# 		total=total+totsize
# 		print "{:2} {:5}".format(N,totsize)
# 	print 
# 	print total
	# for N in range(0,Nmax+1,2):
	# 	rmaxmax=0
	# 	for w in weight_states[N]:
	# 		if weight_states[N][w]>rmaxmax:
	# 			rmaxmax=weight_states[N][w]
	# 			wmax=w
	# 	print "{:2}: {} {:3}".format(N,wmax,rmaxmax)



# if __name__ == '__main__':
# 	s=U3State([6,6,4])
# 	Spstates=Sp_irrep_red(s,6)
# 	for N in Spstates:
# 		for state in Spstates[N]:
# 			print state

# if __name__ == '__main__':
# 	s=U3State(836,SU3State(78,10))
# 	#print u3.couple(s,SU3State(4,2))
# 	n=U3State(8,SU3State(4,2))
# 	w=U3State(844,SU3State(83,7))
# 	np=U3State(10,SU3State(2,4))
# 	wp=U3State(846,SU3State(84,6))
# 	rmax=u3.mult(s,n,w)
# 	rpmax=u3.mult(s,np,wp)
# 	for r in range(1,rmax+1):
# 		psi=SpStateReduced(s,n,r,w)
# 		for rp in range(1,rpmax+1):
# 			psip=SpStateReduced(s,np,rp,wp)
# 			print psip, psi, vcs.A(psip,psi,1), vcs.A_asymptotic(psip,psi,1), vcs.A_u3boson(psip,psi)




# if __name__ == '__main__':
# 	print	u3.W(SU3State(0,0),1,0,SU3State(2,0),1,0,SU3State(2,0),1,0,1)
# 	print 	u3.W(SU3State(0,0),1,0,SU3State(2,0),1,2,SU3State(2,0),1,2,1)
# 	print 	so3.coupling_9j(
# 		(0,.5,.5),
# 		(2,.5,1.5),
# 		(2,0,2)
# 		)**2
# 	print 	so3.coupling_9j((0,.5,.5),(2,.5,2.5),(2,0,2))**2
# 	print 	so3.coupling_9j(
# 		(0,.5,.5),
# 		(0,.5,.5),
# 		(0,0,0))
# 	print 	so3.coupling_9j((2,0,2),(0,0,0),(2,0,2))
# 	print 

