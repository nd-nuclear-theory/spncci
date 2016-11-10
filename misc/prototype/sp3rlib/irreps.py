##################################################################################################################################################
##  Generating irreps corresponding to the different symmetries 
##################################################################################################################################################

import math
import sys
import os
import numpy
sys.path.append('/Users/annamccoy/Dropbox/Research/Sp(6,R)/Code/')
from u3states import *
from particle_states import *
from spstates import *
import u3
import so3
import u3boson
import utilities as utils

##################################################################################################################################################
def Sp3RxU2ReducedIrrep(s,S,Nmax):
	""""
	Constructs a canonically sorted list of Sp(3,R)xSU(3) states up to Nmax
	  for given irrep s

	Args: 
		s: symplectic bandhead (U3State)
		S: spin (int)

	Returns:
		[state1,state2 etc]
	"""
	## initializing list
	states=[]	
	## generating dictionary of all possible raising polynomial labels 
	## dictionary of form {N:[n1,n2,..]}
	PolyDict=u3boson.poly(Nmax)
	## iterating over N from 0 to Nmax
	for N in range(0,Nmax+1,2):
		## selecting polynomial labels for each N
		polylist=PolyDict[N]
		## order polynomial labels cannonically 
		polylist_canonical=u3.canonical_sort(polylist)
		## iterating over n in polylist
		for n in polylist_canonical:
			## coupling raising polynomials to bandhead
			CoupledDict=u3.couple(s,n)
			## iterating over generated list of (lbda,mu)
			wlist_canonical=u3.canonical_sort(CoupledDict.keys())
			for lambdamu in wlist_canonical:
				## converting from SU(2) to U(3)
				w=U3State(s.N+N,lambdamu)
				## number of times w occurs in kronecker product sxn
				rmax=CoupledDict[lambdamu]
				## iterating over occurances of w
				for r in range(1,rmax+1):
					## appending SpxU2 state to state list 
					state=SpU2StateRed(SpStateReduced(s,n,r,w),S)
					states.append(state)
	return states



##################################################################################################################################################

def Sp3RxU2ReducedIrrepN(s,S,Nmax):
	"""
	Generates list of states in irrep (s,S) up to Nmax for SxL->J 

	Returns:
		statelist: list of states (SpU2StateJ).  List is ordered according to
					
	"""
	## initializing list of states 
	stateNdic={}
	## generate list of reduced states in canonical ordering  
	redStateList=Sp3RxU2ReducedIrrep(s,S,Nmax)
	## iterating over reduced states 
	for redState in redStateList:
		N=redState.sp.w.N
		if N in stateNdic:
			stateNdic[N]=stateNdic[N]+[redState]
		else:
			stateNdic[N]=[redState]

	return stateNdic





##################################################################################################################################################
def Sp3RxU2Irrep(s,S,J,Nmax):
	"""
	Generates list of states in irrep (s,S) up to Nmax for SxL->J 

	Returns:
		statelist: list of states (SpU2StateJ)
	"""
	## initializing list of states 
	statelist=[]
	## generate list of reduced states 
	redStateList=Sp3RxU2ReducedIrrep(s,S,Nmax)
	## iterating over reduced states 
	for redState in redStateList:
		## extracting SpStateReduced state and Spin 
		stateSp=redState.sp
		S=redState.S
		## branch to L
		Lkset=so3.branching_so3(stateSp.w)
		## sorting L's from smallest to largest
		Lset=sorted(Lkset.keys())
		## iterating over L
		for L in Lset:
			if so3.mult(S,L,J)==0:
				continue
			kmax=Lkset[L]
			for k in range(1,kmax+1):
				statelist.append(SpU2StateJ(stateSp,k,L,S,J))
	return statelist
##################################################################################################################################################
def Sp3RxU2IrrepN(s,S,J,Nmax):
	"""
	Generates list of states in irrep (s,S) up to Nmax for SxL->J 

	Returns:
		statelist: list of states (SpU2StateJ).  List is ordered according to
					
	"""
	## initializing list of states 
	stateNdic={}
	## generate list of reduced states in canonical ordering  
	redStateList=Sp3RxU2ReducedIrrep(s,S,Nmax)
	## iterating over reduced states 
	for redState in redStateList:
		## extracting SpStateReduced state and Spin 
		stateSp=redState.sp
		S=redState.S
		## branch to L
		Lkset=so3.branching_so3(stateSp.w)
		## sorting L's from smallest to largest
		Lset=sorted(Lkset.keys())
		## iterating over L
		for L in Lset:
			if so3.mult(S,L,J)==0:
				continue
			kmax=Lkset[L]
			for k in range(1,kmax+1):
				if stateSp.w.N in stateNdic:
					stateNdic[stateSp.w.N]=stateNdic[stateSp.w.N]+[(SpU2StateJ(stateSp,k,L,S,J))]
				else:
					stateNdic[stateSp.w.N]=[SpU2StateJ(stateSp,k,L,S,J)]

	return stateNdic
##################################################################################################################################################
s=U3State(24,SU3State(0,0))
S=0
J=0
Nmax=4
#print Sp3RxU2IrrepN(s,S,J,Nmax)




def Sp3RxU4ReducedIrrep(s,S,T,Nmax):
	""""
	Constructs a canonically sorted list of Sp(3,R)xSU(3) states up to Nmax
	  for given irrep s

	Args: 
		s: symplectic bandhead (U3State)
		S: spin (int)

	Returns:
		[state1,state2 etc]
	"""
	## initializing list
	states=[]	
	## generating dictionary of all possible raising polynomial labels 
	## dictionary of form {N:[n1,n2,..]}
	PolyDict=u3boson.poly(Nmax)
	## iterating over N from 0 to Nmax
	for N in range(0,Nmax+1,2):
		## selecting polynomial labels for each N
		polylist=PolyDict[N]
		## order polynomial labels cannonically 
		polylist_canonical=u3.canonical_sort(polylist)
		## iterating over n in polylist
		for n in polylist_canonical:
			## coupling raising polynomials to bandhead
			CoupledDict=u3.couple(s,n)
			## iterating over generated list of (lbda,mu)
			wlist_canonical=u3.canonical_sort(CoupledDict.keys())
			for lambdamu in wlist_canonical:
				## converting from SU(2) to U(3)
				w=U3State(s.N+N,lambdamu)
				## number of times w occurs in kronecker product sxn
				rmax=CoupledDict[lambdamu]
				## iterating over occurances of w
				for r in range(1,rmax+1):
					## appending SpxU2 state to state list 
					state=SpU4StateRed(SpStateReduced(s,n,r,w),S,T)
					states.append(state)
	return states

##################################################################################################################################################
def Sp3RxU4Irrep(s,S,J,Nmax):
	"""
	Generates list of states in irrep (s,S) up to Nmax for SxL->J 

	Returns:
		statelist: list of states (SpU2StateJ)
	"""
	## initializing list of states 
	statelist=[]
	##The isospin symmetry is determined by the constraint (J+T)%2=1
	T=(J+1)%2
	## generate list of reduced states 
	redStateList=Sp3RxU4ReducedIrrep(s,S,T,Nmax)
	## iterating over reduced states 
	for redState in redStateList:
		## extracting SpStateReduced state and Spin 
		stateSp=redState.sp
		S=redState.S
		T=redState.T
		## branch to L
		Lkset=so3.branching_so3(stateSp.w)
		## sorting L's from smallest to largest
		Lset=sorted(Lkset.keys())
		## iterating over L
		for L in Lset:
			if so3.mult(S,L,J)==0:
				continue
			if (L+S+T)%2==0:
				continue
			kmax=Lkset[L]
			for k in range(1,kmax+1):
				statelist.append(SpU4StateJ(stateSp,k,L,S,J,T))
	return statelist


##################################################################################################################################################
def singleState_set(Nmax):
	"""
	Returns a list of single particle states. 
	 Note that here N is not the radial quantum number n, its N=2n+L
	"""
	Stateset=[]
	for N in range(Nmax+1):
		for L in range(N%2,N+1,2):
			jset=so3.couple_so3(L,.5)
			jset.sort()
			for J in jset:
				Stateset.append(SingleState(N,L,J))
	return Stateset
##################################################################################################################################################
def twobody_states(Nmax1, Nmax2, J,T):
	"""
	Returns a list of two-body states in lexographical order
	 which couple to final angular momentum J
	 H2 format

	Args:
		Nmax1 is the 1Body cutoff (integer)
		Nmax2 is the 2Body cutoff (integer)
		Jmax is the J cutoff (integer or half integer as float)
		Type is [11] for proton-proton, [12] for proton-neutron and [22] for neutron-neutron
	"""
	Twobody=[]
	#Generate a list of single particle states up to Nmax cutoff
	Single_set=singleState_set(Nmax1)
	#Loop over all pairs of single particle states
	for a1 in Single_set:
		for a2 in Single_set:
			#check that N1+N2<N2max 
			if (a1.N+a2.N)>Nmax2:
				continue
			#check for canonical ordering 
			if a1>a2:
				continue
			#check that states can couple to total J
			if so3.mult(a1.J,a2.J,J)==0:
				continue
			#define 2Body state
			state=TBStateJT(a1,a2,J,T)
			#check that 2Body state is symmetry allowed
			if state.is_symmetry_allowed():
				#if symmetry allowed, append to list of 2Body states
				Twobody.append(state)
	return Twobody

##################################################################################################################################################
def irrep_index(s,S,J,Nmax,Nmin=0,indx=0):
	"""
	constructs the simplectic irrep starting with lowest grade U(3) states s

	inputs: 
		s: Lowest grade U(3) irrep, U3 class
		S: Spin, integer or half integer
		J: angular momentum: integer or half integer
		Nmax: max number of excitations, integer
		Nmin: Starting excitation level, integer, optional with defaut Nmin=0
		index: Starting indexing number, option, default is indx=0
	output:
		irrep: Dictionary {N: (n,r,w,k,L)}
		index: Dictionary {(n,r,w,k,L):i}
			where n: Raising polynomial U3 content, class U3
				  w: U3 content of states in Sp irrep, class U3
				  r: Multiplicity index of w in sxn->w, integer>0
				  k: branching multiplicity SU(3)->SO(3), integer>0
				  L: Oribtal angular momentum, integer
	"""
	###initializing dictionarys 
	irrep={}
	index={}
	SO3branch={}
	i=indx # iterator for indexing the basis states
	###intial calculations 
	#generating raising polynomials
	polyset=u3boson.poly(Nmax)
	#generating Symplectic basis states and indexing states
	for N in range(Nmin,Nmax+1,2):
		nset=polyset[N] #generating raising polynomails
		for n in nset:
			wset=u3.couple(s,n)
			for lm in wset:
				if lm not in SO3branch:
					LKset=so3.branching_so3(lm,S,J,Restrict=True)#SO3branching_J(lm,mu,lset)
					SO3branch[lm]=LKset
				#generatoring omega and rhomax
				else:
					LKset=SO3branch[lm]
				##########
				if LKset!={}:
					w=U3State(s.N+N,lm)
					rmax=wset[lm]
					##generating index for matrix states  starting with i=0
					for r in range(1,rmax+1):
						for L in LKset:
							for k in range(1,LKset[L]+1):
								psi=SpState(s,n,r,w,k,L)
								##generating nested dictionary irrep with the form {N:(n,rmax,w,k,L,M)}
								if N in irrep:
									irrep[N].append(psi)
								else:
									irrep[N]=[psi]
								index[psi]=i
								i=i+1
	#print "dim", len(index)
	return irrep, index

def relwST_irrep(Nmax):
	"""
	Returns a list of relative states up to Nmax for all possible S and T
	"""
	stateset=[]
	for Nr in range(0,Nmax+1):
		wr=SU3State(Nr,0)
		for T in [0,1]:
			for S in [0,1]:
				# if Nr%2!=(T+S+1)%2:
				# 	continue
				state=U3SU4State(wr,S,T)
				stateset=stateset+[state]
	return stateset


def relstateT_set(Nmax):
	"""
	Returns a list of relative states up to Nmax for all possible S and T
	"""
	stateset=[]
	for Nr in range(0,Nmax+1):
		wr=SU3State(Nr,0)
		for T in [0,1]:
			for S in [0,1]:
				if (Nr+T+S)%2!=1:
					continue
				state=U3SU4State(wr,S,T)
				stateset=stateset+[state]
	return stateset






if __name__ == '__main__':
	pass
# 	################################################################################################
# 	##  bandheads={(lbda,mu):[(S,num)]}	
# 	##  Ns is the number of oscillator quanta in the Nex=0 space 
# 	##  Nspace is max excitation for bandheads, Nspace=2 would mean bandheads in the Nex=0 and Nex=2 spaces
# 	################################################################################################
# 	Z=3
# 	NN=3
# 	bandheads={(0,1):[(0,1),(1,1)],(2,0):[(0,1),(1,1)]}
# 	Ns=11
# 	################################################################################################
# 	Z=4
# 	NN=4
# 	bandheads={(4,0):[(0,1)],(2,1):[(0,2),(1,2)],(0,2):[(0,1),(1,1),(2,1)],(1,0):[(0,1),(1,4),(2,1)]}
# 	Ns=16
# 	#################################################################################################
# 	Z=4
# 	NN=5
# 	bandheads={
# 		(3,1):[(.5,2)],
# 		(1,2):[(.5,4),(1.5,3)],
# 		(2,0):[(.5,2),(1.5,4)],
# 		(0,1):[(.5,3),(1.5,2),(2.5,1)]
# 	}
# 	Ns=14.5
# 	################################################################################################
# 	Z=5
# 	NN=5
# 	bandheads={(0,0):[(0,2),(1,2),(2,1),(3,1)],(0,3):[(0,2),(1,2)],(1,1):[(0,1),(1,3),(2,2)],(2,2):[(0,1),(1,1)],(3,0):[(0,1),(1,1)]}
# 	Ns=18
# 	################################################################################################
# 	Z=6
# 	NN=6
# 	bandheads={
# 	(0,4):[(0,1)],
# 	(1,2):[(0,1),(1,1)],
# 	(2,0):[(0,2),(1,1),(2,1)],
# 	(0,1):[(0,1),(1,2),(2,1)],
# 	(1,2):[(1,2)]
# 	}
# 	Ns=20
# 	#################################################################################################
# 	Z=7
# 	NN=7
# 	bandheads={
# 	(0,2):[(0,1),(1,1)],
# 	(1,0):[(0,1),(1,1)],
# 	}
# 	Ns=22
	# #################################################################################################
	Z=8
	NN=8
	bandheads={(0,0):[(0,1)]}
	Ns=24
# 	##################################################################################################
# 	Z=10
# 	NN=10
# 	Ns=50
# 	bandheads={
# 	(0,4):[(0,3),(2,1),(4,1)],
# 	(1,2):[(0,3),(2,6),(4,2)],
# 	(2,0):[(0,4),(2,3),(4,1)],
# 	(2,3):[(0,2),(2,5),(4,1)],
# 	(3,1):[(0,4),(2,6),(4,2)],
# 	(4,2):[(0,4),(2,3)],
# 	(6,1):[(0,1),(2,2)],
# 	(8,0):[(0,1)],
# 	(0,1):[(0,1),(2,3),(4,1)],
# 	(5,0):[(0,1),(2,3),(4,1)]
# 	}
# 	##################################################################################################
# 	Z=12
# 	NN=12
# 	Ns=64
# 	bandheads={
# 	(0,2):[(0,24),(2,35),(4,27),(6,9),(8,9)],
# 	(0,5):[(0,13),(2,31),(4,15),(6,5),(8,5)],
# 	(0,8):[(0,3),(2,1),(4,1)],
# 	(1,0):[(0,10),(2,26),(4,16),(6,5),(8,5)],
# 	(1,3):[(0,40),(2,74),(4,46),(6,14),(8,14)],
# 	(1,6):[(0,11),(2,20),(4,10),(6,2),(8,2)],
# 	(2,1):[(0,36),(2,75),(4,49),(6,16),(8,16)],
# 	(2,4):[(0,40),(2,61),(4,37),(6,9),(8,9)],
# 	(2,7):[(0,2),(2,5),(4,1)],
# 	(3,2):[(0,46),(2,91),(4,53),(6,13),(8,13)],
# 	(3,5):[(0,16),(2,26),(4,12),(6,2),(8,2)],
# 	(4,0):[(0,30),(2,44),(4,30),(6,9),(8,9)],
# 	(4,3):[(0,32),(2,60),(4,30),(6,6),(8,6)],
# 	(4,6):[(0,4),(2,3),(4,1)],
# 	(5,1):[(0,32),(2,56),(4,30),(6,6),(8,6)],
# 	(5,4):[(0,11),(2,31),(4,7)],
# 	(6,2):[(0,23),(2,9),(4,15),(6,2),(8,2)],
# 	(6,5):[(0,1),(2,2)],
# 	(7,0):[(0,8),(2,18),(4,6)],
# 	(7,3):[(0,6),(2,8),(4,2)],
# 	(8,1):[(0,7),(2,11),(4,3)],
# 	(8,4):[(0,1)],
# 	(9,2):[(0,1),(2,2)],
# 	(10,0):[(0,2),(2,1),(4,1)],
# 	}




	#################################################################################################
	Nmax=30
	Jmax=10
	Nspace=0
	#################################################################################################
	# #################################################################################################
	#file1=open('Dimensions/dimJp_{}_{}_{}_{}_{}'.format(Z,NN,Nspace,Jmax,Nmax),'a')
	file2=open('Dimensions/dim-{}-{}.dat'.format(Z,NN),'a')
	dimtotal=numpy.zeros(Nmax/2+1)
	for J in range(0,Jmax+1):
		#file3=open('Dimensions/DimJMg-{}.dat'.format(J),'a')
		#file1.write("J={}\n".format(J))
		Dims=numpy.zeros(Nmax/2+1)
		Dimeven=numpy.zeros(Nmax/2+1)
		Dimodd=numpy.zeros(Nmax/2+1)
		
		for (l,m) in bandheads:
			s=U3State(Ns,SU3State(l,m))
			Sset=bandheads[l,m]
			for (S,num) in Sset:
				#print "{} {:3}".format(s,S)
				irrep0,index0=irrep_index(s,S,J,Nmax)
				if irrep0=={}:
					continue

				for N in range(0,Nmax+1,2):
					e=0
					o=0
					if N not in irrep0:
						continue
					for psi in irrep0[N]:
						if psi.L%2==0:
							e=e+1
						else:
							o=o+1
					Ndimeven=num*e
					Ndimodd=num*o
					Ndim=num*len(irrep0[N])
					Dimeven[N/2]=Dimeven[N/2]+Ndimeven
					Dimodd[N/2]=Dimodd[N/2]+Ndimodd
					Dims[N/2]=Dims[N/2]+Ndim
		#print J	
		eventotal=0
		oddtotal=0
		for i in range(Nmax/2+1):
			eventotal=eventotal+Dimeven[i]
			oddtotal=oddtotal+Dimodd[i]
			#file3.write("{:2}  {:12}\n".format(2*i,eventotal))
			#file1.write("{:2}  {:12} {:12} {:12} {:12}\n".format(2*i,int(Dimeven[i]),int(eventotal),int(Dimodd[i]),int(oddtotal)))

			dimtotal[i]=dimtotal[i]+eventotal+oddtotal	
		#file3.close()
	for i in range(Nmax/2+1):
		file2.write("{:2}  {:15}\n".format(2*i,int(dimtotal[i])))
			#file1.write("{:2}  {:12} {:12}\n".format(2*i,int(Dims[i])**2,int(Dims[i])**2-zer0MEint(total)**2))
	#file1.close()
	file2.close()
	#			for psi in irrep0[N]:
				# 	print psi











	# file1=open('Dimensions/MEp_{}_{}_{}_{}_{}'.format(Z,NN,Nspace,Jmax,Nmax),'a')
	# for J in range(0,Jmax+1):
	# 	file1.write("J={}\n".format(J))
	# 	Dims=numpy.zeros(Nmax/2+1)
	# 	Dimeven=numpy.zeros(Nmax/2+1)
	# 	Dimodd=numpy.zeros(Nmax/2+1)
		
	# 	for (l,m) in bandheads:
	# 		s=U3State(Ns,SU3State(l,m))
	# 		Sset=bandheads[l,m]
	# 		for (S,num) in Sset:
	# 			#print "{} {:3}".format(s,S)
	# 			irrep0,index0=irrep_index(s,S,J,Nmax)
	# 			if irrep0=={}:
	# 				continue

	# 			for N in range(0,Nmax+1,2):
	# 				e=0
	# 				o=0
	# 				if N not in irrep0:
	# 					continue
	# 				for psi in irrep0[N]:
	# 					if psi.L%2==0:
	# 						e=e+1
	# 					else:
	# 						o=o+1
	# 				Ndimeven=num*e
	# 				Ndimodd=num*o
	# 				Ndim=num*len(irrep0[N])
	# 				Dimeven[N/2]=Dimeven[N/2]+Ndimeven
	# 				Dimodd[N/2]=Dimodd[N/2]+Ndimodd
	# 				Dims[N/2]=Dims[N/2]+Ndim
	# 	#print J	
	# 	eventotal=0
	# 	oddtotal=0
	# 	for i in range(Nmax/2+1):
	# 		eventotal=eventotal+Dimeven[i]
	# 		oddtotal=oddtotal+Dimodd[i]

	# 		file1.write("{:2}  {:12} {:12} {:12} {:12}\n".format(2*i,int(Dimeven[i])**2,int(eventotal)**2,int(Dimodd[i])**2,int(oddtotal)**2))
	# 		#file1.write("{:2}  {:12} {:12}\n".format(2*i,int(Dims[i])**2,int(Dims[i])**2-zer0MEint(total)**2))
	# file1.close()
	# 			# for psi in irrep0[N]:
	# 			# 	print psi


















# # if __name__ == '__main__':
# # 	s=U3State(24,SU3State(0,0))
# # 	S=1
# # 	Nmax=6
# # 	for statered in Sp3RxU2ReducedIrrep(s,S,Nmax):
# # 		print statered
# # 	print 
# # 	for state in Sp3RxU2Irrep(s,S,1,Nmax):
# # 		print state

# # ##################################################################################################################################################
# # if __name__ == '__main__':
# # 	s=U3State(24,SU3State(0,0))
# # 	S=1
# # 	T=1
# # 	Nmax=6
# # 	print 
# # 	for statered in Sp3RxU4ReducedIrrep(s,S,T,Nmax):
# # 		print statered
# # 	print 
# # 	for state in Sp3RxU4Irrep(s,S,2,T,Nmax):
# # 		print state










