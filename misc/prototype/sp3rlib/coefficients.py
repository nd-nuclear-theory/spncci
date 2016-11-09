import math
from u3states import *
import so3
import u3
import u3boson
import utilities as utils

Verbose=True
#######################################################################################################
##Coefficients 
#######################################################################################################
### Caching dictionaries 
global Cdict
Cdict={}

global bdict
bdict={}

global Fdict
Fdict={}
#######################################################################################################
def bcoef(n1,n2,n3,r):
	"""
	Calculates the B coefficient needed for the coupling of two symplectic raising polynomials

	Args: 
		n1, n2, n3: U3State class variable
		r: integers, coupling multiplicity of n1xn2->n3

	Returns:
		b coefficient.  For example,
			n1=U3State([2,2,0])
			n2=U3State([2,2,0])
			n3=U3State([4,2,2])	
			bcoef(n1,n2,n3,1)
		returns 
			1.63299316186
	"""
#Should probably check that the elements in the su3lib table are correct if we get multiplicity
	if (n1,n2,n3,r) in bdict:
		coef=bdict[n1,n2,n3,r]
	else:
		nset=[]
		coef=0.0
		n0=U3State([0,0,0])
		#I'm not sure why I had if n3p==[0,0,0], so I'm eliminating it for now. 
		#if n1 or n2 is 0(0,0) then coefficient is trivial
		if n1==n0 or n2==n0:
			coef=1
		else: 
			rmax=u3.mult(n1,n2,n3)
			m1=n3.one()
			m2=n3.two()
			m3=n3.three()
			if (m1-2)>=m2:
				nset.append(U3State([m1-2,m2,m3]))
			if (m2-2)>=m3:
				nset.append(U3State([m1,m2-2,m3]))
			if (m3-2)>=0:
				nset.append(U3State([m1,m2,m3-2]))
			t1=n1.one()
			t2=n1.two()
			t3=n1.three()
			if (t1-2)>=t2:
				n1p=U3State([t1-2,t2,t3])
			else:
				if (t2-2)>=t3:
					n1p=U3State([t1,t2-2,t3])
				else:
					if (t3-2)>=0:
						n1p=U3State([t1,t2,t3-2])
					else:
						print("n1 not valid")
						n1p=0
			for np in nset:
				rpmax=u3.mult(n1p,n2,np)
				for rp in range(1,rpmax+1):
					#print SU3State(2,0),n1p,n3,n2,n1,1,r,np,rp,1
					coef=coef+u3boson.A_n(n3,np)*bcoef(n1p,n2,np,rp)*u3.U(SU3State(2,0),n1p,n3,n2,n1,1,r,np,rp,1)
			coef=coef/u3boson.A_n(n1,n1p)
		bdict[n1,n2,n3,r]=coef
	return coef

########################################################################################################
def ccoef(n,q,np):
	"""
	Calculates the C coefficients in the expansion of the raising polynomial into a single jacobi particle and remaining jacobi particles 
	nbx(q,0)->n
	Args: 
		q: even int
		n, nb: U3State or SU3State class variables 

	Returns:
		Coef, for example 
			CCoef(10,6,2,2,4,2)
		returns
			1.96396101212
	"""
	if (n,q,np) in Cdict:
		Coef=Cdict[n,q,np]
	else:
	####################################################################
		N=n.N
		lbdaN=n.su3.lbda
		muN=n.su3.mu
		Coef=0.0
		## base case
		if q==0:
			if n==np:
				Coef=1
		##recursion relation
		else:
			## special case
			if n.su3==SU3State(N,0): 
				if SU3State(N-q,0)==np.su3:
					Coef=math.sqrt(math.factorial(N/2)*math.factorial(q)/(math.factorial(N/2-q/2)*2**(q/2)))/math.factorial(q/2)
			## recurance
			else:
				[m1,m2,m3]=n.cartesian()
				n2list=[]
				if (m1-2)>=m2:
					nb=U3State([m1-2,m2,m3])
				elif (m2-2)>=m3:
					nb=U3State([m1,m2-2,m3])
				elif (m3-2)>=0:
					nb=U3State([m1,m2,m3-2])
				else:
					return 0.0
				coef1=ccoef(nb,q-2,np)*u3.U(SU3State(2,0),SU3State(q-2,0),n,np,SU3State(q,0),1,1,nb,1,1)*math.sqrt(q*(q-1)/2.)
				## coef2
				coef2=0
				[m1,m2,m3]=np.cartesian()
				## obtain list of possible n2's
				if (m1-2)>=m2:
					n2list.append(U3State([m1-2,m2,m3]))
				if (m2-2)>=m3:
					n2list.append(U3State([m1,m2-2,m3]))
				if (m3-2)>=0:
					n2list.append(U3State([m1,m2,m3-2]))
				## sum over possible n2's
				for n2 in n2list:
					coef2=coef2+u3boson.A_n(np,n2)*ccoef(nb,q,n2)*u3.Z(SU3State(q,0),n2,n,SU3State(2,0),nb,1,1,np,1,1)
				#print n,nb
				Coef=(coef1+coef2)/u3boson.A_n(n,nb)
			#################################################################
		Cdict[n,q,np]=Coef
	return Coef


#######################################################################################################
def F(w0,qbp,qb,wt,rt,tb,t,w0p,qp,q):
	"""
	Calculates the F coefficient for full calculation
	rbp, rb come from initila unit tensor 
"""
	f1=0
	f2=0
	if (qp,q)==(qbp+tb,qb-t):
		f1=(
			math.sqrt(
				utils.choose(qp,tb)
				*utils.choose(qb,t)
				#(1.
				#*math.factorial(qp)
				#*math.factorial(qb)
				#/math.factorial(qbp)
				#/math.factorial(tb)
				#/math.factorial(q)
				#/math.factorial(t)
				)

			*u3.coupling_9lm(
				SU3State(0,tb),SU3State(0,t), wt.conj(),1,
				SU3State(qp,0),SU3State(0,q), w0p,1,
				SU3State(qbp,0),SU3State(0,qb), w0, 1,
				1,1,rt
				)
			)
	if ((qp,q)==(qbp+tb+t,qb)) and (wt.su3==SU3State(t+tb,0)):
		f2=(
			math.sqrt(1.
				*math.factorial(qp)
				/math.factorial(t)
				/math.factorial(tb)
				/math.factorial(qbp)
			)
			*u3.U(SU3State(0,t+tb),SU3State(qp,0),w0,SU3State(0,qb),SU3State(qbp,0),1,1,w0p,1,1)
			)
	fcoef=f1-f2
	return fcoef



if __name__ == '__main__':
	print "F"
	print F(SU3State(0,0),0,0,SU3State(0,0),1,0,0,SU3State(0,0),0,0)
	print F(SU3State(2,2),2,2,SU3State(2,0),1,2,0,SU3State(2,0),2,0)

#######################################################################################################
def mosh_SU3(w1,w2,wr,wc,w):
	mosh=(
			#utils.d(w.lbda/2., (wr.lbda-wc.lbda)/2.,(w1.lbda-w2.lbda)/2.)
			utils.d(w.lbda/2.,(w1.lbda-w2.lbda)/2.,(wr.lbda-wc.lbda)/2.)
			#*utils.parity((w1.pt()-w2.pt())/2.)
			)
	return mosh

#####fishy

#######################################################################################################
def mosh_SO3((n1,L1),(n2,L2),(nr,Lr),(nc,Lc),L):
	mosh=0
	##defining the SU(3) content of states
	w1=SU3State(2*n1+L1,0)
	w2=SU3State(2*n2+L2,0)
	wr=SU3State(2*nr+Lr,0)
	wc=SU3State(2*nc+Lc,0)
	# generating list of w's
	wset=u3.couple(w1,w2)
	## summing over w in wset
	for w in wset:
		## check that w is also in {wrxwc}
		if u3.mult(wr,wc,w)==0:
			continue 
		## calculate max branching multiplicity of L in w
		kmax=so3.multk(w,L)
		## summing over branching multplicity from 1 to kmax
		for k in range(1,kmax+1):
			## calculating moshinksy coefficient
			mosh=mosh+(
				u3.W(wr,1,Lr,wc,1,Lc,w,k,L,1)
				*u3.W(w1,1,L1,w2,1,L2,w,k,L,1)
				*mosh_SU3(w1,w2,wr,wc,w)
				*utils.parity(n1+n2+nr+nc)
				)
	return mosh
#######################################################################################################

def moshinsky(state1,state2,stater,statec,state,type):
	"""
	calculates the moshinsky coefficients used to transform between the two-particle and relative-center of mass frames.  
	w1 x w2->w
	wr x wc->w
	type can be either "SU3" or "SO3"
	if type=="SU3":
		state1,state2,stater,statec,state are simply SU3State class 
	if type=="SO3" then state1 etc are (n1,L1) where N1=2*n1+L1 and w1=SU3State(N1,0)
					state=L
	"""
	mosh=0
	if type=="SU3":
		#check SU3 coupling
		if u3.mult(state1,state2,state)*u3.mult(stater,statec,state)!=0:
			mosh=mosh_SU3(state1,state2,stater,statec,state)
	if type=="SO3":
		#check that angular momentum coupling is allowed and that N1+N2=Nr+Nc
		(n1,L1)=state1
		(n2,L2)=state2
		(nr,Lr)=stater
		(nc,Lc)=statec
		L=state
		if (so3.mult(L1,L2,L)!=0) and (so3.mult(Lc,Lr,L)!=0) and ((2*n1+L1+2*n2+L2)==(2*nr+Lr+2*nc+Lc)):
			mosh=mosh_SO3((n1,L1),(n2,L2),(nr,Lr),(nc,Lc),L)
	return mosh 


if __name__ == '__main__':
	## Testing SO(3) branched moshinsky coefficients against MAC's code.  
	filename="moshinsky_bracket_table_Nmax6.out"
	moshinsky_file=open(filename,'r')
	for line in moshinsky_file:
		values=line.split()
		[Nr,Lr,Nc,Lc,N1,L1,N2,L2,L]=map(int,values[:-1])
		mosh_MAC=float(values[-1])

		mosh_AEM=moshinsky(((N1-L1)/2,L1),((N2-L2)/2,L2),((Nr-Lr)/2,Lr),((Nc-Lc)/2,Lc),L,"SO3")
 		print "{}{}{}{}{}  {:4.10f}".format((Nr,Lr),(Nc,Lc),(N1,L1),(N2,L2),L,abs(mosh_AEM-mosh_MAC))







# if __name__ == '__main__' and (Verbose==True):

# 	print "moshinsky checked against Table of Transformation Brackets, T.A. Brody and M. Moshinsky"
# 	(n1,n2)=(1,0)
# 	(L1,L2,L)=(0,1,1)
# 	(nr,Lr,nc,Lc)=(0,1,1,0)
# 	listn=[
# 	(0,0,0,0,0,0,0,0,0,1),
# 	(0,0,0,4,4,0,0,0,4,.25),
# 	(0,0,0,4,4,0,1,0,3,-.5),
# 	(0,0,0,4,4,0,2,0,2,.61237245),
# 	(0,0,0,4,4,0,3,0,1,-.5),
# 	(0,0,0,4,4,0,4,0,0,.25),
# 	(0,0,2,5,4,0,2,0,5,.1178511),
# 	(0,0,2,5,4,0,2,1,3,-.39086799),
# 	(0,0,2,5,4,0,3,0,4,-.05976141),
# 	(0,0,2,5,4,0,3,1,2,.44320265),
# 	(0,0,2,5,4,0,4,0,3,-.03976142),
# 	(0,0,2,5,4,0,4,1,1,-.30276504),
# 	(0,0,2,5,4,0,5,0,2,.11785114),
# 	(0,0,2,5,4,1,1,0,4,-.30275503),
# 	(0,0,2,5,4,1,2,0,3,.44320263),
# 	(0,0,2,5,4,1,3,0,2,-.39086801),
# 	(0,0,2,5,4,1,4,0,1,.20412414),
# 	(1,0,1,6,5,0,6,0,3,-.07364478),
# 	(1,0,1,6,5,0,6,1,1,-.18481203),
# 	(1,0,1,6,5,1,5,0,2,-.17496462),
# 	(1,0,1,6,5,2,2,0,3,-.12580908),
# 	(1,0,1,6,5,2,5,0,0,.11306677),
# 	(0,1,1,2,1,0,0,2,1,.28867514),
# 	(0,1,1,2,1,0,1,2,0,-.44095855),
# 	(0,1,1,2,1,1,1,1,0,0.0),
# 	(0,1,1,2,1,2,1,0,0,.28867514),
# 	(0,1,1,2,1,0,3,0,2,-.18257419),
# 	(2,0,0,4,4,0,1,2,3,-.17840601),
# 	(2,0,0,4,4,0,2,2,2,-.00498010),
# 	(2,0,0,4,4,0,3,0,5,-.23821426)
# 	]
# 	for n1,n2,L1,L2,L,nr,Lr,nc,Lc,tabmosh in listn:
# 		print "{}{}{}{}{}{:18.10f}{:18.10f}".format((n1,L1),(n2,L2),(nr,Lr),(nc,Lc),L,tabmosh,moshinsky((n1,L1),(n2,L2),(nr,Lr),(nc,Lc),L,"SO3"))
# 	#################################################################################################################
	
# 	print "bcoef checked against npa-455-1986-315-Suzuki"
# 	print bcoef(U3State([2,0,0]), U3State([2,0,0]),U3State([4,0,0]),1), math.sqrt(2)
# 	print bcoef(U3State([2,0,0]), U3State([2,0,0]),U3State([2,2,0]),1), math.sqrt(2)
# 	print bcoef(U3State([4,0,0]), U3State([2,0,0]),U3State([6,0,0]),1), math.sqrt(3)
# 	print bcoef(U3State([4,0,0]), U3State([2,0,0]),U3State([4,2,0]),1), math.sqrt(4./3)
# 	print bcoef(U3State([2,2,0]), U3State([2,0,0]),U3State([4,2,0]),1), math.sqrt(5./3)
# 	print bcoef(U3State([2,2,0]), U3State([2,0,0]),U3State([2,2,2]),1), math.sqrt(3)
# 	print bcoef(U3State([6,0,0]), U3State([2,0,0]),U3State([8,0,0]),1), math.sqrt(4)
# 	print bcoef(U3State([6,0,0]), U3State([2,0,0]),U3State([6,2,0]),1), math.sqrt(6./5)
# 	print bcoef(U3State([4,2,0]), U3State([2,0,0]),U3State([6,2,0]),1), math.sqrt(14./5)
# 	print bcoef(U3State([4,2,0]), U3State([2,0,0]),U3State([4,2,2]),1), math.sqrt(5./2)
# 	print bcoef(U3State([4,2,0]), U3State([2,0,0]),U3State([4,4,0]),1), math.sqrt(4)
# 	print bcoef(U3State([4,0,0]), U3State([2,2,0]),U3State([6,2,0]),1), math.sqrt(7./3)
# 	print bcoef(U3State([4,0,0]), U3State([2,2,0]),U3State([4,2,2]),1), math.sqrt(5./3)
# 	print bcoef(U3State([4,0,0]), U3State([2,2,0]),U3State([4,4,0]),1), math.sqrt(0)
# 	print bcoef(U3State([2,2,0]), U3State([2,2,0]),U3State([4,2,2]),1), math.sqrt(8./3)
# 	print bcoef(U3State([2,2,0]), U3State([2,2,0]),U3State([4,4,0]),1), math.sqrt(10./3)
# 	print 
# 	#print "hi"
# 	#print bcoef(U3State(0,SU3State(0,0)),U3State(8,SU3State(8,0)),U3State(8,SU3State(8,0)),1)
# 	#(2,0) 0(0,0) 10(6,2) 8(8,0) 2(2,0) 1 1 8(8,0) 1 1
# 	#################################################################################################################
# 	print "ccoef, checked against npa-455-1986-315-Suzuki"
# 	print ccoef(U3State(4,SU3State(4,0)),4,U3State(0,SU3State(0,0))), math.sqrt(3)
# 	print ccoef(U3State(4,SU3State(4,0)),2,U3State(2,SU3State(2,0))), math.sqrt(2)
# 	print ccoef(U3State(4,SU3State(0,2)),2,U3State(2,SU3State(2,0))), math.sqrt(2)
# 	print ccoef(U3State(6,SU3State(6,0)),6,U3State(0,SU3State(0,0))), math.sqrt(15)
# 	print ccoef(U3State(6,SU3State(6,0)),4,U3State(2,SU3State(2,0))), math.sqrt(9)
# 	print ccoef(U3State(6,SU3State(6,0)),2,U3State(4,SU3State(4,0))), math.sqrt(3)
# 	print ccoef(U3State(6,SU3State(2,2)),4,U3State(2,SU3State(2,0))), math.sqrt(4)
# 	print ccoef(U3State(6,SU3State(2,2)),2,U3State(4,SU3State(4,0))), math.sqrt(4./3)
# 	print ccoef(U3State(6,SU3State(2,2)),2,U3State(4,SU3State(0,2))), math.sqrt(5./3)
# 	print ccoef(U3State(6,SU3State(0,0)),2,U3State(4,SU3State(0,2))), math.sqrt(3)
# 	print ccoef(U3State(8,SU3State(8,0)),8,U3State(8,SU3State(0,0))), math.sqrt(105)
# 	print ccoef(U3State(8,SU3State(8,0)),6,U3State(2,SU3State(2,0))), math.sqrt(4*15)
# 	print ccoef(U3State(8,SU3State(8,0)),4,U3State(4,SU3State(4,0))), math.sqrt(18)
# 	print ccoef(U3State(8,SU3State(8,0)),4,U3State(4,SU3State(0,2))), math.sqrt(0)
# 	print ccoef(U3State(8,SU3State(8,0)),2,U3State(6,SU3State(6,0))), math.sqrt(4)
# 	print ccoef(U3State(8,SU3State(4,2)),6,U3State(2,SU3State(2,0))), math.sqrt(18)
# 	print ccoef(U3State(8,SU3State(4,2)),4,U3State(4,SU3State(4,0))), math.sqrt(4)
# 	print ccoef(U3State(8,SU3State(4,2)),4,U3State(4,SU3State(0,2))), math.sqrt(7)
# 	print ccoef(U3State(8,SU3State(4,2)),2,U3State(6,SU3State(6,0))), math.sqrt(6./5)
# 	print ccoef(U3State(8,SU3State(4,2)),2,U3State(6,SU3State(2,2))), math.sqrt(14./5)
# 	print ccoef(U3State(8,SU3State(0,4)),4,U3State(4,SU3State(4,0))), math.sqrt(8)
# 	print ccoef(U3State(8,SU3State(0,4)),2,U3State(6,SU3State(2,2))), math.sqrt(4)
# 	print ccoef(U3State(8,SU3State(2,0)),4,U3State(4,SU3State(0,2))), math.sqrt(5)
# 	print ccoef(U3State(8,SU3State(2,0)),2,U3State(6,SU3State(2,2))), math.sqrt(5./2)
# 	print ccoef(U3State(8,SU3State(2,0)),2,U3State(6,SU3State(0,0))), math.sqrt(3./2)


# # sets=u3.couple(SU3State(76,0),SU3State(2,0))
# # for w in sets:
# # 	print w.pr(), sets[w]

# # Ccoef_dic={}
# # Nmax=20
# # polynomials=u3boson.poly(Nmax)
# # for N in range(0,Nmax+1,2):
# # 	nset=polynomials[N]
# # 	for n in nset:
# # 		for q in range(0,N+1,2):
# # 			nbset=u3.couple(n,SU3State(0,q))
# # 			for lmnb in nbset:
# # 				nb=U3State(N-q,lmnb)
# # 				if nb not in polynomials[N-q]:
# # 					continue
# # 				coef=ccoef(n,q,nb)
# # 				Ccoef_dic[n,q,nb]=coef
# # print 
# # print len(Ccoef_dic)
# # print sys.getsizeof((n,q,nb))*len(Ccoef_dic)
# ## at 20 50 kB





