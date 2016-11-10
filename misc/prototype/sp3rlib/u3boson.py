"""
Matrix elements in the U3-boson basis and laddering functions 
"""
#######################################################################################################
import numpy
from u3states import *
from spstates import *
import u3
import utilities as utils
#######################################################################################################
def nlist(s,w):
    """
    Creates a dictionary of symplectic multiplicities (n,rmax) 

    Args: 
        s: Symplectic irrep labes (u3State class) 
        w: u3 label of symplectic weight (u3State class) 

    Returns: 
      nl: dictionary
    """
    N=int(w.N-s.N)
    nl={}
    if N%2!=0:
        pass
        #If om is a weight state in the irrep s then their N value must differ by a multiple of 2. 
    else: 
        #this is generating a list of possible n's such that sxn->omega
        if N==2:
            mlt=u3.mult(s,SU3State(2,0),w) 
            if mlt>0:
                n=U3State([2,0,0])
                nl[n]=1    
        else:                                               
            for a in range(0,N/3*2+1,2):                  #N/4 comes from the fact that for n to be a u3 state, n1>=n2>=n3 and even
                for b in range(0,N/3+1,2):
                    if (N-a)>=(a-b)>=b:
                      n=U3State([N-a,a-b,b])
                    #Calculate the multiplicity of (lambda_omega,mu_omega) in the coupling (lambda_s,mu_s)X(lambda_n,mu_n) 
                      mlt=u3.mult(s,n,w) 
                      if mlt>=1:
                          nl[n]=mlt
    return nl

##################################################################################################
def poly(Nmax):
    """
    Creates a dictionary with the excitation levels Nn as keys and the possible U(3) numbers for n.  

    Args:
        Nmax: the max number of excitation.  Or how high you want to ladder in the irrep. 

    Returns:
        Ply: dictionary of the form {Nn:[n]}.  For example, 
        poly(4)
        returns 
            {0:U3State([0,0,0]),2:U3State([2,0,0]), 4:[U3State([4,0,0]), U3State([2,2,0])}
    """
    Ply={}
    Ply={0:[U3State([0,0,0])], 2:[U3State([2,0,0])]}
    for N in range(4,Nmax+1,2):
        nlist=[]
        for a in range(0, 2*N/3+1, 2):
            for b in range(0, N/3+1, 2):
                if (N-a)>=(a-b)>=b:
                    n=U3State([N-a,a-b,b])
                    nlist.append(n)
        Ply[N]=nlist
    return Ply



# ####################################################################
#Functions

def A_n(a,b):
  """
  Calculates (a||a^\dagger||b). 

  Args: 
    a: U(3) final state.  Variable class u3.
    b: U(3) initial state. Variable class u3. 

  Returns:
    The matrix element (a||a^\dagger||b).  For example:  
      np=u3([20,15,10])
      n=u3([20,13,10])
      A_n(np,n)
    returns
        2.78886675511. 
  """
  matrix=0
  [a1,a2,a3]=a.cartesian()
  [b1,b2,b3]=b.cartesian()
  
  if (a1==b1+2 and a2==b2 and a3==b3):
      matrix=numpy.sqrt((b1+4)*(b1-b2+2)*(b1-b3+3)/(2.*(b1-b2+3)*(b1-b3+4)))
  elif (a1==b1 and a2==b2+2 and a3==b3):
      matrix=numpy.sqrt((b2+3)*(b1-b2)*(b2-b3+2)/(2.*(b1-b2-1)*(b2-b3+3)))
  elif(a1==b1 and a2==b2 and a3==b3+2):
      matrix=numpy.sqrt((b3+2)*(b2-b3)*(b1-b3+1)/(2.*(b1-b3)*(b2-b3-1)))
  return matrix

if __name__ == '__main__':
  print "A_n"
  np=U3State([20,15,10])
  n=U3State([20,13,10])
  print A_n(np,n)
########################################################################

def A(psipred, psired,r0=1):
  """
  Calculates (s np rp wp||a^\dagger|s n r w)

  Args:
    s: symplectic bandhead.  u3 class variable. 
    np: U(3) label for raising polynomial of final state. u3 class var. 
    rp: SU(3) outer multiplicity of s x np -> wp.  Integer var. 
    wp: U(3) label for symplectic weight of final state. u3 class var
    n: U(3) label for raising polynomail of initial state. U2 class var. 
    r: SU(3) outer mult. of s x n ->w. Integer var. 
    w: U(3) label for symplectic weight of initial state. u3 class var.

  Returns:
    Matrix element (s np rp wp||a^\dagger|s n r w).  For example 
      wp=u3([22,15,10])
      np=u3([4,0,0])
      rp=1
      w=u3([21,14,10])
      n=u3([2,0,0])
      r=1
      s=u3([20,13,10])
      A_weylboson(s,np,rp,wp,n,r,w)
    returns 
      1.12687233964
  """
  matrix=0.0
  if r0==1: 	 	
    s=psired.s
    n=psired.n
    r=psired.r
    w=psired.w
    np=psipred.n
    rp=psipred.r
    wp=psipred.w
    matrix=u3.U(s,n,wp,SU3State(2,0),w,r,1,np,1,rp)*A_n(np,n)
  return matrix

# if __name__ == '__main__':
#   print "A"
#   wp=U3State([22,15,10])
#   np=U3State([4,0,0])
#   rp=1
#   w=U3State([21,14,10])
#   n=U3State([2,0,0])
#   r=1
#   s=U3State([20,13,10])
#   psi=SpStateReduced(s,n,r,w)
#   psip=SpStateReduced(s,np,rp,wp)
#   print A(psip,psi,1)
#   wp=U3State([22.5,15.5,10.5])
#   np=U3State([4,0,0])
#   rp=1
#   w=U3State([21.5,14.5,10.5])
#   n=U3State([2,0,0])
#   r=1
#   s=U3State([20.5,13.5,10.5])
#   psi=SpStateReduced(s,n,r,w)
#   psip=SpStateReduced(s,np,rp,wp)
#   print A(psip,psi,1)
#   # # returns 
#     # # 1.12687233964



##################################################################################################
def B(psipred,psired,r0=1):
  """
  Calculates (s np rp wp||a||s n r w). 

  Args:

  Returns: 
    Matrix element (s np rp wp||a||s n r w).  For example
      w=U3State([22,15,10])
      n=U3State([4,0,0])
      r=1
      wp=U3State([21,14,10])
      np=U3State([2,0,0])
      rp=1
      s=U3State3([20,13,10])
      B_weylboson(s,np,rp,wp,n,r,w)
    returns
      -1.28102523044
  """
  matrix=0.0
  if r0==1:
	  w=psired.w
	  wp=psipred.w
	  matrix=utils.parity(wp.pt()-w.pt())*numpy.sqrt(w.dim()/wp.dim())*A(psired,psipred)
  return matrix

# if __name__ == '__main__':
#   print "B"
#   wp=U3State([22,15,10])
#   np=U3State([4,0,0])
#   rp=1
#   w=U3State([21,14,10])
#   n=U3State([2,0,0])
#   r=1
#   s=U3State([20,13,10])
#   psi=SpStateReduced(s,n,r,w)
#   psip=SpStateReduced(s,np,rp,wp)
#   print B(psi,psip,1)

##################################################################################################
def C(psipred,psired,r0=1):
    """ 
    w=U3State([22,15,10])
    n=U3State([4,0,0])
    r=1
    wp=w
    np=n
    rp=r
    s=U3State([20,13,10])
    print C_weylboson(s,np,rp,wp,n,r,w)
    -24.0831891576
    """
    matrix=0.0
    if (psipred.n==psired.n) and (psipred.r==psired.r) and (psipred.w==psired.w) and (r0==1):
      matrix=u3.C_su3(psired.w)
    return matrix
# if __name__ == '__main__':
#   print "C"
#   w=U3State([22,15,10])
#   n=U3State([4,0,0])
#   r=1
#   s=U3State([20,13,10])
#   psi=SpStateReduced(s,n,r,w)
#   print C(psi,psi)
##################################################################################################

def C_n(psipred,psired,r0=1):
  """
  Matrix element of C acting on the raising polynomial space
  """
  matrix=0.0
  n=psired.n
  np=psipred.n
  if np==n:
  	s=psired.s
  	r=psired.r
  	w=psired.w
  	rp=psipred.r
  	wp=psipred.w
  	matrix=u3.U(s,n,wp,SU3State(1,1),w,r,r0,n,1,rp)*u3.C_su3(n)
  return matrix
# if __name__ == '__main__':
#   s=U3State([6,4,4])
#   n2=U3State([4,0,0])
#   r=1
#   n1=U3State([2,2,0])
#   w=U3State([9,5,4])
#   psi=SpStateReduced(s,n2,r,w)
#   # print C_n(psi,psi)

##################################################################################################
####Not sure you need this function...I think it's self conjugate

def C_n_conj(psipred,psired,r0=1):
  """
  Conjugate matrix element of Cn_weylboson
  """
  matrix=0.0
  n=psired.n
  np=psipred.n
  if np==n:
  	s=psired.s
  	r=psired.r
  	w=psired.w
  	rp=psipred.r
  	wp=psipred.w
  	matrix=utils.parity(w.pt()-wp.pt())*numpy.sqrt(w.dim()/wp.dim())*u3.U(s,n,w,SU3State(1,1),wp,rp,1,n,1,r)*u3.C_su3(n)
  return matrix
# if __name__ == '__main__':
#   print "C_n and C_n_conj"
#   w=U3State([22,15,10])
#   n=U3State([4,0,0])
#   r=1
#   wp=U3State([23,14,10])
#   np=n
#   rp=r
#   s=U3State([20,13,10])
#   psi=SpStateReduced(s,n,r,w)
#   psip=SpStateReduced(s,np,rp,wp)
#   print C_n_conj(psi,psi)
#   print C_n(psi, psi)
################################################################################################# 
def C_s1(psi2red,psi1red,r0=1):
  matrix=0.0
  if psi2red.check() and psi1red.check():
    n1=psi1red.n
    n2=psi2red.n
    if n2==n1:
    	s=psi1red.s
    	w1=psi1red.w
    	r1=psi1red.r
    	w2=psi2red.w
    	r2=psi2red.r
    	mult1max=u3.mult(s,n1,w1)
    	mult2max=u3.mult(s,n2,w2)
    	for r1p in range(1,mult1max+1):
    		for r2p in range(1,mult2max+1):
    			matrix=matrix+(u3.phi(s,n1,w1,r1,r1p)
    				*u3.phi(n1,s,w2,r2,r2p)
    				*u3.U(n1,s,w2,SU3State(1,1),w1,r1p,r0,s,1,r2p)
    				*u3.C_su3(s))
  return matrix

# ##################################################################################################
def C_s2(psi2red,psi1red,r0=1):
  matrix=C(psi2red,psi1red,r0)-C_n(psi2red,psi1red,r0)
  return matrix 




##################################################################################################
def C_s3(psipred,psired,r0=1):
  """
  Incorrect derivation reflected in order causing minus sign error.  Need to check with Escher thesis. 
  """
  matrix=0.0
  n=psired.n
  np=psipred.n
  if np==n:
  	s=psired.s
  	r=psired.r
  	w=psired.w
  	rp=psipred.r
  	wp=psipred.w
  	r1max=u3.mult(s,n,wp)
  	for r1 in range(1,r1max+1):
  		matrix=matrix+(u3.phi(n,s,wp,r1,rp)
        *u3.Z(n,s,wp,SU3State(1,1),w,r,r0,s,1,r1)
        *u3.C_su3(s)
        )
  return matrix

# if __name__ == '__main__':
#   print "C_s variations"
#   s=U3State([20,13,10])
#   w=U3State([22,15,10])
#   n=U3State([4,0,0])
#   r=1
#   wp=U3State([23,14,10])
#   np=U3State([4,0,0])
#   rp=1
#   psip=SpStateReduced(s,np,rp,wp)
#   psi=SpStateReduced(s,n,r,w)
#   print "C_s1", C_s1(psip,psi)
#   print "C_s2", C_s2(psip,psi)
#   print "C_s3", C_s3(psip,psi)

# if __name__ == '__main__':
#   print "nlist"
#   s=U3State([6,4,4])
#   w=U3State([8,4,4])
#   nlist1=nlist(s,w)
#   s=U3State([6.5,4.5,4.5])
#   w=U3State([8.5,4.5,4.5])
#   nlist2=nlist(s,w)
#   print nlist1
#   print nlist1==nlist2


#print nlist.__doc__
# if __name__ == '__main__':

#   print "dimensions testing "
#   Nmax=20
#   mmax=0
#   s=U3State([16,13,7])
#   for i in range(Nmax):
#     #print "i", i
#     for j in range(i+1):
#       #print i,j
#       if (s.two()+j)>(s.one()+Nmax-i):
#         continue
#       if (s.three()+(i-j))>(s.two()+j):
#         continue
#       w=U3State([s.one()+Nmax-i,s.two()+j,s.three()+i-j])
#      # print w
#       nset=nlist(s,w)
#       #if nset !={}:
#       m=sum(nset.values())
#       if m>mmax:
#         mmax=m

#   print mmax

# s=U3State(836,SU3State(78,10))
# w=U3State(852,SU3State(77,11))
# print nlist(s,w)

# if __name__ == '__main__':
#   print "poly"
#   Poly=poly(4)
#   for N in Poly:
#       setn=Poly[N]
#       for n in setn:
#           print N,n.cartesian()
