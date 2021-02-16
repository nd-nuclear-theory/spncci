##################################################################################
# Anna E. McCoy
# University of Notre Dame
#
# SPDX-License-Identifier: MIT
##################################################################################


import numpy
import numpy.linalg
from u3states import *
from spstates import *
import u3
import utilities as utils
import u3boson 
import coefficients as coef
import time

#############################################################################
global Kdictionary
Kdictionary={}
global Sdictionary
Sdictionary={}
global SAdictionary
SAdictionary={}
global S2dictionary
S2dictionary={}
global Kdic
Kdic={}
global KAdictionary
KAdictionary={}
global Poly
Poly=u3boson.poly(20)
global Pdic
Pdic={}
#############################################################################
def omega(n,w):
  """
  Calculates Omega(n,w)

  Args: 
    n is the symplectic multiplicity (class U3State)
    w is the symplectic weight (class U3State) 
  
  Returns:
    Omega (float). For example: 
      om=U3State([4.5,.5,.5])
      n=U3State([4,0,0])
      omega(n,om)
     returns
      4.375.
  """
  Omga=0
  nn=n.cartesian()
  ww=w.cartesian()
  #suming over cartesian componenets of nn and ww
  for i in range(3):    
    Omga=Omga+(2*(ww[i]**2)-(nn[i]**2)+8*(ww[i]-nn[i])-2*(i+1)*(2*ww[i]-nn[i]))
    i=i+1
  Omga=Omga/4.0
  return Omga

if __name__ == '__main__':
  w=U3State([8,3,2])
  n=U3State([4,0,0])
  print "omega", w, n, omega(n,w)
  w=U3State([8.5,3.5,2.5])
  n=U3State([4,0,0])
  print "omega", w, n, omega(n,w)
  print 



######################################################################################
def smatrix(s,wp,n1p,r1p,n2p,r2p):
    """
    Calculates the matrix element of S=K^2 and stores the calculated value in the global dictionary Sdictioanry
    
    Args: 
      s-Symplectic irrep label (U3state class)
      w-symplectic weight (U3State class)
      (n1,r1p) and (n2,r2p) are symplectic multiplicities of w in s (U3state class, int). 
    Examples 
      s=U3State([20,13,10])
      w=U3State([22,15,10])
      n1=U3State([2,2,0])
      r1=1
      n2=U3State([4,0,0])
      r2=1
      smatrix(s,w,n1,r1,n1,r1) returns 
        1014.0
      Smatrix(s,w,n1,r1,n2,r2) returns 
        -28.9827534924
    """#results are given in Rowe 1984 jmp 25
    np1={}
    np2={}
    wnp1={}
    wnp2={}
    Sm=0
    #############################################################################################
    #Check if the matrix element of S is already calculated and stored in the global dictioanry Sdictionary
    # if stored return matrix element otherwise calcuate the matrix element and store it.
    #############################################################################################
    if (s,wp,n1p,r1p,n2p,r2p) in Sdictionary: 
        Sm=Sdictionary[s,wp,n1p,r1p,n2p,r2p]
        return Sm
    #############################################################################################
    #Smatrix Calculation
    #############################################################################################     
    # If wp is a bandhead the the matrix element is trivial 1
    if wp==s:
        Sm=1
    #Otherwise, use recursion relationship to calculate S
    else: 
        psi1pred=SpStateReduced(s,n1p,r1p,wp)
        psi2pred=SpStateReduced(s,n2p,r2p,wp)
        ############################################################################################
        #Generating sets of n1 and n2 which will have nonzero u3boson.A_n matrix elements 
        ############################################################################################
        #set of allowed n1
        n1set=[]
        [m1,m2,m3]=n1p.cartesian()
        if m1-2>=m2:
            n1set.append(U3State([m1-2,m2,m3]))
        if m2-2>=m3:
            n1set.append(U3State([m1,m2-2,m3]))
        if m3-2>=0:
            n1set.append(U3State([m1,m2,m3-2]))
        #generating list of allowed n2's
        n2set=[]
        [m1,m2,m3]=n2p.cartesian()
        if m1-2>=m2:
            n2set.append(U3State([m1-2,m2,m3]))
        if m2-2>=m3:
            n2set.append(U3State([m1,m2-2,m3]))
        if m3-2>=0:
            n2set.append(U3State([m1,m2,m3-2]))
        ############################################################################################
        #summing over w,n1,r1,n2,r2 to obtain the matrix element of S
        ############################################################################################
        Nw=wp.N-2
        #compute list of w labels for states that may have nonzero matrix elements of u3boson.A_n
        lmset=u3.couple(wp,SU3State(0,2))
        #summing over w
        for lm in lmset:
          #converging SU3State lm to U3State w
          w=U3State(Nw,lm)
          #summing over possible n1
          for n1 in n1set:
            #checking for each s, n1, and w if (s, n1,r, w) are a state in the symplectic irrep s
            r1max=u3.mult(s,n1,w)
            if r1max==0:
              continue
            #summing over possible 
            for n2 in n2set:
              #checking for each s, n2, and w if (s, n1,r, w) are a state in the symplectic irrep s
              r2max=u3.mult(s,n2,w)
              if r2max==0:
                continue
              #suming over r1
              for r1 in range(1,r1max+1):
                psi1red=SpStateReduced(s,n1,r1,w)
                #summing over r2
                for r2 in range(1,r2max+1):
                  psi2red=SpStateReduced(s,n2,r2,w)
                  #calcualting the matrix element of S via a recursion relation
                  Sm=Sm+(omega(n1p,wp)-omega(n1,w))*u3boson.A(psi1pred,psi1red)*u3boson.A(psi2pred,psi2red)*smatrix(s,w,n1,r1,n2,r2)#this will be the calculation
       ## 
        Sm=2*Sm/(wp.N-s.N)
        #Storing the calculated matrix element of S in a global dictionary
        Sdictionary[(s,wp,n1p,r1p,n2p,r2p)]=Sm
    return Sm


#a few examples
if __name__ == '__main__':
  print "smatrix"
  s=U3State([20,13,10])
  w=U3State([22,15,10])
  n1=U3State([2,2,0])
  r1=1
  n2=U3State([4,0,0])
  r2=1
  print "diagonal    ", smatrix(s,w,n1,r1,n1,r1)
  # #         #1014.0
  print"off-diagonal", smatrix(s,w,n1,r1,n2,r2)
  # #         #-28.9827534924      
  print 
  s=U3State([20.5,13.5,10.5])
  w=U3State([22.5,15.5,10.5])
  n1=U3State([2,2,0])
  r1=1
  n2=U3State([4,0,0])
  r2=1
  print "diagonal    ", smatrix(s,w,n1,r1,n1,r1)
  # #         #1014.0
  print"off-diagonal", smatrix(s,w,n1,r1,n2,r2)
  # #         #-28.9827534924      
  print 
######################################################################################
def smatrix_asymptotic(s,wp,n1p,r1p,n2p,r2p):
    """
    Calculates the matrix element of S=K^2 using the asymptotic approxiamtion .

    Args: 
      s is the symplectic irrep label (U3state class)
      w is the symplectic weight (U3state class)
      (n1,r1p) and (n2,r2p) are symplectic multiplicities of w in s. (U3State class,int)

    Examples 
      s=U3State([20,13,10])
      w=U3State([22,15,10])
      n1=U3State([2,2,0])
      r1=1
      n2=U3State([4,0,0])
      r2=1
      smatrix_asymptotic(s,w,n1,r1,n1,r1) returns
        1014.0
    """
    #############################################################################################
    # Check if the matrix element of S is already calculated and stored in the global dictioanry Sdictionary
    # if stored return matrix element otherwise calcuate the matrix element and store it.
    #############################################################################################
    if (s,wp,n1p,r1p,n2p,r2p) in SAdictionary: 
        SmA=SAdictionary[s,wp,n1p,r1p,n2p,r2p]
    #############################################################################################
    #Smatrix Calculation
    #############################################################################################     
    else:
      np1={}
      np2={}
      wnp1={}
      wnp2={}
      SmA=0

      #In asymptotic limit, off diagonal matrix elements are set to zero. 
      if (n1p!=n2p) or (r1p!=r2p):
        SmA=0.0
      else:
        # If wp is a bandhead the the matrix element is trivial 1
        if wp==s:
            SmA=1
        #Otherwise, use recursion relationship to calculate S
        else: 
            ##############################################################################################
            #Smatrix Calculation
            ##############################################################################################  
            psipred=SpStateReduced(s,n1p,r1p,wp)
            #Generating sets of n which will have nonzero u3boson.A_n matrix elements
            nset=[]
            nn=U3State([n1p.one()-2,n1p.two(),n1p.three()])
            if nn.check():
                nset.append(nn)
            nn=U3State([n1p.one(),n1p.two()-2,n1p.three()])
            if nn.check():
                nset.append(nn)
            nn=U3State([n1p.one(),n1p.two(),n1p.three()-2])
            if nn.check():
                nset.append(nn)
           
            #summing over w,n,r
            ##############################################
            Nw=wp.N-2
            lmset=u3.couple(wp,SU3State(0,2))
            #generating list of w's and summing over them
            for lm in lmset:
              w=U3State(Nw,lm)
              #summing over n
              for n in nset:
                #checking state (s,n,r,w) is allowed for  some r
                rmax=u3.mult(s,n,w)
                if rmax==0:
                  continue
                for r in range(1,rmax+1):
                    psired=SpStateReduced(s,n,r,w)
                    #Calculating the asymptotic S matrix 
                    SmA=SmA+(
                    (omega(n1p,wp)-omega(n,w))
                    *u3boson.A(psipred,psired)
                    *u3boson.A(psipred,psired)
                    *smatrix_asymptotic(s,w,n,r,n,r)
                    )
           ## 
            SmA=2*SmA/(n1p.N)
      #caching S matrix element in SAdictionary
      SAdictionary[s,wp,n1p,r1p,n2p,r2p]=SmA
    return SmA
if __name__ == '__main__':
  print "smatrix_asymptotic"
  s=U3State([20,13,10])
  w=U3State([22,15,10])
  n1=U3State([2,2,0])
  r1=1
  n2=U3State([4,0,0])
  r2=1
  print "diagonal    ", smatrix_asymptotic(s,w,n1,r1,n1,r1)
  # #         #1014.0
  print "off-diagonal", smatrix_asymptotic(s,w,n1,r1,n2,r2)
  # #         #-28.9827534924      
  print
######################################################################################
#Calculates the K matrix
def kmatrix(s,om):
  """
  K(s,om) calculates the K matrix for irrep "s" and weight state "om"
  Input "s" and "om" are type class U3
  Calls functions "nlist", "Smatrix"
  Output is matrix with entry key nitr which has the form {(ni,ri)} where ni variables type class U3 and ri is multiplicity of om in sxni su3 product, the key indexes which entry in the matrix we are looking at  
  Example: 
    s=U3State([20,13,10])
    w=U3State([22,15,10])
    print Kmatrix(s,w)
      matrix([[ 31.83998951,  -0.46375447],[ -0.46375447,  30.65591186]]), {(<u3.U3 instance at 0x101034170>, 1): 1, (<u3.U3 instance at 0x101034950>, 1): 0}
  """
  #############################################################################################
    # Check if the matrix element of S is already calculated and stored in the global dictioanry Sdictionary
    # if stored return matrix element otherwise calcuate the matrix element and store it.
  #############################################################################################
     
  #Check if value in Kdictionary
  if (s,om) in Kdictionary:
    return Kdictionary[s,om]
  #if not in Kdictionary, calculate Kmatrix
  else:
    nitr={}
    n1list=u3boson.nlist(s,om)
    nlen=sum(n1list.values())
    if nlen==0:
      K=numpy.matrix(0)
    else:
      #Ceate a empy matrix to store the Smatrix values in
      Smatrix=numpy.identity(nlen)
      i=0    
      for n1 in n1list:
        r1max=n1list[n1]    
        for r1 in range(1,r1max+1):
            #nitr is a list that indexes the row (collumn) that corresponds to the multiplicity pair (n,rho) starting with i=0
            nitr[n1,r1]=i
            j=0
            if j>i:
              continue
            for n2 in n1list:
                r2max=n1list[n2]
                for r2 in range(1,r2max+1):
                    #Inputing Smatrix values into matrix
                    sm=smatrix(s,om,n1,r1,n2,r2)
                    Smatrix[i,j]=sm
                    Smatrix[j,i]=sm
                    j=j+1
            i=i+1

    #diagonalizing S and solving for the square root of the matrix.  
      #Turn smatrix into a matrix
      Smatrix=numpy.matrix(Smatrix)
      #Get eigenvalues and eigenvectors S=Uk^2U^\dagger
      ks,Umatrix=numpy.linalg.eig(Smatrix)
      for i in range(len(ks)):
        if abs(ks[i])<10**(-7):
          ks[i]=0.0 
      k=numpy.sqrt(ks)
      #Generate diagonal matrix with k along diagonal
      km=numpy.identity(len(k))*k
      #KK^dagger=S=Uk^2U^\dagger=(UkU^\dagger)(UkU^\dagger)^\dagger
      K=Umatrix*km*Umatrix.H
      #storing calculated K matrix with corresponding n's in global Kdictionary as value under key (s,Om)    
      Kdictionary[s,om]=(K,nitr)
    return K,nitr


######################################################################################
#Calculates the K matrix
def kmatrix_asymptotic(s,om):
  """
  K(s,om) calculates the K matrix for irrep "s" and weight state "om"
  Input "s" and "om" are type class U3
  Calls functions "nlist", "Smatrix"
  Output is matrix with entry key nitr which has the form {(ni,ri)} where ni variables type class U3 and ri is multiplicity of om in sxni su3 product, the key indexes which entry in the matrix we are looking at  
  Example: 
    s=U3State([20,13,10])
    w=U3State([22,15,10])
    print Kmatrix(s,w)
      matrix([[ 31.83998951,  -0.46375447],[ -0.46375447,  30.65591186]]), 
  """
  #Check if value has already been calculated and stored in Kdictionary and calculates if it isn't
  if (s,om) in KAdictionary:
    K,nitr=KAdictionary[s,om]
  else:
    nitr={}
    nList=u3boson.nlist(s,om)
    nlen=sum(nList.values())
    if nlen==0:
      K=numpy.matrix(0)
    else:
      #Ceate a empy matrix to store the Smatrix values in
      Smatrix=numpy.identity(nlen)
      i=0    
      for n in nList:
        rmax=nList[n]    
        for r in range(1,rmax+1):
            #nitr is a list that indexes the row (collumn) that corresponds to the multiplicity pair (n,rho) starting with i=0
            nitr[n,r]=i
              #Inputing Smatrix values into matrix
            Smatrix[i,i]=smatrix_asymptotic(s,om,n,r,n,r)  
            i=i+1

    #diagonalizing S and solving for the square root of the matrix.  
      #Turn smatrix into a matrix
      Smatrix=numpy.matrix(Smatrix)
      #Get eigenvalues and eigenvectors S=Uk^2U^\dagger
      ks,Umatrix=numpy.linalg.eig(Smatrix)
      for i in range(len(ks)):
        if abs(ks[i])<10**(-7):
          ks[i]=0.0 
      k=numpy.sqrt(ks)
      #Generate diagonal matrix with k along diagonal
      km=numpy.identity(len(k))*k
      #KK^dagger=S=Uk^2U^\dagger=(UkU^\dagger)(UkU^\dagger)^\dagger
      K=Umatrix*km*Umatrix.H
      #storing calculated K matrix with corresponding n's in global Kdictionary as value under key (s,Om)    
      KAdictionary[s,om]=(K,nitr)
  return K,nitr

if __name__ == '__main__':
  print "kmatrices and index"
  s=U3State([20,13,10])
  w=U3State([22,15,10])
  t1=time.time()
  K,ns=kmatrix(s,w)
  t2=time.time()
  t3=time.time()
  KA,nsA=kmatrix_asymptotic(s,w)
  t4=time.time()
  for (n,r) in ns:
    print n, r, ns[n,r]

  numpy.set_printoptions(precision=3,formatter={"all": lambda x:format(x,"8.4f")})   
  print "kmatrix           ", t2-t1
  print K
  print 
  print "kmatrix_asymptotic", t4-t3
  print KA
  print 
    
  s=U3State([20.5,13.5,10.5])
  w=U3State([22.5,15.5,10.5])
  K,ns=kmatrix(s,w)
  print "half int K"
  print K


########################################################################

def A(psi2,psi1,r0=1):
  """
  Exact matrix elements for Adagger

  Args:
    psi1,psi1 are reduced basis states (SpStateReduced class)
    r0 is outer multiplicity of coupling psi1xA->psi2

  Example:
    s=U3State([20,13,10])
    wp=U3State([22,15,10])
    np=U3State([4,0,0])
    rp=1
    w=U3State([21,14,10])
    n=U3State([2,0,0])
    r=1
    psip=SpStateReduced(s,np,rp,wp)
    psi=SpStateReduced(s,n,r,w)

    A(psip,psi) returns 
      6.272527136548157
  """
#Calculate k and U for w1 and s, then k and U for s and w2.
#sum over all n states for w1 but only need to sum over 3 possible n states for w2.  
  matrix=0.0
  if r0==1:
    s=psi2.s
    n2=psi2.n
    p2=psi2.r
    w2=psi2.w
    n1=psi1.n
    p1=psi1.r
    w1=psi1.w
    K2,nitr2=kmatrix(s,w2)
    K1,nitr1=kmatrix(s,w1)
    if K1.all()==0 or K2.all()==0:
      return matrix
    K1inv=K1**(-1)
    #position of (n2,p2) in K2 matrix
    j2=nitr2[n2,p2] 
    #position of (n1,p1) in K1 matrix      
    j1=nitr1[n1,p1]       
    #sum over all (n,p) multiplicities of w1  ######problem
    for m1, r1 in nitr1: 
        i1=nitr1[m1,r1]
        psim1=SpStateReduced(s,m1,r1,w1)
        #first possibility for n2
        m2set=[U3State([m1.one()+2,m1.two(),m1.three()]),U3State([m1.one(),m1.two()+2,m1.three()]),U3State([m1.one(),m1.two(),m1.three()+2])]
        for m2 in m2set:
          if m2.check()==True:
              r2max=u3.mult(s,m2,w2) 
              # this is summing over su3 multiplicities of w2 multiplicities 
              for r2 in range(1, r2max+1):   
                  # index for picking element of K matrix for w2                              
                  i2=nitr2[m2,r2]     
                  psim2=SpStateReduced(s,m2,r2,w2)                                     
                  matrix=matrix+K2[j2,i2]*u3boson.A(psim2,psim1)*K1inv[j1,i1]
  return matrix

###################################################################################################################################
def A_asymptotic(psi2red,psi1red,r0=1):
  #a,b2,p2,c2,b1,p1,c1,r0=1):
  """
  Assymptotic approximation for the RME of A

  args:
    psi2red and psi1red are the basis states (class SpStateReduced)
    r0 is the outer multiplicity of psi1redxA->psi2red (int=1)

  returns: 
    matrix-the RME of A (float) 
  """
  matrix=0.0
  #check that outer multipicity is 1 else, RME is zero
  if r0==1:
    #check that the basis states are good Sp reduced states 
    if psi1red.check() and psi2red.check():
      n2=psi2red.n
      w2=psi2red.w
      n1=psi1red.n
      w1=psi1red.w
      #calculate the omega factor
      om=omega(n2,w2)-omega(n1,w1)
      #In some cases om<0 and so matrix elements would not be real.  
      #But in these cases, matrix elements are approximately zero, so we jut set to zero. 
      if om>0:
        matrix=numpy.sqrt(om)* u3boson.A(psi2red,psi1red)
      #print numpy.sqrt(om)
  return matrix
 
###################################################################################################################################
def A_u3boson(psi2red,psi1red,r0=1):
  """
  Calculates RME of A operator calculated using the u3-boson approximation 

  Args:
    psi2red, psi1red are the basis states (class SpStateReduced)
    r0 is the outer multiplicity of the couping of psi1red with A to psi2red (int).  Should be 1

  Returns:
    matrix which is the reduced matrix element of A (float)
  """
  matrix=0.0
  #check that outer multiplicity is 1 and that psi1red and psi2red are good reduced symplectic states
  if r0==1:# and psi1red.check() and psi2red.check():
    s=psi1red.s
    #check that both states belong to the same irrep
    if s==psi2red.s:
      #calculate the matrix element
      matrix=numpy.sqrt(2.*s.N/3)*u3boson.A(psi2red,psi1red)
  return matrix  

if __name__ == '__main__':
    #setting up the basis states
    s=U3State([20,13,10])
    wp=U3State([22,15,10])
    np=U3State([4,0,0])
    rp=1
    w=U3State([21,14,10])
    n=U3State([2,0,0])
    r=1
    psip=SpStateReduced(s,np,rp,wp)
    psi=SpStateReduced(s,n,r,w)
    print "Comparison of exact RME of A with asymptotic and U(3)-boson approx."
    print "A           ", A(psip,psi)
    print "A_asymptotic", A_asymptotic(psip,psi)
    print "A_u3boson   ", A_u3boson(psip,psi)
    print 
    #checking that it still works for projective symplectic representations 
    s=U3State([20.5,13.5,10.5])
    wp=U3State([22.5,15.5,10.5])
    np=U3State([4,0,0])
    rp=1
    w=U3State([21.5,14.5,10.5])
    n=U3State([2,0,0])
    r=1
    psip=SpStateReduced(s,np,rp,wp)
    psi=SpStateReduced(s,n,r,w)
    print "half int A  ", A(psip,psi)

####################################################################################################################
def P(n0,psipred,psired,r0):
  """
  Calculates the matrix element of the polynomial of A's with u3 character n0

  args:
    n0 is U(3) symmetry of polynomial (U3State class)
    psipred,psired are reduced basis states (SpStateReduced class)
    r0 is outer multiplicity of psired.wxA->psipred.w

  returns:
    prme which is the reduced matrix element of P (float)
  """
  # if (n0,psipred,psired,r0) in Pdic:
  #   prme=Pdic[n0,psipred,psired,r0]
  # else:
  prme=0
  ## Calculate the K matrix element and generating index list
  K,nitr=kmatrix(psired.s,psired.w)
  ## inverting the K matrix 
  Kinv=K**(-1)
  ## calculating the Kp matrix and generating index list 
  Kp, npitr=kmatrix(psipred.s,psipred.w)
  ## obtaining index of psired.n,psired.r in K matrix
  i=nitr[psired.n,psired.r]
  ## obtaining index of psipred.n,psipred.r in Kp matrix
  ip=npitr[psipred.n, psipred.r]
  # suming over index values of n1,r1 in K matrix
  for n1,r1, in nitr:
    ## obtaining index of n1,r1 in K matrix
    j=nitr[n1,r1]
    # suming over index values of n1p,r1p in Kp matrix
    for n1p, r0p in npitr:
      ## obtaining index of n1p,r1p in Kp matrix
      jp=npitr[n1p,r0p]
      ## Calculating the max multiplicity of sxn1p->wp 
      r1pmax=u3.mult(n0,n1,n1p)
      ## summing over multiplicity
      for r1p in range(1,r1pmax+1):
        ##calculating RME of P
        prme=prme+(
          Kinv[i,j]
          *u3.U(n0,n1,psipred.w,psired.s,n1p,r1p,r0p,psired.w,r1,r0)
          *coef.bcoef(n0,n1,n1p,r1p)
          *Kp[ip,jp])
    #Pdic[n0,psipred,psired,r0]=prme
  return prme


if __name__ == '__main__':
  s=U3State([20,13,10])
  w=U3State([24,13,10])
  n=U3State([4,0,0])
  r=1
  wp=U3State([22,13,10])
  np=U3State([2,0,0])
  rp=1
  n0=U3State([2,0,0])
  psip=SpStateReduced(s,np,rp,wp)
  psi=SpStateReduced(s,n,r,w)
  print "P            ", P(n0,psi,psip,1)

  # n0=U3State(12, SU3State(2,2))
  # psip=SpStateReduced(U3State(16,SU3State(2,1)),U3State(18,SU3State(0,0)),1,U3State(34,SU3State(2,1)))
  # psi=SpStateReduced(U3State(16,SU3State(2,1)),U3State(6,SU3State(2,2)),1,U3State(22,SU3State(3,2)))
  #print u3.mult(psi.w,n0,psip.w)
  #print P(n0,psip,psi,2)
####################################################################################################################
def P_asymptotic(n0,psipred,psired,r0):
  """
  Calculates the asymptotic matrix element of the polynomial of A's with u3 character n0

  args:
    n0 is U(3) symmetry of polynomial (U3State class)
    psipred,psired are reduced basis states (SpStateReduced class)
    r0 is outer multiplicity of psired.wxA->psipred.w

  returns:
    prme which is the reduced matrix element of P (float)
  """
  prme=0
  ## Calculate the K matrix element and generating index list
  K,nitr=kmatrix_asymptotic(psired.s,psired.w)
  ## inverting the K matrix 
  Kinv=K**(-1)
  ## calculating the Kp matrix and generating index list 
  Kp, npitr=kmatrix_asymptotic(psipred.s,psipred.w)
  ## obtaining index of psired.n,psired.r in K matrix
  i=nitr[psired.n,psired.r]
  ## obtaining index of psipred.n,psipred.r in Kp matrix
  ip=npitr[psipred.n, psipred.r]
  # suming over index values of n1,r1 in K matrix
  for n1,r1, in nitr:
    ## obtaining index of n1,r1 in K matrix
    j=nitr[n1,r1]
    # suming over index values of n1p,r1p in Kp matrix
    for n1p, r0p in npitr:
      ## obtaining index of n1p,r1p in Kp matrix
      jp=npitr[n1p,r0p]
      ## Calculating the max multiplicity of sxn1p->wp 
      r1pmax=u3.mult(n0,n1,n1p)
      ## summing over multiplicity
      for r1p in range(1,r1pmax+1):
        ##calculating RME of P
        prme=prme+(
          Kinv[i,j]
          *u3.U(n0,n1,psipred.w,psired.s,n1p,r1p,r0p,psired.w,r1,r0)
          *coef.bcoef(n0,n1,n1p,r1p)
          *Kp[ip,jp])
  return prme



#######################################################################################################
def B(psi2,psi1,r0=1):
  """
  Exact RME B in the symplectic basis calculated using VCS theory
  
  args:
    psi2,psi1 are reduced symplectic basis states (class SpStateReduced)
    r0 is outer mutplicity of psi1xB->psi2
 
  Example:
    s=U3([20,13,10])
    wp=U3([22,15,10])
    np=U3([4,0,0])
    rp=1
    w=U3([21,14,10])
    n=U3([2,0,0])
    r=1
    psip=SpStateReduced(s,np,rp,wp)
    psi=SpStateReduced(s,n,r,w)

    B(psip,psi,1) returns 
      -7.13059078469 
  """
#Calculate k and U for w1 and s, then k and U for s and w2.
#sum over all n states for w1 but only need to sum over 3 possible n states for w2.
  matrix=0.0
  if psi1.check() and psi2.check():
    s=psi2.s
    n2=psi2.n
    p2=psi2.r
    w2=psi2.w
    n1=psi1.n
    p1=psi1.r
    w1=psi1.w  
    K2,nitr2=kmatrix(s,w2)
    K1,nitr1=kmatrix(s,w1)
    if K1.all()==0 or K2.all()==0:
      return matrix
    K2=K2**(-1)
    K1inv=K1
    #position of (n2,p2) in K2 matrix
    i2=nitr2[n2,p2]  
    #position of (n1,p1) in K1 matrix     
    j1=nitr1[n1,p1]       
    #sum is restricted by the coupling of n2x(2,0)->n1
    for m2,r2 in nitr2:
      j2=nitr2[m2,r2] 
      psim2=SpStateReduced(s,m2,r2,w2)
      #first possibility for n2 given that n1X(2,0)->n2
      m1set=[U3State([m2.one()+2,m2.two(),m2.three()]),U3State([m2.one(),m2.two()+2,m2.three()]), U3State([m2.one(),m2.two(),m2.three()+2])]
      for m1 in m1set: 
        if m1.check()==True:
          #multiplicity of the sxm1 coupling
          r1max=u3.mult(s,m1,w1)
          #summing over m1,r1  
          for r1 in range(1, r1max+1): 
              # index for picking element of K matrix for w2                                
              i1=nitr1[m1,r1]   
              psim1=SpStateReduced(s,m1,r1,w1)                                       
              matrix=matrix+K2[i2,j2]*u3boson.B(psim2,psim1)*K1inv[j1,i1]
        #second possibility for n2       
  return matrix

##################################################
def B2(psi2red,psi1red,r0=1):
  """
  Alternative method for calculating exact RME of B in the symplectic basis using exact RME of A in symplectic basis
  
  args:
    psi2red, psi1red are reduced symplectic basis states (SpStateReduced class)
    r0 is outer multiplicity of psi1redxB->psi2red (int=1)
 
  Example:
    s=U3([20,13,10])
    wp=U3([22,15,10])
    np=U3([4,0,0])
    rp=1
    w=U3([21,14,10])
    n=U3([2,0,0])
    r=1
    psipred=SpStateReduced(s,np,rp,wp)
    psired=SpStateReduced(s,n,r,w)
    
    B(psipred,psired) returns
      -5.51771905958
  """
  matrix=0.0
  #check that outer mutliplicity is 1, else RME=0
  if r0==1:
    w1=psi1red.w
    w2=psi2red.w
    #calculate RME
    matrix=utils.parity(w1.pt()+w2.pt())*numpy.sqrt(w1.dim()/w2.dim())*A(psi1red,psi2red,r0)
  return matrix
  
##########################################################
def B_asymptotic(psi2red,psi1red,r0=1):
  """
  Calculates the asymptotic approximation of the RME of B in the symplectic basis
  
  args:
    psi2red,psi1red are reduced symplectic basis states (SpStateReduced class)
    r0 is outer multiplicity of psi1redxB->psi2red

  Example:
    s=U3State([20,13,10])
    wp=U3State([22,15,10])
    np=U3State([4,0,0])
    rp=1
    w=U3State([21,14,10])
    n=U3State([2,0,0])
    r=1
    psipred=SpStateReduced(s,np,rp,wp)
    psired=SpStateReduced(s,n,r,w)
    B(psired, psipred) returns 
      -7.13059078469  
    """
  matrix=0.0
  #check that outer multiplicity is 1, else RME=0
  if r0==1:
    w1=psi1red.w
    w2=psi2red.w
    #calculate the RME
    matrix=utils.parity(w1.pt()+w2.pt())*numpy.sqrt(w1.dim()/w2.dim())*A_asymptotic(psi1red,psi2red)
  return matrix

###################################################################################################################################
def B_u3boson(psi2red,psi1red,r0=1):
  """
  Calculates the U(3)-boson approximation of B in the symplectic basis

  args:
    psi2red,psi1red are reduced symplectic basis state (SpStateReduced class)
    r0 is the outer multiplicity of psi1redxB->psi2red (int=1)

  returns:
    matrix (float)
  """
  matrix=0.0
  #check that outer multiplicity is 1
  if r0==1:
    #check that psi1red and psi2red are good reduced symplectic states 
    if psi1red.check() and psi2red.check():
      s=psi1red.s
      #check that basis states belong to same symplectic irrep
      if s==psi2red.s:
        #calculate the RME
        matrix=numpy.sqrt(2.*s.N/3)*u3boson.B(psi2red,psi1red)
      #print numpy.sqrt(om)
  return matrix  


if __name__ == '__main__':
  s=U3State([20,13,10])
  w=U3State([22,15,10])
  n=U3State([4,0,0])
  r=1
  wp=U3State([21,14,10])
  np=U3State([2,0,0])
  rp=1
  psip=SpStateReduced(s,np,rp,wp)
  psi=SpStateReduced(s,n,r,w)
  print "B           ", B(psip,psi)
  print "B2          ", B2(psip,psi)
  print "B_asymptotic", B_asymptotic(psip,psi)
  print "B_u3boson   ", B_u3boson(psip,psi)
  print 
####################################################################################################################
def C(psi2red,psi1red,r0=1):
  """
  Calculates the exact matrix elements of C in the symplectic basis 
  
  args:
    psi2red,psi1red are reduced symplectic basis state (SpStateReduced class)
    r0 is the outer multiplicity of psi1redxC->psi2red (int=1)

  Example:
    s=U3State([20,13,10])
    wp=U3State([22,15,10])
    np=U3State([4,0,0])
    rp=1
    w=U3State([21,14,10])
    n=U3State([2,0,0])
    r=1
    psipred=SpStateReduced(s,np,rp,wp)
    psired=SpStateReduced(s,n,r,w)

    C(psipred,psired,1) returns 
      0

    C(psired,psired) returns
      -22.4499443206
  """
  matrix=0.0
  #check that outer mutliplicity is 1, else RME=0
  if r0==1:
    #check that psi2red and psi1red are good reduced symplectic states 
    if psi2red.check() and psi1red.check():
      #calculate matrix element 
      matrix=u3boson.C(psi2red,psi1red,r0)
      # set nonzer0 matrix elements with abs value less than 10^(-13) to zero
      if abs(matrix)<10**(-13):
        matrix=0.0
  return matrix


if __name__ == '__main__':
  s=U3State([20,13,10])
  wp=U3State([22,15,10])
  np=U3State([4,0,0])
  rp=1
  w=U3State([21,14,10])
  n=U3State([2,0,0])
  r=1
  print "C", C(psip,psip)
  print 

###############################################################################################
def C_n(psi2red,psi1red,r0=1):
    """
    Calculate the reduced matrix element of the n-space C operator in the symplectic basis

    args:
      psi2red,psi1red are reduced symplectic basis state (SpStateReduced class)
      r0 is the outer multiplicity of psi1redxB->psi2red (int=1)
    """
    matrix=0.0
    #check that psi2red and psi1red are good reduced symplectic states 
    if psi2red.check() and psi1red.check():
      s=psi1red.s
      n2=psi2red.n
      p2=psi2red.r
      w2=psi2red.w
      n1=psi1red.n
      p1=psi1red.r
      w1=psi1red.w
      #compute K matrices for psi2red and psi1red 
      K2,n2itr=kmatrix(s,w2)
      K1,n1itr=kmatrix(s,w1)
      #invert K matrix of psi1red
      K1inv=K1**(-1)
      #obtain indexes of desired K-matrix elements 
      i2=n2itr[n2,p2]
      j1=n1itr[n1,p1]
      #sum over symplectic multipliciites 
      for nb,rb in n1itr:
        i1=n1itr[nb,rb]
        psi1bred=SpStateReduced(s,nb,rb,w1)
        rbpmax=u3.mult(s,nb,w2)
        for rbp in range(1,rbpmax+1):
          j2=n2itr[nb,rbp]
          psi2bred=SpStateReduced(s,nb,rbp,w2)
          matrix=matrix+K2[i2,j2]*u3boson.C_n(psi2bred,psi1bred,r0)*K1inv[i1,j1]
      if abs(matrix)<(10**(-13)):
        matrix=0.0
    return matrix

#################################################################################################################
def C_n_asymptotic(psi2red,psi1red,r0=1):
    """
    Calculated reduced matrix element of the n-space C operator using the asymptotic approximation in the symplectic basis

    args:
      psi2red,psi1red are reduced symplectic basis state (SpStateReduced class)
      r0 is the outer multiplicity of psi1redxB->psi2red (int)

    """
    matrix=0.0
    #check that psi2red and psi1red are good reduced symplectic states 
    if psi2red.check() and psi1red.check():
      s=psi1red.s
      n2=psi2red.n
      p2=psi2red.r
      w2=psi2red.w
      n1=psi1red.n
      p1=psi1red.r
      w1=psi1red.w
      #Calculate the K matrices for psi1red and psi2red
      K1,n1itr=kmatrix_asymptotic(s,w1)
      K2,n2itr=kmatrix_asymptotic(s,w2)
      #inver the psi1red K matrix
      K1inv=K1**(-1)
      #look up indices of desired K-matrix elements 
      i2=n2itr[n2,p2]
      j1=n1itr[n1,p1]
      #sum over symplectic multiplicities 
      for nb,rb in n1itr:
        i1=n1itr[nb,rb]
        psi1bred=SpStateReduced(s,nb,rb,w1)
        rbpmax=u3.mult(s,nb,w2)
        for rbp in range(1,rbpmax+1):
          j2=n2itr[nb,rbp]
          psi2bred=SpStateReduced(s,nb,rbp,w2)
          matrix=matrix+K2[i2,j2]*u3boson.C_n(psi2bred,psi1bred,r0)*K1inv[i1,j1]
      # setting matrix elements whose absolute value is less than 10^(-13) to zero
      if abs(matrix)<(10**(-13)):
        matrix=0.0
    return matrix  


if __name__ == '__main__':
    
  s=U3State([20,13,10])
  w=U3State([22,15,10])
  n=U3State([4,0,0])
  r=1
  wp=U3State([23,14,10])
  np=U3State([4,0,0])
  rp=1
  psi1red=SpStateReduced(s,np,rp,wp)
  psi2red=SpStateReduced(s,n,r,w)
  print "C_n           ", C_n(psi2red,psi1red,r0=1)
  print "C_n_asymptotic", C_n_asymptotic(psi2red,psi1red,r0=1)
  print 
#############################################################################################
# def C_n_conj(psi2red,psi1red,r0=1):
#     """
#     su3 conjuage matrix element of C_n
#     """
#     s=psi1red.s
#     w1=psi1red.w
#     n1=psi1red.n
#     p1=psi1red.r
#     w2=psi2red.w
#     n2=psi2red.n
#     p2=psi2red.r
#     matrix=0.0
#     K2,n2itr=kmatrix(s,w2)
#     K1,n1itr=kmatrix(s,w1)
#     #note that inverted matrix here is opposite of Cn, here K2 is inverted instead of K1
#     K2inv=K2**(-1)
#     i2=n2itr[n2,p2]
#     j1=n1itr[n1,p1]
#     for nb,rb in n1itr:
#       i1=n1itr[nb,rb]
#       psi1bred=SpStateReduced(s,nb,rb,w1)
#       rbpmax=u3.mult(s,nb,w2)
#       for rbp in range(1,rbpmax+1):
#         j2=n2itr[nb,rbp]
#         psi2bred=SpStateReduced(s,nb,rbp,w2)
#         matrix=matrix+K2inv[i2,j2]*u3boson.C_n_conj(psi2bred,psi1bred)*K1[i1,j1]
#     return matrix
# if __name__ == '__main__':
#   s=U3State([20,13,10])
#   wp=U3State([22,13,10])
#   np=U3State([2,0,0])
#   rp=1
#   w=U3State([21,14,10])
#   n=U3State([2,0,0])
#   r=1
#   psip=SpStateReduced(s,n,r,w)
#   print "C_n     ", C_n(psip,psip)#*numpy.sqrt(dim(wp)/dim(w))*parity(wp.lbda()+wp.mu()-w.lbda()-w.mu())
#   print "C_n_conj", C_n_conj(psip,psip)
#   print 
#   #print C_n_conj(s,n,r,w,np,rp,wp)
#   # #print U_su3(7,3,2,0,9,3,1,1,7,4,1,1,2,0,1,1)
#   # #print numpy.sqrt(Dim(7,4)/Dim(9,3))*U_su3(7,3,2,0,7,4,1,1,9,3,1,1,2,0,1,1)


############################################################################################
def C_s2(psi2red,psi1red,r0=1):
  """
  Reduced matrix element of SU(3) generators C^(1,1) restrictued to the lowest grade SU(3) irrep s calculated by taking the difference of
  the unrestricted generator and the generator restricted to the raising A_asymptoticnomial space n 

  args:
    psi2red,psi1red are reduced symplectic basis state (SpStateReduced class)
    r0 is the outer multiplicity of psi1redxB->psi2red (int)
  """
  matrix=0.0
  n2=psi2red.n
  n1=psi1red.n
  matrix=C(psi2red,psi1red,r0)-C_n(psi2red,psi1red,r0)
  if abs(matrix)<10**(-13):
    matrix=0.0
  return matrix



def C_s2_asymptotic(psi2red,psi1red,r0=1):
  """
  Asymptotic reduced matrix element of SU(3) generators C^(1,1) restrictued to the lowest grade SU(3) irrep s calculated by taking the difference of
  the unrestricted generator and the generator restricted to the raising A_asymptoticnomial space n 

  args:
    psi2red,psi1red are reduced symplectic basis state (SpStateReduced class)
    r0 is the outer multiplicity of psi1redxB->psi2red (int)
  """
  matrix=0.0
  n2=psi2red.n
  n1=psi1red.n
  matrix=C(psi2red,psi1red,r0)-C_n_asymptotic(psi2red,psi1red,r0)
  if abs(matrix)<10**(-13):
    matrix=0.0
  return matrix  

def C_s1(psi2red,psi1red,r0=1):
    """
    Reduced matrix elements of restricted C operator caculated via VCS theory 

    args:
      psi2red,psi1red are reduced symplectic basis state (SpStateReduced class)
      r0 is the outer multiplicity of psi1redxB->psi2red (int=1)
    """
    #Defining the U3 variables 
    matrix=0.0
    if psi1red.check() and psi2red.check():
      s=psi1red.s
      n2=psi2red.n
      p2=psi2red.r
      w2=psi2red.w
      n1=psi1red.n
      p1=psi1red.r
      w1=psi1red.w
      K2,n2itr=kmatrix(s,w2)
      K1,n1itr=kmatrix(s,w1)
      K1inv=K1**(-1)
      i2=n2itr[n2,p2]
      j1=n1itr[n1,p1]
      for nb,rb in n1itr:
        i1=n1itr[nb,rb]
        psi1bred=SpStateReduced(s,nb,rb,w1)
        rbpmax=u3.mult(s,nb,w2)
        for rbp in range(1,rbpmax+1):
          j2=n2itr[nb,rbp]
          psi2bred=SpStateReduced(s,nb,rbp,w2)
          matrix=matrix+K2[i2,j2]*u3boson.C_s2(psi2bred,psi1bred,r0)*K1inv[i1,j1]
      if abs(matrix)<(10**(-13)):
        matrix=0.0
    return matrix

if __name__ == '__main__':
  s=U3State([8,4,4])
  n=U3State([2,2,0])
  r=1
  np=U3State([2,2,0])
  w=U3State([10,6,4])
  wp=U3State([10,6,4])
  psi1red=SpStateReduced(s,n,r,w)
  psi2red=SpStateReduced(s,np,r,wp)
  print "C_s1           ", C_s1(psi1red,psi2red,2)
  print "C_s2           ", C_s2(psi1red,psi2red,2)
  print "C_s2_asymptotic", C_s2(psi1red,psi2red,2)

# ################################################################################################################
def T(psip,psi,r0=1):
  """
  Kinetic energy operator in symplectic basis
  """  
  Kinetic=0
  if psip.w.N==psi.w.N:
    Kinetic=psi.w.N+3./2
  if psip.w.N==psi.w.N+2:
    Kinetic=-numpy.sqrt(3./2)*A(psip,psi,1)
  if psip.w.N==psi.w.N-2:
    Kinetic=-numpy.sqrt(3./2)*B(psip,psi,1)
  return Kinetic

if __name__ == '__main__':
  s=U3State([8,4,4])
  n=U3State([2,2,0])
  r=1
  np=U3State([2,2,0])
  w=U3State([10,6,4])
  wp=U3State([10,6,4])
  psi1red=SpStateReduced(s,n,r,w)
  psi2red=SpStateReduced(s,np,r,wp)
  print "T           ", T(psi1red,psi2red,1)
  np=U3State([4,2,0])
  wp=U3State([12,6,4])
  psi2red=SpStateReduced(s,np,r,wp)
  print "T           ", T(psi1red,psi2red,1)
  np=U3State([2,0,0])
  wp=U3State([8,6,4])
  psi2red=SpStateReduced(s,np,r,wp)
  print "T           ", T(psi1red,psi2red,1)




