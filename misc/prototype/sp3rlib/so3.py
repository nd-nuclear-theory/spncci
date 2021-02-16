##################################################################################
#   Anna E. McCoy
#   University of Notre Dame
#
#   SPDX-License-Identifier: MIT
##################################################################################


import math
import sys
import os
sys.path.append('/Users/annamccoy/Dropbox/Research/Sp(6,R)/Code/')
from u3states import *
import pygsl.sf
import utilities as utils

def mult(J1,J2,J3):
	mult=0
	if abs(J1-J2)<=J3<=(J1+J2):
		mult=1
	return mult
if __name__ == '__main__':
	print mult(1,2,3)
	print mult(1,2,4)
	print mult(1,1,2)

###################################################################################################

def hat(J):
    """
    Calculates the sqrt(2*J+1) factor from angular momentum theory
    args:
        J (int or half int)
    returns:
        hat(J) (float)
    """
    hat=math.sqrt(2*J+1)
    return hat

###################################################################################################
def couple_so3(L1,L2,orbital=False):
    """
    Couples to angular momenta together according to the triangle rule 
    args:
        L1, L2 are orbital angular momenta, only takes integers 
    returns:
        List of total orbital angular momenta
    """
    imax=int(2*min(L1,L2))
    if orbital==True:
        if (L1+L2)%1==0:
            L3set=[int(L1+L2-i) for i in range(0,imax+1)]
        else:
            print "error message"
            L3set=[]
    else:
        L3set=[L1+L2-i for i in range(0,imax+1)] 
    return L3set
if __name__ == '__main__':
    print "couple_so3"
    print couple_so3(2,1,True)

################################################################################################################################################################################################################

def clebsch_gordan((J1,M1),(J2,M2),(J3,M3)):
    """
    Calculates the SO(3) coupling coefficient, or clebsch-gordan coefficient.  Calls pygsl special functions library.

    args:
        J1,J2,J3 are angualr momentum (int or half int)
        M1,M2,M3 are angular momentum projections (int or half int)

    """
    cg=(-1)**(J2-J1-M3)*math.sqrt(2*J3+1)*pygsl.sf.coupling_3j(int(2*J1), int(2*J2), int(2*J3), int(2*M1), int(2*M2),int(-2*M3))[0]
    #
    return cg

if __name__ == '__main__':
    print "clebsch_gordan"
    print clebsch_gordan((.5,.5),(.5,.5),(1,1)), 1.0
    print clebsch_gordan((.5,.5),(.5,-.5),(1,0)), 1/math.sqrt(2)
    print clebsch_gordan((.5,.5),(.5,-.5),(0,0)), 1/math.sqrt(2)
    print clebsch_gordan((.5,-.5),(.5,.5),(1,0)), 1/math.sqrt(2)
    print clebsch_gordan((.5,-.5),(.5,.5),(0,0)), -1/math.sqrt(2)
    print clebsch_gordan((.5,-.5),(.5,-.5),(1,-1)), 1.0
    print clebsch_gordan((3./2,3./2),(1,0),(3./2,3./2)), math.sqrt(3./5)
    print 
########################################################################################################################################################################################
def coupling_6j((J1,J2,J3),(J4,J5,J6)):
    """
    Unitary 6j symbol. Calls pygsl special functions library.

    Args:
        J1, J2, J3, J4, J5, J6 (int or float) : angular momenta (not twice)


    Returns:
        (float) : value of unitary 6j sybol
    """
    cj6=pygsl.sf.coupling_6j(int(2*J1),int(2*J2), int(2*J3), int(2*J4),int(2*J5), int(2*J6))[0]
    cj6=(-1)**(J1+J2+J4+J5)*hat(J3)*hat(J6)*cj6 
    return cj6 
if __name__ == '__main__':
    print "coupling_6j", coupling_6j((0,.5,.5),(.5,1,1))
########################################################################################################################################################################################
def coupling_9j((J1,J2,J3), (J4,J5,J6), (J7,J8,J9)):
    """
    Unitary 9j symbol.  Calls pygsl special functions library. 

    args: 
        J1, J2, J3, J4, J5, J6 are angular momenta (int or half int)
    """
    cj9=pygsl.sf.coupling_9j(
        int(2*J1),int(2*J2), int(2*J3),
        int(2*J4),int(2*J5), int(2*J6),
        int(2*J7),int(2*J8), int(2*J9)
        )[0]
    cj9=cj9*hat(J3)*hat(J6)*hat(J7)*hat(J8)
    if abs(cj9)<10**(-13):
        cj9=0.0
    return cj9
if __name__ == '__main__':
    print "coupling_9j",coupling_9j((1,1,1),(2,3,4),(3,4,5))
    print "coupling_9j",coupling_9j((0,.5,.5),(1,.5,.5),(1,1,0))#, coupling_6j((1,1,1),(3,2,3))
    print "coupling_9j",coupling_9j((1,1,1),(1,1,1),(2,2,0)), coupling_6j((1,1,1),(1,1,2))*utils.parity(1+1+1)
    print "coupling_9j",coupling_9j((1,1,0),(1,1,0),(2,2,0)), coupling_6j((1,1,0),(1,1,2))*utils.parity(1+1+0),1.*hat(2)/hat(1)/hat(1)
    print "basis testing"
    print "coupling_9j",coupling_9j((1,.5,.5),(1,.5,.5),(0,0,0))**2
    print "coupling_9j",coupling_9j((1,.5,.5),(1,.5,.5),(1,1,0))**2
    print "coupling_9j",coupling_9j((1,1,0),(0,0,0),(1,1,0))**2
    print "coupling_9j",coupling_9j((1,1,0),(1,1,0),(1,1,0))**2
    print "coupling_9j",coupling_9j((1,1,0),(2,2,0),(1,1,0))**2
    print "test coupling_9j",coupling_9j((2,0,2),(1,0,1),(1,0,1))





############################################################################################################################################################################################
def multJ(J,M):
    """
    Calculates the branching multiplicity of J to M.  Multiplicity is 1 if M if valid projection and 0 if not. 
    """
    mult=0
    if abs(M)<=J:
        #check if both are half integers or integers
        if (2*J)%2==(2*M)%2:
            mult=1
    return mult

 ############################################################################################################################################################################
def multk(w,l):
    """
    Calculates branching multiplicity for SU(3) to SO(3) (kappa) for a given l. 

    Args:
        lm: SU3State instance  
        l: SO(3) quantum number.  

    returns:
        mk: SU(3)->SO(3) branching multiplicity.  For example,
            multk(SU3State(4,3),3)
            returns:
                2 
    """
    lm=w.su3
    ##will need to change for python 3
    mk=max(0,int((lm.lbda+lm.mu+2.-l)/2))-max(0,int((lm.lbda+1.-l)/2))-max(0,int((lm.mu+1.-l)/2))
    return mk
if __name__ == '__main__':
    print "multk"
    print multk(SU3State(4,3),3)
    print 

############################################################################


def branching_so3(w,S=0, J=0, Restrict=False):
    """
    Determines the allowed SO(3) irreps (L) in a a given SU(3) irrep and their branching multiplicity kappa accordin to the branching rule:
    kappa=mubar,mubar-2,...,1 or 0    mubar=min(lbda,mu)
    
    L=
        kappa,kappa+1,...,lbdabar          if kappa !=0         lbdabar=max(lbda,mu)
        lbdabar, lbdabar-2,...,1 or 0      if kappa=0
        kappamax=max(...)+max(...)+max(...) formula

    Args:
        w: SU(3) or U(3) irrep labels, SU3State class

    Returns:
        A dictionary {L:kappa}.  
        For example, 
            branching_so3(SU3State(4,5))
        returns 
            {1: 1, 2: 1, 3: 2, 4: 2, 5: 3, 6: 2, 7: 2, 8: 1, 9: 1}
    """
    Lkset={}
    lm=w.su3
    lbdabar=max(lm.lbda,lm.mu)
    if min(lm.lbda,lm.mu)==0: 
        for L in range(lbdabar,-1,-2):
            kmax=multk(lm,L)
            if kmax>0:
                Lkset[L]=kmax
    else:
        for L in range(lm.lbda+lm.mu+1):
            kmax=multk(lm,L)
            if kmax>0:
                Lkset[L]=kmax
    if Restrict==True:
        Lset=couple_so3(S,J)
        Lintersect=Lset&Lkset.viewkeys()
        restrictset={L:Lkset[L] for L in Lintersect}
        Lkset=restrictset
    return Lkset
if __name__ == '__main__':
    "print branching_so3"
    w=SU3State(4,5)
    print branching_so3(SU3State(4,5))
    print "restricted to S=1, J=0", branching_so3(w,S=0, J=2, Restrict=True)




