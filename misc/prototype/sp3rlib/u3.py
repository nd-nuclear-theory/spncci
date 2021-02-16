##################################################################################
# U3 and SU(3) coefficients and functions
#   Anna E. McCoy
#   University of Notre Dame
#
#   SPDX-License-Identifier: MIT
##################################################################################

import numpy
import sys
import os
#sys.path.append('/afs/crc.nd.edu/user/a/amccoy/python/su3lib/')
sys.path.append('/Users/annamccoy/Dropbox/Research/Sp(6,R)/Code/su3lib/')
sys.path.append('/Users/annamccoy/Dropbox/Research/Sp(6,R)/Code/')
import su3lib
su3lib.blocks()
from u3states import *
import utilities as utils
import so3

#####################################################################################################################################################################################################


########################################################################################################################################################################################
def mult(w1,w2,w3):
    """
    Calculates the outer multipicity of the SU(3) coupling of w1xw2 to w3.  Calls su3lib library
    
    args:
        w1,w2,w3 (U3State or SU3State class)
    returns 
        mult (integer)
    """
    mult=su3lib.multu3(w1.su3.lbda,w1.su3.mu,w2.su3.lbda,w2.su3.mu,w3.su3.lbda,w3.su3.mu)
    return mult

def mult2(w1,w2,w3):
    """
    Calculate outer multiplicity without calling su3lib
    """
    multu3=0
    Nx=w1.su3.lbda+w2.su3.lbda-w3.su3.lbda-w1.su3.mu-w2.su3.mu+w3.su3.mu
    Mx=Nx/3
    if Nx==Mx*3:
        if Mx>=0:
            L1=w1.su3.lbda
            L2=w2.su3.lbda
            L3=w3.su3.lbda
            M1=w1.su3.mu
            M2=w2.su3.mu
            M3=w3.su3.mu
        else:
            L1=w1.su3.mu
            L2=w2.su3.mu
            L3=w3.su3.mu
            M1=w1.su3.lbda
            M2=w2.su3.lbda
            M3=w3.su3.lbda
            Mx=-Mx
        Nx=Mx+M1+M2-M3
        Mu=min(L1-Mx,M2)
        if Mu>=0:
            Nu=min(L2-Mx,M1)
            if Nu>=0:
                multu3=max(min(Nx,Nu)-max(Nx-Mu,0)+1,0)
    return multu3

if __name__ == '__main__':
    print "mult"
    w1=SU3State(5,4)
    w2=SU3State(1,1)
    w3=U3State([12,7,3])
    print mult(w1,w2,w3),mult2(w1,w2,w3)
    w1=SU3State(5,4)
    w2=SU3State(1,1)
    w3=U3State([12.5,7.5,3.5])
    print mult(w1,w2,w3), mult2(w1,w2,w3)

# #####################################################################################################################################################################################################
def couple(w1,w2):
    """
    Generates a dictionary of SU(3) coupled states of '(labda1,mu1)' and '(lbda2,mu2)' with their multiplicities. 
        
    Args: 
        w1,w2  (SU3State or U3State class)
        canonical: if True, order output to be in canonical ordering and return a list of tuples (w,rmax) rather than dictionary

    Returns: 
        A dict of the SU(3) coupled states with their mutiplicities  {SU3(lbda,mu):rho}.  
    
    For example:
        lm1=SU3State(2,1)
        lm2=SU3State(1,1)
        couple_SU3()   
        returns 
            {SU3State(3, 2): 1, SU3State(1, 3): 1, SU3State(2, 1): 2, SU3State(1, 0): 1, SU3State(0,2):1, SU3State(4,0):1}
    """
    #Extracting the SU(3) content
    lbda1=w1.su3.lbda
    mu1=w1.su3.mu
    lbda2=w2.su3.lbda
    mu2=w2.su3.mu
    lm1=w1.su3
    lm2=w2.su3
    SU3_coupled={}
    canonicallist=[]
    #Special case
    if lm1==SU3State(0,0) or lm2==SU3State(0,0):
        SU3_coupled[SU3State(lbda1+lbda2,mu1+mu2)]=1
    else:
        for k in range(0,min(mu2,lbda1+mu1)+1):
            for j in range(min(mu1,lbda2,lbda1+mu1-k)+1):
                for i in range(max(0,j-mu1+k),min(lbda2-j+k,lbda1)+1):
                    lbda=lbda1+lbda2-j-2*i+k
                    mu=mu1+mu2+i-j-2*k
                    w=SU3State(lbda,mu)
                    if w in SU3_coupled:
                        SU3_coupled[w]=SU3_coupled[w]+1
                    else:
                        SU3_coupled[w]=1              
    return SU3_coupled



def canonical_sort(su3list):
    switch=True
    while switch:
        switch=False
        for i in range(len(su3list)-1):
            if su3list[i].su3>su3list[i+1].su3:
                tmp=su3list[i+1]
                su3list[i+1]=su3list[i]
                su3list[i]=tmp
                switch=True
    return su3list

    # for w,rmax in sets2:
    #     print w,rmax
#come back and check by hand
if __name__ == '__main__':
    print "couple"
    sets=couple(SU3State(2,1),SU3State(1,1))
    for lm in sets:
        print lm,sets[lm]
    print canonical_sort(sets.keys())

##########################################################################################


def U(w1,w2,w,w3,w12,r12,r12_3,w23,r23,r1_23):
    """
    Racah unitary coefficient. 
    input:
        w1,w2,w,w2,w12,w23 are U3State or SU3State class 
        r12,r12_3,r23,r1_23 are integers 
    output:
        Coef: float
    E
    """
    r12max=mult(w1,w2,w12)
    r12_3max=mult(w12,w3,w)
    r23max=mult(w2,w3,w23)
    r1_23max=mult(w1,w23,w)
    #checking that all couplings are allowed
    rmax=r12max*r12_3max*r23max*r1_23max
    #check that function is calling for valid r
    r=r12*r12_3*r23*r1_23
    if (rmax==0) or (r==0):
        Coef=0
    else:
        index=r12+r12max*(r12_3-1)+r12max*r12_3max*(r23-1)+r12max*r12_3max*r23max*(r1_23-1)-1
        #index=r1_23+r1_23max*(r12-1)+r1_23max*r12max*(r12_3-1)+r1_23max*r12max*r12_3max*(r23-1)-1
        coef=su3lib.wru3optimized(w1.su3.lbda,w1.su3.mu,w2.su3.lbda,w2.su3.mu,w.su3.lbda,w.su3.mu,w3.su3.lbda,w3.su3.mu,w12.su3.lbda,w12.su3.mu,w23.su3.lbda,w23.su3.mu,r12max,r12_3max,r23max,r1_23max,rmax)
        try: 
            Coef=coef[index]
        except IndexError: 
            print w1,w2,w,w3,w12,r12,r12_3,w23,r23,r1_23, "index=",index,"indexmax=",len(coef)
    if abs(Coef)<10**(-13):
        Coef=0.0
    return Coef

if __name__ == '__main__':
    print "U"
    print U(SU3State(6,4), SU3State(2,0), SU3State(5,5), SU3State(2,0), SU3State(7,3),1,1, SU3State(2,1),1,1)
    print U(SU3State(6,4), SU3State(2,0), SU3State(5,5), SU3State(2,0), SU3State(7,3),1,1, SU3State(2,1),1,2)
    print 
    print U(SU3State(6,4), SU3State(2,0), SU3State(5,5), SU3State(2,0), SU3State(4,6),1,1, SU3State(4,0),1,1)
    a=numpy.sqrt((SU3State(4,6).dim()*SU3State(4,0).dim())/(SU3State(2,0).dim()*SU3State(5,5).dim()))
    print U(SU3State(4,6), SU3State(4,6), SU3State(4,0), SU3State(2,0), SU3State(2,0),1,1, SU3State(5,5),1,1)*a
    print a
    print 
    print U(SU3State(5,2), SU3State(1,4), SU3State(5,4), SU3State(4,0), SU3State(6,3),1,1, SU3State(0,2),1,1)

    print mult(SU3State(1,4), SU3State(4,0), SU3State(0,2))
    #36(5,2) (1,4) 34(5,4) (4,0) (6,3) 1 1 (0,2) 2 1 
    #(1,4) (4,0) (0,2) (0,0) (0,2) 1 1 (4,0) 1 2
    # print U(SU3State(6,4), SU3State(2,0), SU3State(5,5), SU3State(2,0), SU3State(4,6),1,1, SU3State(0,2),1,1)
    # print U(SU3State(6,4), SU3State(2,0), SU3State(5,5), SU3State(2,0), SU3State(4,6),1,1, SU3State(2,1),1,1)
    # print U(SU3State(6,4), SU3State(2,0), SU3State(5,5), SU3State(2,0), SU3State(4,6),1,1, SU3State(2,1),1,2)
    # print 
    # print U(SU3State(6,6), SU3State(0,4), SU3State(6,6), SU3State(4,0), SU3State(8,3),1,1, SU3State(4,4),1,1)
    # print U(SU3State(6,6), SU3State(0,4), SU3State(6,6), SU3State(4,0), SU3State(8,3),1,1, SU3State(4,4),1,2)
    # print U(SU3State(6,6), SU3State(0,4), SU3State(6,6), SU3State(4,0), SU3State(8,3),1,1, SU3State(4,4),1,3)
    # print U(SU3State(6,6), SU3State(0,4), SU3State(6,6), SU3State(4,0), SU3State(8,3),1,1, SU3State(4,4),1,4)
    # print U(SU3State(6,6), SU3State(0,4), SU3State(6,6), SU3State(4,0), SU3State(8,3),1,1, SU3State(4,4),1,5)
    # print 
    # print U(SU3State(9,3), SU3State(1,2), SU3State(8,2), SU3State(2,1), SU3State(8,3),1,1, SU3State(2,2),1,1)
    # print U(SU3State(9,3), SU3State(1,2), SU3State(8,2), SU3State(2,1), SU3State(8,3),1,1, SU3State(2,2),1,2)
    # print U(SU3State(9,3), SU3State(1,2), SU3State(8,2), SU3State(2,1), SU3State(8,3),1,1, SU3State(2,2),2,1)
    # print U(SU3State(9,3), SU3State(1,2), SU3State(8,2), SU3State(2,1), SU3State(8,3),1,1, SU3State(2,2),2,2)
    # print U(SU3State(9,3), SU3State(1,2), SU3State(8,2), SU3State(2,1), SU3State(8,3),1,2, SU3State(2,2),1,1)
    # print U(SU3State(9,3), SU3State(1,2), SU3State(8,2), SU3State(2,1), SU3State(8,3),1,2, SU3State(2,2),1,2)
    # print U(SU3State(9,3), SU3State(1,2), SU3State(8,2), SU3State(2,1), SU3State(8,3),1,2, SU3State(2,2),2,1)
    # print U(SU3State(9,3), SU3State(1,2), SU3State(8,2), SU3State(2,1), SU3State(8,3),1,2, SU3State(2,2),2,2)
    # print U(SU3State(9,3), SU3State(1,2), SU3State(8,2), SU3State(2,1), SU3State(8,3),2,1, SU3State(2,2),1,1)
    # print U(SU3State(9,3), SU3State(1,2), SU3State(8,2), SU3State(2,1), SU3State(8,3),2,1, SU3State(2,2),1,2)
    # print U(SU3State(9,3), SU3State(1,2), SU3State(8,2), SU3State(2,1), SU3State(8,3),2,1, SU3State(2,2),2,1)
    # print U(SU3State(9,3), SU3State(1,2), SU3State(8,2), SU3State(2,1), SU3State(8,3),2,1, SU3State(2,2),2,2)
    # print U(SU3State(9,3), SU3State(1,2), SU3State(8,2), SU3State(2,1), SU3State(8,3),2,2, SU3State(2,2),1,1)
    # print U(SU3State(9,3), SU3State(1,2), SU3State(8,2), SU3State(2,1), SU3State(8,3),2,2, SU3State(2,2),1,2)
    # print U(SU3State(9,3), SU3State(1,2), SU3State(8,2), SU3State(2,1), SU3State(8,3),2,2, SU3State(2,2),2,1)
    # print U(SU3State(9,3), SU3State(1,2), SU3State(8,2), SU3State(2,1), SU3State(8,3),2,2, SU3State(2,2),2,2)
    print 

def Z(w1,w2,w,w3,w12,r12,r12_3,w23,r23,r1_23):
    """
    Racah Z coefficient
    Args: w1,w2,w,w3,w12,w23 are SU3State or U3State instances 
    Returns: Coef, float
    """
    r12max=mult(w1,w2,w12)
    r12_3max=mult(w12,w3,w)
    r23max=mult(w2,w3,w23)
    r1_23max=mult(w1,w23,w)
    rmax=r12max*r12_3max*r23max*r1_23max
    r=r12*r12_3*r23*r1_23
    if (rmax>0)and (r>0):
        index=r12+r12max*(r12_3-1)+r12max*r12_3max*(r23-1)+r12max*r12_3max*r23max*(r1_23-1)-1
        #index=r1_23+r1_23max*(r12-1)+r1_23max*r12max*(r12_3-1)+r1_23max*r12max*r12_3max*(r23-1)-1
        coef=su3lib.wzu3optimized(w1.su3.lbda,w1.su3.mu,w2.su3.lbda,w2.su3.mu,w.su3.lbda,w.su3.mu,w3.su3.lbda,w3.su3.mu,w12.su3.lbda,w12.su3.mu,w23.su3.lbda,w23.su3.mu,r12max,r12_3max,r23max,r1_23max,rmax)
    #print coef
        Coef=coef[index]
    else:
        Coef=0
    return float(Coef)
if __name__ == '__main__':
    print "Z"
    print Z(SU3State(2,0),SU3State(0,0),SU3State(4,0),SU3State(2,0),SU3State(2,0),1,1,SU3State(2,0),1,1)
    print 
##################################################################################################   
def phi(w1,w2,w3,r,rp):
    """
    Phi phase factor that arrises in changing the coupling order of SU(3) irreps 
    """
    phi=0
    mlt=mult(w1,w2,w3)
    if (mlt<r) or (mult<rp):
        print "invalid multiplicity"
        return phi
    if mlt==1:            
        phi=utils.parity(w1.pt()+w2.pt()-w3.pt())
    else:
        phi=Z_SU3(w1,SU3State(0,0),w3,w2,w1,1,r,w2,1,rp)
    return phi
if __name__ == '__main__':
    print "phi"
    w1=U3State([2,0,0])
    w2=U3State([1,0,-1])
    w3=U3State([2,0,0])
    print phi(w2,w1,w3,1,1)
    print 

##################################################################################################
def coupling_9lm(w1,w2,w12,r12,w3,w4,w34,r34,w13,w24,w,r13_24,r13,r24,r12_34):
    """
    Unitary 9lm symbol
    """
    
    lm9=0.0
    # w13=SU3xU1_to_U3(lbda13+2*mu13,lbda13,mu13)
    # w2=SU3xU1_to_U3(lbda2+2*mu2,lbda2,mu2)
    w132set=couple(w13,w2)
    for w123 in w132set: 
        r12_3max=mult(w12,w3,w123)
        r13_2max=mult(w13,w2,w123)
        r04max=mult(w123,w4,w)
        for r12_3 in range(1,r12_3max+1):
            for r13_2 in range(1,r13_2max+1):
                for r04 in range(1,r04max+1):
                    #print w12,w3,w,w4,w123,r12_3,r04,w34,r34,r12_34
                    u1=U(w13,w2,w,w4,w123,r13_2,r04,w24,r24,r13_24)
                    z=Z(w2,w1,w123,w3,w12,r12,r12_3,w13,r13,r13_2)
                    u2=U(w12,w3,w,w4,w123,r12_3,r04,w34,r34,r12_34)
                    #print w123,r12_3,r13_2,r04,u2,z,u2
                    lm9=lm9+u1*z*u2
    return lm9

# w1=SU3State(2,0)
# w2=SU3State(1,1)
# w12=SU3State(3,1)
# w3=SU3State(1,1)
# w4=SU3State(0,2)
# w34=SU3State(1,3)
# w13=SU3State(3,1)
# w24=SU3State(1,3)
# w=SU3State(0,0)

if __name__ == '__main__':
    print "lm9"
    w1=SU3State(0,2)
    w2=SU3State(2,0)
    w12=SU3State(0,0)
    w3=SU3State(1,1)
    w4=SU3State(2,0)
    w34=SU3State(3,1)
    w13=SU3State(0,2)
    w24=SU3State(4,0)
    w=SU3State(3,1)



    print coupling_9lm(
        SU3State(2,0), SU3State(0,2), SU3State(0,0), 1,
        SU3State(2,1), SU3State(1,1), SU3State(2,1), 1,
        SU3State(0,3), SU3State(1,3), SU3State(2,1), 1,
        1, 1, 1)

    # print coupling_9lm(w1,w2,w12,1,
    #     w3,w4,w34,1,
    #     w13,w24,w,1,
    #     1,1,1), numpy.sqrt(w13.dim()/w1.dim()/w3.dim())*U(w13,w1.conj(),w,w4,w3,1,1,w24,1,1)
    # print 
    ###need to test after testing Z_SU3
##############################################################################################
def Wreduced_su3_dic(w1,l1,w2,l2,w3,l3):
    """
    Calculates dictionary of reduced SU(3) coupling coefficient for all rho, k1,k2,k3,. 

    Args:
        lbda1: SU(3) quantum number of (lbda1,mu1)
        mu1: SU(3) quantum number of (lbda1,mu1)
        l1: SO(3) quantum number in reduction chain of SU(3) state (lbda1,mu1)
        lbda2: SU(3) quantum number of (lbda2,mu2)
        mu2: SU(3) quantum number of (lbda2,mu2)
        l2: SO(3) quantum number in reduction chain of SU(3) state (lbda2,mu2)
        lbda3: SU(3) quantum number of (lbda3,mu3)
        mu3: SU(3) quantum number of (lbda3,mu3)
        l3: SO(3) quantum number in reduction chain of SU(3) state (lbda3,mu3)

    returns:
        coef: dictionary of coupling coefficients with key k0,k1,k2,k3, where:
            k0 is SU(3) outer multiplicity of (lbda1,mu1)x(lbda2,mu2)->(lbda3,mu3)
            ki, i=1..3 are SU(3)->SO(3) braching multiplicities: (lbdai,mui)->li
            For example, u3r3(4,3,3,2,0,2,5,2,2) returns: 
                {(1, 2, 1, 1): -0.68567430587075617, (1, 1, 1, 1): -0.32432095153102564}

    """
    k0max=mult(w1,w2,w3) #SU(3) outer multiplicity of the (lbda1,mu1)x(lbda2,mu2)->(lbda3,mu3) coupling
    k1max=so3.multk(w1,l1)
    k2max=so3.multk(w2,l2)
    k3max=so3.multk(w3,l3)
    dummy=numpy.float64(1) #dummy variable to satisfy Fortran code.  Not used in Fotran code
    coef={}
    #Note that k1max, k2max and k3max are calculated in the Fortran code.  Here they just serve as dummy variables 
    coefm=su3lib.wu3r3w(w1.su3.lbda,w1.su3.mu,w2.su3.lbda,w2.su3.mu,w3.su3.lbda,w3.su3.mu,l1,l2,l3,k0max,k1max,k2max,k3max,dummy)
    #coefm is the (9,9,9,9) array returns by fortran.  Dimensions cannot be adjusted. The following converts output into dictionary with key
        #k0,k1,k2,k3
    for k0 in range(1,k0max+1):
        for k1 in range(1,k1max+1):
            for k2 in range(1,k2max+1):
                for k3 in range(1,k3max+1):
                    coef[k0,k1,k2,k3]=coefm[k0-1,k1-1,k2-1,k3-1]
    #print k3max
    return coef

# ################################################################################################################################

#print "hi"
#wu39lm(la,ma,lb,mb,lc,mc,ld,md,le,me,lf,mf,lg,mg,lh,mh,li,mi,s,n9lmax,e__er)
#print su3lib.wu39lm(1,1,0,0,1,1,0,0,2,0,2,0,1,1,2,0,2,0,1)[0]



########################
def W(w1,k1,L1,w2,k2,L2,w3,k3,L3,r,c1=0,c2=0,c3=0):
    """
    Calculates the SO(3) reduced SU(3) wigner coupling coefficients
    Note 
    
    """
    W=0.0
    if (c1+c2+c3)==0:
        wdic=Wreduced_su3_dic(w1,L1,w2,L2,w3,L3)
        ######add conjugate
    #print wdic
        if (r,k1,k2,k3) in wdic:
            W=wdic[r,k1,k2,k3]
    if (c1,c2,c3)==(1,0,0):
        rpmax=mult(w1,w2,w3)
        wdic=Wreduced_su3_dic(w3,L3,w1.conj(),L1,w2,L2)
        for rp in range(1,rpmax+1):
            W=W+(utils.parity(w2.pt()-w3.pt())
                *phi(w1,w2,w3,r,rp)
                *numpy.sqrt(w3.dim()/w2.dim())
                *wdic[rp,k3,k1,k2])
    if (c1,c2,c3)==(0,1,0):
        wdic=Wreduced_su3_dic(w3,L3,w2.conj(),L2,w1,L1)
        if(r,k3,k2,k1) in wdic:
            W=(utils.parity(w1.pt()-w3.pt()+L1+L2-L3)
                *numpy.sqrt(w3.dim()/w1.dim())
                *wdic[r,k3,k2,k1])
    if (c1,c2,c3)==(0,0,1):
        print "unfinished code"
    if (c1+c2+c3)>=2:
        print "unfinished code"
    if abs(W)<10**(-10):
            W=0.0

    return W
if __name__ == '__main__':
    print "W"
    print W(SU3State(2,0),1,2,SU3State(0,2),1,2,SU3State(0,0),1,0,1),W(SU3State(2,0),1,2,SU3State(0,2),1,2,SU3State(0,0),1,0,1,c2=1),numpy.sqrt(2*(2+1))
    print W(SU3State(4,3),1,2,SU3State(2,0),1,2,SU3State(5,2),1,2,1)
    print W(SU3State(4,3),1,3,SU3State(2,0),1,2,SU3State(5,2),1,2,1)
    print W(SU3State(4,3),2,3,SU3State(2,0),1,2,SU3State(5,2),1,2,1)
    print "basis testing"
    print W(SU3State(1,0),1,1,SU3State(1,0),1,1,SU3State(0,1),1,1,1)
    print W(SU3State(0,1),1,1,SU3State(1,1),1,2,SU3State(0,1),1,1,1)**2
    # print W(SU3State(4,3),1,2,SU3State(2,0),1,2,SU3State(5,2),1,3,1)
    # print W(SU3State(4,3),1,2,SU3State(2,0),1,2,SU3State(5,2),2,3,1)
    # print W(SU3State(4,3),1,3,SU3State(2,0),1,2,SU3State(5,2),1,3,1)
    # print W(SU3State(4,3),2,3,SU3State(2,0),1,2,SU3State(5,2),1,3,1)
    # print W(SU3State(4,3),1,3,SU3State(2,0),1,2,SU3State(5,2),2,3,1)
    # print W(SU3State(4,3),2,3,SU3State(2,0),1,2,SU3State(5,2),2,3,1)

    # print W(SU3State(8,5),2,3,SU3State(4,4),2,2,SU3State(8,5),2,3,1)
    # print W(SU3State(8,5),2,3,SU3State(4,4),2,2,SU3State(8,5),2,3,2)
    # print W(SU3State(8,5),2,3,SU3State(4,4),2,2,SU3State(8,5),2,3,3)
    # print W(SU3State(8,5),2,3,SU3State(4,4),2,2,SU3State(8,5),2,3,4)
    # print W(SU3State(8,5),2,3,SU3State(4,4),2,2,SU3State(8,5),2,3,5)
    
####################################################################################
def casimir2_su3(w): 
  """
  Second order SU(3) Casimir operator
  inputs: 
    lbda: integer
    mu: integer 
  returns:
    matrix: value of the 2nd order Casimir operator, float
  """
  lm=w.su3
  matrix=2./3.*(lm.lbda**2+lm.lbda*lm.mu+lm.mu**2+3*lm.lbda+3*lm.mu)
  return matrix
# if __name__ == '__main__':
#     print "casimir2_su3"
#     w=SU3State(4,0)
#     print C_su3(w)**2, casimir2_su3(w)


####################################################################################
def C_su3(w):
  """
  SU3 reduced matix element of the SU3 generator C^(1,1)
  """
  lm=w.su3
  matrix=0
  lm_max=max(lm.lbda,lm.mu)
  if lm_max>0:
      l=lm_max%2
      if l==0:
        l=l+2 
      sign=numpy.sign(W(lm,1,l,SU3State(1,1),1,1,lm,1,l,1))
      matrix=sign*numpy.sqrt(casimir2_su3(w))
      #matrix=numpy.sqrt(lm.lbda**2+lm.mu**2+lm.lbda*lm.mu+3*lm.lbda+3*lm.mu)
      #matrix=numpy.sqrt(4*(lbda**2+lbda*mu+mu**2+3*lbda+3*mu))
      #matrix=2*matrix*sign
  return matrix
if __name__ == '__main__':
    print "C_su3"
    w=SU3State(4,0)
    print C_su3(w), numpy.sqrt(4*(16+12))


####################################################################################
def casimir3_su3(w): 
  """
  Third order SU(3) Casimir operator
  inputs: 
    lbda: integer
    mu: integer 
  returns:
    matrix: value of the 3rd order Casimir operator, float
  """
  lm=w.su3
  matrix=1./9.*(lm.su3.lbda-lm.su3.mu)*(lm.su3.lbda+2*lm.su3.mu+3)*(2*lm.su3.lbda+lm.su3.mu+3)
  return matrix
if __name__ == '__main__':
    print "casimir3_su3"
    w=SU3State(4,0)
    print casimir3_su3(w)



if __name__ == '__main__':
    ### Checking orthogonality
    coupled=0
    w1=SU3State(4,0)
    w2=SU3State(0,2)
    L1set=so3.branching_so3(w1)
    L2set=so3.branching_so3(w2)
    for L1 in L1set:
        k1max=L1set[L1]
        for L2 in L2set:
            k2max=L2set[L2]
            for k1 in range(1,k1max+1):
                for k2 in range(1,k2max+1):
                    coupled=coupled+W(w1,k1,L1,w2,k2,L2,SU3State(4,2),2,2,1)*W(w1,k1,L1,w2,k2,L2,SU3State(4,2),2,2,1)

    print "coupled",coupled 
