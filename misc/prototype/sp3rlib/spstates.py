##################################################################################
# Classes for Sp(3,R) states, reduced states 
# Operator class for SU(3) tensor operators
# 
#   Anna E. McCoy
#   University of Notre Dame
#
#   SPDX-License-Identifier: MIT
##################################################################################


from u3states import *
import u3
#########################################################################################################
#class definitions
class SpStateReduced(U3State):
    """
    unbranched symplectic state |snrw>
    s,n,w: U3State class
    r: int
    """
    def __init__(self,s,n,r,w):
        self.s=s
        self.n=n          
        self.r=r
        self.w=w            

    def __repr__(self):
        """
        print form
        """
        return "[{},{},{},{}]".format(self.s, self.n, self.r, self.w)
    
    def U3(self):
        return (self.s,self.n,self.r,self.w)
    def check(self):
        rmax=u3.mult(self.s,self.n,self.w)
        check=False
        if self.r in range(1,rmax+1):
            if self.s.check() and self.n.check() and self.w.check():
                if self.n.N%2==0:
                    check=True
        return check

    def __hash__(self): 
        return hash(self.U3())

    def __eq__(self,other):
        return self.U3()==other.U3()
        ## add raise exception for mismatch class types 

    def __ne__(self,other):
        return not self.__eq__(other)
    
    def __gt__(self,other):
        """
        Ordering defined by s then n then r then w
        """
        gt=False
        #ordering by s
        if self.s>other.s:
            gt=True
        else:
            ## if s's are the same
            if self.s==other.s:
                ## check n
                if self.n>other.n:
                    gt=True
                else:
                    #if s and n are the same
                    if self.n==other.n:
                        #check r
                        if self.r>other.r:
                            gt=True
                        else:
                            # if s,n,r are the same check w
                            if self.r==other.r and self.w>other.w:
                                gt=True
        return gt 

    def __lt__(self,other):
        return not (self==other or self>other)


# #########################################################################################################

class SpState(SpStateReduced):
    def __init__(self,s,n,r,w,k,L):
        self.s=s
        self.n=n
        self.r=r
        self.w=w
        self.k=k
        self.L=L
        
    def red(self):
        return SpStateReduced(self.s,self.n,self.r,self.w)

    def check(self):
        check=False
        if self.red().check():
            kmax=u3.multk(self.w,self.L)
            if self.k in range(1,kmax+1):
                check=True
        return check
    def __repr__(self):
        return "[{}{}{}{}{}{}]".format(self.s, self.n,self.r,self.w,self.k,self.L)

    def __hash__(self):
        return hash((self.red(),self.k,self.L))

    def __eq__(self,other):
        return (self.red(),self.k,self.L)==(other.red(),other.k,other.L)
        
    def __ne__(self,other):
        return not self.__eq__(other)

    def __gt__(self,other):
        gt=False
        ## compare reduced states 
        if self.red()>other.red():
            gt=True
        else:
            ## if reduced states are the same 
            if self.red()==self.ref():
                ## check branching mult.
                if self.k>other.k:
                    gt=True
                else:
                    ## if k's are equal, check L
                    if self.k==other.k and self.L>other.L:
                        gt=True
        return gt

    def __lt__(self,other):
        return not (self==other or self>other)

if __name__ == '__main__':

    s=U3State([6,4,4])
    n=U3State([4,0,0])
    r=1
    w=U3State([10,4,4])
    np=U3State([2,2,0])
    rp=1
    wp=U3State([9,5,4])
    psired=SpStateReduced(s,n,r,w)
    psipred=SpStateReduced(s,np,rp,wp)
    psi1=SpState(s,n,r,w,1,2)
    psi2=SpState(s,np,rp,wp,1,2)
    print psi1.check()

    #print psi1==psi2

##########################################################################################################################################
class SpU2StateRed:
    """
    symplectic state with S 
    spred is SpStateReduced
    """
    def __init__(self,spred,S):
        self.sp=spred
        self.S=S

    def __repr__(self):
        return "[{}{}]".format(self.sp,self.S)
    

    def __eq__(self,other):
        return (self.sp,self.S)==(other.sp,other.S)
        
    def __ne__(self,other):
        return not self.__eq__(other)

    def __gt__(self,other):
        gt=False
        if self.sp.s > other.sp.s:
            gt=True
        else:
            if self.sp.s==other.sp.s and self.S>other.S:
                gt=True
            else:
                if self.S==other.S and self.sp>other.sp:
                    gt=True
        return gt

    def __lt__(self,other):
        return not (self==other or self>other)

    def __hash__(self):
        return hash((self.sp,self.S))

    def Sp(self):
        return SpState(self.sp.s,self.sp.n,self.sp.r,self.sp.w)


##########################################################################################################################################
class SpU4StateRed:
    """
    symplectic state with S 
    spred is SpStateReduced
    """
    def __init__(self,spred,S,T):
        self.sp=spred
        self.S=S
        self.T=T

    def __repr__(self):
        return "[{}{}]".format(self.sp,self.S,self.T)
    

    def __eq__(self,other):
        if isinstance(other,self.__class__):
            return (self.sp,self.S,self.T)==(other.sp,other.S,other.T)
        else:
            return False    ####maybe change to error message

    def __ne__(self,other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.sp,self.S,self.T))

##########################################################################################################################################


class SpU2State:
    """
    formerly Sp3RxU4State
    symplectic state with S
    spred is SpStateReduced
    """
    def __init__(self,spred,k,L,S):
        self.sp=spred
        self.k=k
        self.L=L
        self.S=S

    def __repr__(self):
        return "[{}{}{}{}{}]".format(self.sp,self.k,self.L,self.S)
    

    def red(self):
        return SpU2StateRed(self.sp,self.S)

    def __eq__(self,other):
        return (self.sp,self.k,self.L,self.S)==(other.sp,other.k,other.L,other.S)
        
    def __ne__(self,other):
        return not self.__eq__(other)

    def __gt__(self,other):
        gt=False
        if self.sp.s>other.sp.s:
            gt=True
        else:
            if self.sp.s==other.sp.s and self.S>other.S:
                gt=True
            else:
                if self.S==other.S and self.sp>other.sp:
                    gt=True
        return gt


    def __hash__(self):
        return hash((self.sp,self.k,self.L,self.S))

    def SpkL(self):
        return SpState(self.sp.s,self.sp.n,self.sp.r,self.sp.w,self.k,self.L)

##########################################################################################################################################
class SpU4State:
    """
    formerly Sp3RxU4State
    symplectic state with S and T
    spred is SpStateReduced
    """
    def __init__(self,spred,k,L,S,T):
        self.sp=spred
        self.k=k
        self.L=L
        self.S=S
        self.T=T

    def __repr__(self):
        return "[{}{}{}{}{}]".format(self.sp,self.k,self.L,self.S,self.T)
    

    def __eq__(self,other):
        return (self.sp,self.k,self.L,self.S,self.T)==(other.sp,other.k,other.L,other.S,other.T)
        
    def __ne__(self,other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.sp,self.k,self.L,self.S,self.T))

    def SpkL(self):
        return SpState(self.sp.s,self.sp.n,self.sp.r,self.sp.w,self.k,self.L)
##########################################################################################################################################
class SpU2StateJ:
    def __init__(self,spred,k,L,S,J):
        self.sp=spred
        self.k=k
        self.L=L
        self.S=S
        self.J=J

    def __repr__(self):
        """
        print form 
        """
        return "[{},{},{},{},{}]".format(self.sp,self.k,self.L,self.S,self.J)

    def spu2(self):
        """
        LS coupled Sp3RxU2 state 
        """
        return SpU2State(self.sp,self.k,self.L,self.S)

    def spkL(self):
        """
        branched symplectic states 
        """
        return SpState(self.sp.s,self.sp.n,self.sp.r,self.sp.w,self.k,self.L)

    def red(self):
        return SpU2StateRed(self.sp,self.S)

    def __eq__(self,other):
        return (self.sp,self.k,self.L,self.S,self.J)==(other.sp,other.k,other.L,other.S,other.J)
        
    def __ne__(self,other):
        return not self==other

    def __gt__(self,other):
        gt=False
        if self.spu2()>other.spu2():
            gt=True
        else:
            if self.spu2()==other.spu2() and self.J>other.J:
                gt=True
        return gt

    def __lt__(self,other):
        return not (self==other or self>other)
    def __hash__(self):
        return hash((self.sp,self.k,self.L,self.S,self.J))

##########################################################################################################################################
class SpU4StateJ:
    def __init__(self,spred,k,L,S,J,T):
        self.sp=spred
        self.k=k
        self.L=L
        self.S=S
        self.J=J
        self.T=T

    def __repr__(self):
        """
        print form 
        """
        return "[{},{},{},{},{},{}]".format(self.sp,self.k,self.L,self.S,self.J,self.T)

    def spu4(self):
        """
        LS coupled Sp3RxU4 state 
        """
        return SpU4State(self.sp,self.k,self.L,self.S,self.T)

    def spkL(self):
        """
        branched symplectic states 
        """
        return SpState(self.sp.s,self.sp.n,self.sp.r,self.sp.w,self.k,self.L)

    def __eq__(self,other):
        return (self.sp,self.k,self.L,self.S,self.J,self.T)==(other.sp,other.k,other.L,other.S,other.J,other.T)
         
################################################################################################################################
class Sp3RSU3Operator:
    def __init__(self,w0,RME,su3Conserve):
        self.RME=RME
        self.u3=w0
        self.su3Conserved=su3Conserve

    def SU3RME(self,psi2red,psi1red,r):
        ME=0.0
        if (self.su3Conserved==True) and (r!=1):
            return ME
        rmax=u3.mult(psi1red.w,self.u3,psi2red.w)
        if r>rmax:
            return ME
        if psi1red.s==psi2red.s:
            ME=self.RME(*(psi2red,psi1red,r))
        return ME


# #########################################################################################################    
if __name__ == '__main__':

    s=U3State([6,4,4])
    n=U3State([4,0,0])
    r=1
    w=U3State([10,6,4])
    np=U3State([4,0,0])
    rp=1
    wp=U3State([8,6,4])
    psired=SpStateReduced(s,n,r,w)
    psipred=SpStateReduced(s,np,rp,wp)
    print psired.check()
    print psired==psipred
    print psipred




