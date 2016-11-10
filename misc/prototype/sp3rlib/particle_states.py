
from u3states import *
import so3
import u3
###########


class OBStateLS:
    """
    One-body state with LS symmetry and upstream label N
    """
    def __init__(self,N,l,s):
        self.N=N
        self.L=l
        self.S=s
        
    def __repr__(self):
        return "[{},{},{}]".format(self.N, self.L, self.S)
    
    def __eq__(self,other):
        if isinstance(other,self.__class__):
              return (self.N,self.L,self.S)==(other.N,other.L,other.S)
        else:
              return False    ####maybe change to error message

    def __ne__(self,other):
        return not self.__eq__(other)   

    def __hash__(self): 
        return hash((self.N,self.L,self.S))

if __name__ == '__main__':
    a1=OBStateLS(4,0,.5)
    a2=OBStateLS(4,2,1.5)
    print a1==a2
    print a1!=a2
    print a1.N, a1.L, a1.S

####################################################################################################################
class OBStateLST:
    """
    One-body state with LST symmetry and upstream lable N
    """
    def __init__(self,N,l,s,t):
        self.N=N
        self.L=l
        self.S=s
        self.T=t
        
    def __repr__(self):
        return "[{},{},{},{}]".format(self.N, self.L, self.S,self.T)
    
    def __eq__(self,other):
        if isinstance(other,self.__class__):
              return (self.N,self.L,self.S,self.T)==(other.N,other.L,other.S,other.T)
        else:
              return False    ####maybe change to error message

    def __ne__(self,other):
        return not self.__eq__(other)   

    def __hash__(self): return hash((self.N,self.L,self.S, self.T))

if __name__ == '__main__':
    b1=OBStateLST(4,0,.5,.5)
    b2=OBStateLST(4,2,1.5,.5)
    print b1==b2
    print b1==b1
    print b1!=b2
    print b2.N, b2.L, b2.S, b2.T

####################################################################################################################

class OBStateJ:
    """
    one-body state with symmetry J and upstream labels N,L,S
    """
    def __init__(self,N,L,S,J):
        self.N=N
        self.L=L
        self.S=S
        self.J=J
    def __repr__(self):
        return "[{},{},{},{}]".format(self.N, self.L,self.S,self.J)

    def __eq__(self,other):
        if isinstance(other,self.__class__):
            return (self.N, self.L,self.S,self.J)==(other.N,other.L,other.S,other.J)
        else:
            return False
    def __ne__(self,other):
        return not self.__eq__(other)
    def __hash__(self):
        return hash((self.N,self.L,self.S,self.J))
    def sp(self):
        if self.S==.5:
            return SingleState(self.N,self.L,self.J)
        else:
            print "not a single state "
####################################################################################################################
class OBStateJT:
    """
    one-body state with JT symmetry and known N,L and S upstream labels
    """
    def __init__(self,N,L,S,J,T):
        self.N=N
        self.L=L
        self.S=S
        self.J=J
        self.T=T
    def __repr__(self):
        return "[{},{},{},{},{}]".format(self.N, self.L,self.S,self.J,self.T)

    def __eq__(self,other):
        if isinstance(other,self.__class__):
            return (self.N, self.L,self.S,self.J,self.T)==(other.N,other.L,other.S,other.J,other.T)
        else:
            return False
    def __ne__(self,other):
        return not self.__eq__(other)
    def __hash__(self):
        return hash((self.N,self.L,self.S,self.J,self.T))
    def sp(self):
        if self.S==.5 and self.T==.5:
            return SingleState(self.N,self.L,self.J)
        else:
            print "not a single state "


####################################################################################################################
class SingleState:
    """
    single particle state
    """
    cutoff=20
    global H2_index
    H2_index={}
    i=1
    for N in range(cutoff+1):
        for L in range(N%2,N+1,2):
            jset=so3.couple_so3(L,.5)
            jset.sort()
            for J in jset:
                H2_index[N,L,J]=i
                i=i+1

    def __init__(self,N,l,j):
        self.N=N
        self.L=l
        self.J=j
        self.index=H2_index[N,l,j]

    def __repr__(self):
        return "[{},{},{}]".format(self.N, self.L, self.J)
    
    def __eq__(self,other):
        if isinstance(other,self.__class__):
              return (self.N,self.L,self.J)==(other.N, other.L, other.J)
        else:
              return False    ####maybe change to error message

    def __gt__(self,other):
        return (self.N,self.J)>(other.N, other.J)


    def __ne__(self,other):
        return not self.__eq__(other)   

    def __hash__(self): 
        return hash((self.N,self.L, self.J))


####################################################################################################################

class TBStateJ:
    """
    Two-body state J1xJ2->J coupling

    args:
        a1,a2 are single particle states (SingleState)
        J: total angular momentum (int)
        Type:
            if proton-proton: Type=11
            if neutron-neutron: Type=22
            if proton-neutron: Type=12
    """
    def __init__(self,a1,a2,J,Type):
        self.a1=a1
        self.a2=a2
        self.J=J
        self.Type=Type
    
    def __repr__(self):
        return "([{},{}],{},{})".format(self.a1,self.a2, self.J, self.Type)

    ## Defining state==stateb
    def __eq__(self,other):
        if isinstance(other,self.__class__):
            return (self.a1,self.a2,self.J, self.Type)==(other.a1, other.a2,other.J, other.Type)
        else:
            return False
    ## Defining state>stateb
    def __gt__(self,other):
        return (self.a1,self.a2)>(other.a1,other.a2)

    def __ne__(self,other):
        return not self.__eq__(other)

    def is_symmetry_allowed(self):
        sym=True
        if self.a1==self.a2 and self.J%2==1:
            sym=False
        return sym

    def is_canonical(self):
        return self.a1.index<=self.a2.index

    def h2(self):
        return "{:3} {:3} {:3} {:3}".format(self.a1.index, self.a2.index, self.J, self.Type)
####################################################################################################################

class TBStateJT:
    """
    Two-body state J1xJ2->J coupling and T symmetry

    args:
        a1,a2 are single particle states (SingleState)
        J: total angular momentum (int)
        T: isospin (int) can be 0 or 1

    """
    def __init__(self,a1,a2,J,T):
        self.a1=a1
        self.a2=a2
        self.J=J
        self.T=T
    
    def __repr__(self):
        return "([{},{}],{},{})".format(self.a1,self.a2, self.J, self.T)

    ## Defining state==stateb
    def __eq__(self,other):
        if isinstance(other,self.__class__):
            return (self.a1,self.a2,self.J, self.T)==(other.a1, other.a2,other.J, other.T)
        else:
            return False
    ## Defining state>stateb
    def __gt__(self,other):
        return (self.a1,self.a2)>(other.a1,other.a2)

    def __ne__(self,other):
        return not self.__eq__(other)

    def is_symmetry_allowed(self):
        sym=True
        if (self.a1==self.a2) and (self.J+self.T)%2==0:
            sym=False
        return sym

    def is_canonical(self):
        return self.a1.index<=self.a2.index

    def h2(self):
        return "{:3} {:3} {:3} {:3}".format(self.a1.index, self.a2.index, self.J, self.T)
#####################################################################################################################
class TBStateLS:
    """
    Two-body state LS coupled
     if proton-proton: Type=11
     if neutron-neutron: Type=22
     if proton-neutron: Type=12
     a1 and a2 are OBStateLS
    """
    def __init__(self,a1,a2,L,S,Type):
        self.a1=a1
        self.a2=a2
        self.L=L
        self.S=S
        self.Type=Type
    
    def __repr__(self):
        return "([{},{}],{},{},{})".format(self.a1,self.a2, self.L, self.S, self.Type)

    def __eq__(self,other):
        if isinstance(other,self.__class__):
            return (self.a1,self.a2,self.L,self.S,self.Type)==(other.a1, other.a2,other.L,other.S, other.Type)
        else:
            return False
#####################################################################################################################
class TBStateLST:
    """
    Two-body state LS-coupled 
     if Isopin T=1, Type=1
     if Isospin T=0: Type=0
     a1 and a2 are OBStateLST
    """
    def __init__(self,a1,a2,L,S,T):
        self.a1=a1
        self.a2=a2
        self.L=L
        self.S=S
        self.T=T
    
    def __repr__(self):
        return "([{},{}],{},{},{})".format(self.a1,self.a2, self.L, self.S, self.T)

    def __eq__(self,other):
        if isinstance(other,self.__class__):
            return (self.a1,self.a2,self.L,self.S,self.T)==(other.a1, other.a2,other.L,other.S, other.T)
        else:
            return False

    def is_symmetry_allowed(self):
        sym=True
        if self.a1==self.a2 and (self.L+self.S+self.L)%2==0:
            sym=False
        return sym






##################################################################################################

class TBStateWS:
    """
    Two-body state
    args:
        q1,q2: SU3State 
        state: U3SU2State
    """
    def __init__(self,q1,q2,state,Type):
        if isinstance(q1,SU3State):
            self.q1=q1
        else:
            self.q1=SU3State(q1,0)
        if isinstance(q2,SU3State):
            self.q2=q2
        else:
            self.q2=SU3State(q2,0)
        self.u3su2=state
        self.Type=Type
    
    def __repr__(self):
        return "[[{},{}] {} {}]".format(self.q1,self.q2, self.u3su2.w, self.u3su2.S,self.Type)

    def __eq__(self,other):
        if isinstance(other,self.__class__):
            return (self.q1,self.q2,self.u3su2)==(other.q1, other.q2,other.u3su2,self.Type)
        else:
            return False
    def __ne__(self,other):
        return not self.__eq__(other)

    def pt(self):
        """
        The 1 is coming from S1+S2=.5+.5=1
        """
        parity=self.q1.pt()+self.q2.pt()+self.u3su2.w.pt()+self.u3su2.S+1
        return parity

##################################################################################################
class TBStateWST:
    """
    Two-body state
    args:
        q1,q2: SU3State
        state: U3SU4State
    """
    def __init__(self,q1,q2,state):
        if isinstance(q1,SU3State):
            self.q1=q1
        else:
            self.q1=SU3State(q1,0)
        if isinstance(q2,SU3State):
            self.q2=q2
        else:
            self.q2=SU3State(q2,0)
        self.u3su4=state
    
    def __repr__(self):
        return "[[{},{}],{} {} {}]".format(self.q1,self.q2, self.u3su4.w, self.u3su4.S, self.u3su4.T)


    # def pr(self):
    #     return (self.a1, self.a2,self.J, self.Type)

    def __eq__(self,other):
        if isinstance(other,self.__class__):
            return (self.q1,self.q2,self.u3su4)==(other.q1, other.q2,other.u3su4)
        else:
            return False
    def __ne__(self,other):
        return not self.__eq__(other)

    def pt(self):
        parity=self.q1+self.q2+self.u3su4.w.pt()+self.u3su4.S+self.u3su4.T
        return parity

##################################################################################################



