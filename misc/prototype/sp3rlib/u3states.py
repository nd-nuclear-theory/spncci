##Classes for U3 and SU(3) states 

##################################################################################################################
class SU3State:
    """
    state with SU(3) symmetry, not branched to L
    lbda,mu: positive int
    """
    def __init__(self,x,y):
        self.lbda=x
        self.mu=y
        self.su3=self
        
    def __repr__(self):
        return "({},{})".format(self.lbda, self.mu)

    def __eq__(self,other):
        if isinstance(other,self.__class__):
              return (self.lbda,self.mu)==(other.lbda,other.mu)
        else:
              return False    ####maybe change to error message

    def __ne__(self,other):
        return not self.__eq__(other)   

    def __hash__(self): 
        return hash((self.lbda,self.mu))

    def __gt__(self,other):
        """
        defines a canonical ordering of SU(3) irreps, first by lbda then by mu
        """
        gt=False
        if self.lbda>other.lbda:
            gt=True
        else:
            if (self.lbda==other.lbda) and (self.mu>other.mu):
                gt=True
        return gt

    def __lt__(self,other):
        return (not self.__gt__(other)) and (not self==other)

    def __ge__(self,other):
        return not self<other
    def __le__(self,other):
        return not self>other

    def dim(self):
        """
        returns dimension of SU(3) irrep
        """
        d=(self.lbda+1)*(self.mu+1)*(self.lbda+self.mu+2)/2.
        return d

    def conj(self):
        """
        returns the conjugate of SU(3) irrep, e.g. (lbda,mu)->(mu,lbda)
        """
        return SU3State(self.mu,self.lbda)

    def pt(self):
        """
        returns exponenet of (-1), e.g. (-1)^w=(-1)^(lbda+mu)
        """
        #this is not the usual parity but the number that goes into the parity function 
        p=self.lbda+self.mu
        return p

###########################################################################################
class U3State:
    def __init__(self,x,y=False):
        """
        Defines a U(3) basis state.  
        Args: takes either N,SU3state(lambda,mu) or [w1,w2,w3] as an argument, 
                but stores internally in the Elliot notation N(lambda,mu)
        """
        #determine which form of U3 var is input
        if y!=False:
            self.N=x
            self.su3=y
        else: 
            self.N=x[0]+x[1]+x[2]
            self.su3=SU3State(int(x[0]-x[1]),int(x[1]-x[2]))
    
    ###################################################
    ## Converting to cartesian representation
    ###############################################
    def three(self):
        three=(self.N-self.su3.lbda-2*self.su3.mu)/3.
        if int(2*three)!=2*three:
            print "raise exception", self, three
        else:
            return three

    def one(self):
        return self.su3.lbda+self.su3.mu+self.three()

    def two(self):
        return self.su3.mu+self.three()

    def cartesian(self):
        """
        returns list of cartesian labels 
        """
        return [self.one(), self.two(), self.three()]
    ###############################################

    def __repr__(self):
        return "{}{}".format(self.N, self.su3)

    def check(self):
        """
        checks if variable satisfies the U(3) condition that w1>=w2>=w3>=0
        """
        ck=False
        if (self.one()>=self.two()>=self.three()>=0):
        # t3=self.three()
        # if self.N==0:
        #     if self.su3==SU3State(0,0):
        #         ck=True
        # elif t3%3==0 and t3>=0:
            ck=True
        return ck

    def __eq__(self,other):
        if isinstance(other,self.__class__):
              return (self.N,self.su3)==(other.N,other.su3)
        else:
              return False    ####maybe change to error message

    def __ne__(self,other):
        return not self.__eq__(other)

    def __gt__(self,other):
        """
        defines a canonical ordering of U(3) states, based on canonical ordering of SU(3) states  
        """
        gt=False
        if self.N>other.N:
            gt=True
        else:
            if (self.N==other.N) and self.su3>other.su3:
                gt=True
        return gt

    def __lt__(self,other):
        lt=False
        if (not self>other) and (not self==other):
            lt=True
        return lt

    def __ge__(self,other):
        return not self<other

    def __le__(self,other):
        return not self>other

    def __hash__(self): 
        return hash((self.N,self.su3))

    def dim(self):
        """
        returns the dimension of U(3) irrep, same as SU(3) dim, since dim(U(1))=1
        """
        return self.su3.dim()

    def conj(self):
        """
        returns conjugate of U(3) state 
        """
        return U3State(-self.N, self.su3.conj())

    def pt(self):
        """
        returns exponenet of (-1) in (-1)^w
        """
        return self.su3.pt()
#########################################################################################################
#####################################################################################################################
class U3SU2State:
    def __init__(self,w,S):
        """
        Defines a U(3)xSU(2) basis state.  
        Args: 
            w: U3State
            S: spin (int or half int)
        """
        #determin which form of U3 var is input
        self.w=w
        self.S=S
    
    def dim(self):
        dim=self.w.dim()*(2*self.S+1)
        return dim
    def __repr__(self):
        return "[{},{}]".format(self.w,self.S)
    
    def __eq__(self,other):
        if isinstance(other,self.__class__):
              return (self.w,self.S)==(other.w,other.S)
        else:
              return False    ####maybe change to error message

    def __ne__(self,other):
        return not self.__eq__(other)


    def __hash__(self): 
        return hash((self.w,self.S))
        
    def pt(self):
        return self.w.pt()+self.S



#########################################################################################################
class U3SU4State:
    def __init__(self,w,S,T):
        """
        Defines a U(3)xSU(2)xSU(2) basis state with labels w S and T.  
        Args: 
            w: U3State class
            S: int or half int
            T: int or half int
        """
        #determin which form of U3 var is input
        self.w=w
        self.S=S
        self.T=T
    
    def dim(self):
        """
        dimension of U(3)xSU(2)xSU(2) irrep
        """
        dim=self.w.dim()*(2*self.S+1)*(2*self.T+1)
        return dim

    def __repr__(self):
        return "[{},{},{}]".format(self.w,self.S,self.T)
    
    def __eq__(self,other):
        if isinstance(other,self.__class__):
              return (self.w,self.S,self.T)==(other.w,other.S,other.T)
        else:
              return False    ####maybe change to error message

    def __ne__(self,other):
        return not self.__eq__(other)

    def __gt__(self,other):
        """
        defining greater than condition, ording by w, then S, then T
        """
        gt=False
        if self.w>other.w:
            gt=True
        else:
            if (self.w==other.w) and self.S>other.S:
                gt=True
            else:
                if (self.w==other.w) and (self.S==other.S) and (self.T>other.T):
                    gt=True
        return gt

    def __lt__(self,other):
        """
        define less than condition
        """
        lt=False
        ## if self is not equal to other and self is not great than other, self is less than other
        if (not self==other) and (not self>other):
            lt=True
        return lt

    def __ge__(self,other):
        return not self<other

    def __le__(self,other):
        return not self>other

    def __hash__(self): 
        return hash((self.w,self.S,self.T))

    def pt(self):
        return self.w.pt()+self.S+self.T

    def conj(self):
        return U3SU4State(self.w.conj(),self.S, self.T)




if __name__ == '__main__':
    # print "SU3State testing"
    lm1=SU3State(4,5)
    lm2=SU3State(7,2)
    lm3=SU3State(4,5)
    # print lm1          #(4,5)
    # print lm1==lm3          #True
    # print lm2==lm3          #False
    # print lm2!=lm3          #True
    # print lm2.dim()         #132
    # print lm1.conj()   #(5,4)
    # print lm1.pt()          #1
    # print [lm1, lm2, lm3]
    # print lm1<lm2
    # print lm2>lm1
    # print lm1>lm3

    # print "U3State testing"
    # w1=U3State(29,lm1)
    # w2=U3State(32,lm3)
    # w3=U3State(0,SU3State(10,4))
    # print w3.check()
    # print w1.check()            #True
    # print w1.one()              #14
    # print w1.two()              #10
    # print w1.three()            #5
    # print w1                    #29(4,5)
    # print w1.cartesian()        #[14,10,5]
    # print w1.su3                #(4,5)
    # print w2==w1                #False
    # print w2!=w1                #True
    # print w2.su3==w1.su3        #True
    # print w2.conj()             #-32(5,4)
    # print w2>w1
    # print w2.su3>w1.su3
    # print w2.su3<w1.su3

    print "half int"
    w1h=U3State(30.5,lm1)
    w2h=U3State(33.5,lm3)
    print w1h.three()
    print w1h.one()               #14
    print w1h.two()               #10
    print w1h.three()
    print w1h.check()             #True
    print w1h.one()               #14
    print w1h.two()               #10
    print w1h.three()             #5
    print w1h                     #30.5(4,5)
    print w1h.cartesian()         #[14,10,5]
    print w1h.su3                 #(4,5)
    print w2h==w1h                #False
    print w2h!=w1h                #True
    print w2h.su3==w1h.su3        #True
    print w2h.conj()              #-32(5,4)

    print "U3xSU4State testing "
    w=SU3State(4,0)
    wp=SU3State(5,0)
    state1=U3SU4State(w,1,0)
    state2=U3SU4State(w,1,1)
    state3=U3SU4State(wp,1,0)
    state4=U3SU4State(w,1,0)
    print state1
    print state2
    print state1==state2
    print state1==state3
    print state1==state4


