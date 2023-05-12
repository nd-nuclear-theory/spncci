      SUBROUTINE CONMAT(I,J,LAM,MU,KAPMAX,L,DONMAT,N)                   
C     ------------------------------------------------------------------
C     ORTHONORMALIZATION COEFFICIENTS FOR THE U3R3 TRANSFORMATION       
C     ------------------------------------------------------------------
C     UPDATE/MOD: (LSU,05-83)  J.P.DRAAYER        EXTENDED PRECISION    
C                 (ND,11-11)   M.A.CAPRIO         DOUBLE/QUAD SWITCHABLE
C                 (ND,11-22)   A.E.MCCOY          ADDED FLAG SU3QUAD_GNU
C                                                                       
C     REFERENCES--J.P.DRAAYER AND Y.AKIYAMA, J.MATH.PHYS.14(1973)1904   
C                 J.P.DRAAYER, NUCL.PHYS.A129(1969)647-665              
C                 J.D.VERGADOS, NUCL.PHYS.A111(1968)681-754             
C     PARAMETERS--(I,J) : (1,1)=GHW, (1,0)=GHW', (0,0)=GLW, (0,1)=GLW'  
C                 N=MAX(KAPMAX)                                         
C     DIMENSIONS--DONMAT(N*N)                                           
C     ------------------------------------------------------------------
      IMPLICIT REAL*8(D)
#if defined(SU3DBL)
      IMPLICIT REAL*8(Q)
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
      IMPLICIT REAL*16(Q)
                                    
#endif
      DIMENSION DONMAT(1),QONMAT(81)                                    
      KSTART(LAM,MU,L)=MOD(LAM,2)+2*(MAX0(0,(L-MU)/2)+MOD(LAM+1,2)*     
     1MOD(IABS(L-MU),2))                                                
      IF(I.EQ.1)GO TO 10                                                
      KM=KSTART(MU,LAM,L)-2                                             
      IEE=2*LAM+MU                                                      
      JET=MU                                                            
      MET=-MU                                                           
      GO TO 15                                                          
   10 KM=KSTART(LAM,MU,L)-2                                             
      IEE=-(LAM+2*MU)                                                   
      JET=LAM                                                           
      MET=LAM                                                           
   15 IF(I.NE.J)MET=-MET                                                
      K=KM                                                              
      KKAQ=-N                                                           
      DO 35 KKA=1,KAPMAX                                                
      KKAQ=KKAQ+N                                                       
      K=K+2                                                             
      KKAKKA=KKA+KKAQ                                                   
      M=KM                                                              
      KKBQ=-N                                                           
      DO 30 KKB=1,KKA                                                   
      KKBQ=KKBQ+N                                                       
      M=M+2                                                             
      DHOLD=0.D0                                                        
      IF(KKB.EQ.1)GO TO 25                                              
      KB=KKB-1                                                          
      DO 20 KKC=1,KB                                                    
      KKCKKA=KKC+KKAQ                                                   
      KKCKKB=KKC+KKBQ                                                   
   20 DHOLD=DHOLD+DONMAT(KKCKKB)*DONMAT(KKCKKA)                         
   25 DHOLD=DTU3R3(I,J,LAM,MU,IEE,JET,MET,K,L,M)-DHOLD                  
      IF(KKB.EQ.KKA)GO TO 110                                           
      KKBKKA=KKB+KKAQ                                                   
      KKBKKB=KKB+KKBQ                                                   
   30 DONMAT(KKBKKA)=DONMAT(KKBKKB)*DHOLD                               
  110 IF(DHOLD.LE.0.D0)GO TO 115                                        
   35 DONMAT(KKAKKA)=DSQRT(1.D0/DHOLD)                                  
      KKAQ=-N                                                           
      DO 50 KKA=1,KAPMAX                                                
      KKAQ=KKAQ+N                                                       
      IF(KKA.EQ.1)GO TO 50                                              
      KKAKKA=KKA+KKAQ                                                   
      KA=KKA-1                                                          
      KKBQ=-N                                                           
      DO 45 KKB=1,KA                                                    
      KKBQ=KKBQ+N                                                       
      DHOLD=0.D0                                                        
      DO 40 KKC=KKB,KA                                                  
      KKCKKA=KKC+KKAQ                                                   
      KKCKKB=KKC+KKBQ                                                   
   40 DHOLD=DHOLD-DONMAT(KKCKKB)*DONMAT(KKCKKA)                         
      KKAKKB=KKA+KKBQ                                                   
      KKBKKA=KKB+KKAQ                                                   
      DONMAT(KKAKKB)=DONMAT(KKAKKA)*DHOLD                               
   45 DONMAT(KKBKKA)=-DONMAT(KKAKKA)*DONMAT(KKBKKA)                     
   50 CONTINUE                                                          
      RETURN                                                            
  115 K=KM                                                              
      KKAQ=-N                                                           
      DO 135 KKA=1,KAPMAX                                               
      KKAQ=KKAQ+N                                                       
      K=K+2                                                             
      KKAKKA=KKA+KKAQ                                                   
      M=KM                                                              
      KKBQ=-N                                                           
      DO 130 KKB=1,KKA                                                  
      KKBQ=KKBQ+N                                                       
      M=M+2                                                             
      QHOLD=0.D0                                                        
      IF(KKB.EQ.1)GO TO 125                                             
      KB=KKB-1                                                          
      DO 120 KKC=1,KB                                                   
      KKCKKA=KKC+KKAQ                                                   
      KKCKKB=KKC+KKBQ                                                   
  120 QHOLD=QHOLD+QONMAT(KKCKKB)*QONMAT(KKCKKA)                         
  125 QHOLD=QTU3R3(I,J,LAM,MU,IEE,JET,MET,K,L,M)-QHOLD                  
      IF(KKB.EQ.KKA)GO TO 135                                           
      KKBKKA=KKB+KKAQ                                                   
      KKBKKB=KKB+KKBQ                                                   
  130 QONMAT(KKBKKA)=QONMAT(KKBKKB)*QHOLD                               
#if defined(SU3DBL)
  135 QONMAT(KKAKKA)=DSQRT(1.D0/QHOLD)                                  
#elif defined(SU3QUAD)
  135 QONMAT(KKAKKA)=QSQRT(1.D0/QHOLD)                                  
#elif defined(SU3QUAD_GNU)
  135 QONMAT(KKAKKA)=SQRT(1.D0/QHOLD)                                  
#endif
      KKAQ=-N                                                           
      DO 150 KKA=1,KAPMAX                                               
      KKAQ=KKAQ+N                                                       
      IF(KKA.EQ.1)GO TO 150                                             
      KKAKKA=KKA+KKAQ                                                   
      KA=KKA-1                                                          
      KKBQ=-N                                                           
      DO 145 KKB=1,KA                                                   
      KKBQ=KKBQ+N                                                       
      QHOLD=0.D0                                                        
      DO 140 KKC=KKB,KA                                                 
      KKCKKA=KKC+KKAQ                                                   
      KKCKKB=KKC+KKBQ                                                   
  140 QHOLD=QHOLD-QONMAT(KKCKKB)*QONMAT(KKCKKA)                         
      KKAKKB=KKA+KKBQ                                                   
      KKBKKA=KKB+KKAQ                                                   
      QONMAT(KKAKKB)=QONMAT(KKAKKA)*QHOLD                               
      QONMAT(KKBKKA)=-QONMAT(KKAKKA)*QONMAT(KKBKKA)                     
      DONMAT(KKAKKB)=QONMAT(KKAKKB)                                     
  145 DONMAT(KKBKKA)=QONMAT(KKBKKA)                                     
  150 CONTINUE                                                          
      RETURN                                                            
      END                                                               
