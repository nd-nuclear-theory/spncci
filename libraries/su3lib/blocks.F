      SUBROUTINE BLOCKS                                                 
C     ------------------------------------------------------------------
C     BINOMIAL COEFFICIENTS AND FACTORIALS ***** SEE COMMENT BELOW *****
C     ------------------------------------------------------------------
C     UPDATE/MOD: (MTS,06-76)  H.SATO             EXPANDED RANGE        
C                 (LSU,11-78)  J.P.DRAAYER        LOG FACTORIALS        
C                 (LSU,05-80)  J.P.DRAAYER        LOG BINOMIALS         
C                 (LSU,08-81)  J.P.DRAAYER        EXPANDED RANGE        
C                 (LSU,01-83)  J.P.DRAAYER        EXTENDED PRECISION    
C                 (LSU,03-83)  J.P.DRAAYER        D,Q MIX & LOG CUTS    
C                 (LSU,11-84)  J.P.DRAAYER        MODIFIED FACTORIALS   
C                 (LSU,01-88)  J.P.DRAAYER        DLOGF RANGE/INDEX     
C                 (LSU,10-89)  J.P.DRAAYER        BINOMIAL INVERSES     
C                 (LSU,11-89)  J.P.DRAAYER        POWERS OF TWO ARRAY   
C                 (ND,11-11)   M.A.CAPRIO         DOUBLE/QUAD SWITCHABLE
C                 (LSU,11-11)  T. DYTRYCH         MERGE WITH SU3GENBK
C                                                                       
C     BKDB--BINOMIAL (BINO) COEFFICIENTS (EXPANDED 6-76,8-81,1-83)      
C           SCALE: BINO(I,J)=DBINO(I*(I+1)/2+J+1)                       
C           RANGE: BINO(0,0)=DBINO(1) TO BINO(128,128)=DBINO(8385)      
C           ADDED: 2**I = DTWOS(I) WHERE I=-128,128                     
C     BKQB--BINOMIAL (BINO) COEFFICIENTS (EXPANDED 6-76,8-81,1-83)      
C           SCALE: BINO(I,J)=QBINO(I*(I+1)/2+J+1)                       
C           RANGE: BINO(0,0)=QBINO(1) TO BINO(192,192)=QBINO(18721)     
C           ADDED: 2**I = QTWOS(I) WHERE I=-192,192                     
C     BKDF--LOG FACTORIALS (FACT) (INSERTED 11-78, MODIFIED 01-88)      
C           SCALE: LNFACT(I)=DLOGF(2*I)                                 
C           RANGE: LNFACT(0)=DLOGF(0) TO LNFACT(1000)=DLOGF(2000)       
C     ------------------------------------------------------------------
      IMPLICIT REAL*8(D)
#if defined(SU3DBL)
      IMPLICIT REAL*8(Q)
#elif defined(SU3QUAD)
      IMPLICIT REAL*16(Q)                                     
#endif
      COMMON/BKDB/DBINO(8385),DBINV(8385),DTWOS(-128:128)               
      COMMON/BKQB/QBINO(18721),QBINV(18721),QTWOS(-192:192)             
      COMMON/BKDF/DLOGF(0:2000)                                         
      DLOGF(0)=0.D0
C                                                                       
      DO 10 I=2,2000,2                                                  
10       DLOGF(I)=DLOGF(I-2)+DLOG(DFLOAT(I/2))                          
C                                                                       
      DO 30 M=1,193                                                     
         I=M-1                                                          
      DO 30 N=1,M                                                       
         J=N-1                                                          
         L=MIN0(J,I-J)                                                  
         QGEN=1.0                                                      
      IF (L.EQ.0)GO TO 30                                               
C                                                                       
         DO 20 K=1,L
#if defined(SU3DBL)
20          QGEN=DFLOAT(I+1-K)/DFLOAT(K)*QGEN                           
#elif defined(SU3QUAD)
20          QGEN=QFLOAT(I+1-K)/QFLOAT(K)*QGEN                           
#endif
C                                                                       
30       QBINO(I*(I+1)/2+J+1)=QGEN                                      
C                                                                       
      DO 40 I=1,18721                                                   
40       QBINV(I)=1.0/QBINO(I)                                         
C                                                                       
      DO 50 I=1,8385                                                    
         DBINO(I)=QBINO(I)                                              
50       DBINV(I)=QBINV(I)                                              
C                                                                       
         QTWOS(0)=1.0                                                  
      DO 60 I=1,192                                                     
         QTWOS(+I)=QTWOS(I-1)*2.0                                      
60       QTWOS(-I)=QTWOS(1-I)/2.0                                      
C                                                                       
      DO 70 I=-128,128                                                  
70       DTWOS(I)=QTWOS(I)                                              
C                                                                       
      RETURN                                                            
      END                                                               
