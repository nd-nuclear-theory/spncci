      REAL*8 FUNCTION DRR3(J1T,J2T,L2T,L1T,J3T,L3T)
C     ------------------------------------------------------------------
C     RACAH COEFFICIENTS FOR R3--TRIANGLE RELATION CHECKED IN DELTA     
C     REFERENCES--THE 3-J AND 6-J SYMBOLS, M.ROTENBERG, R.BIVINS,       
C                 N.METROPOLIS AND J.K.WOOTEN, MIT PRESS                
C     ------------------------------------------------------------------
      IMPLICIT REAL*8(D)                                                
      COMMON/BKDF/DLOGF(0:2000)                                         
      DRR3=0.D0                                                         
      DX=DELTA(J1T,J2T,J3T)                                             
      IF(DX.EQ.12345D0)GOTO 20                                          
      DC=DX                                                             
      DX=DELTA(L1T,L2T,J3T)                                             
      IF(DX.EQ.12345D0)GOTO 20                                          
      DC=DX+DC                                                          
      DX=DELTA(L1T,J2T,L3T)                                             
      IF(DX.EQ.12345D0)GOTO 20                                          
      DC=DX+DC                                                          
      DX=DELTA(J1T,L2T,L3T)                                             
      IF(DX.EQ.12345D0)GOTO 20                                          
      DC=(DX+DC)/2.D0                                                   
      I1=J3T+L3T-J1T-L1T                                                
      I2=J3T+L3T-J2T-L2T                                                
      I3=J1T+J2T+L1T+L2T+2                                              
      I4=J1T+J2T-J3T                                                    
      I5=L1T+L2T-J3T                                                    
      I6=J1T+L2T-L3T                                                    
      I7=L1T+J2T-L3T                                                    
      ITMIN=MAX0(0,-I1,-I2)                                             
      ITMAX=MIN0(I3,I4,I5,I6,I7)                                        
      IF(ITMIN.GT.ITMAX)GOTO 20                                         
      DO 10 IT=ITMIN,ITMAX,2                                            
      DSUM=DEXP(DC+DLOGF(I3-IT)-(DLOGF(I4-IT)+DLOGF(I5-IT)+             
     1DLOGF(I6-IT)+DLOGF(I7-IT)+DLOGF(IT)+DLOGF(I1+IT)+DLOGF(I2+IT)))   
      IF(BTEST(IT,1))THEN                                               
      DRR3=DRR3-DSUM                                                    
      ELSE                                                              
      DRR3=DRR3+DSUM                                                    
      ENDIF                                                             
   10 CONTINUE                                                          
   20 RETURN                                                            
      END                                                               
