      REAL*8 FUNCTION DWR3(J1T,J2T,J3T,M1T,M2T,M3T)
C     ------------------------------------------------------------------
C     WIGNER COEFFICIENTS FOR R3--TRIANGLE RELATIONS CHECKED IN DELTA   
C     REFERENCES--ELEMENTARY THEORY OF ANGULAR MOMENTUM, M.E.ROSE, WILEY
C     ------------------------------------------------------------------
      IMPLICIT REAL*8(D)                                                
      COMMON/BKDF/DLOGF(0:2000)                                         
      DWR3=0.D0                                                         
      IF(M1T+M2T-M3T.NE.0)GOTO 20                                       
      DC=DELTA(J1T,J2T,J3T)                                             
      IF(DC.EQ.12345D0)GOTO 20                                          
      I1=J3T-J2T+M1T                                                    
      I2=J3T-J1T-M2T                                                    
      I3=J1T+J2T-J3T                                                    
      I4=J1T-M1T                                                        
      IF(BTEST(I4,0))GOTO 20                                            
      I5=J2T+M2T                                                        
      IF(BTEST(I5,0))GOTO 20                                            
      ITMIN=MAX0(0,-I1,-I2)                                             
      ITMAX=MIN0(I3,I4,I5)                                              
      IF(ITMIN.GT.ITMAX)GOTO 20                                         
      DTOP=(DLOG(DFLOAT(J3T+1))+DC+DLOGF(J1T+M1T)+DLOGF(J1T-M1T)+       
     1DLOGF(J2T+M2T)+DLOGF(J2T-M2T)+DLOGF(J3T+M3T)+                     
     2DLOGF(J3T-M3T))/DFLOAT(2)                                         
      DO 10 IT=ITMIN,ITMAX,2                                            
      DBOT=DLOGF(I3-IT)+DLOGF(I4-IT)+DLOGF(I5-IT)+                      
     1DLOGF(IT)+DLOGF(I1+IT)+DLOGF(I2+IT)                               
      DSUM=DEXP(DTOP-DBOT)                                              
      IF(BTEST(IT,1))THEN                                               
      DWR3=DWR3-DSUM                                                    
      ELSE                                                              
      DWR3=DWR3+DSUM                                                    
      ENDIF                                                             
   10 CONTINUE                                                          
   20 RETURN                                                            
      END                                                               
