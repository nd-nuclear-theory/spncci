      REAL*8 FUNCTION DELTA(J1T,J2T,J3T)
C     ------------------------------------------------------------------
C     DELTA FOR R3 ROUTINES--TRIANGLE RELATIONS CHECKED                 
C     ------------------------------------------------------------------
      IMPLICIT REAL*8(D)                                                
      COMMON/BKDF/DLOGF(0:2000)                                         
      DELTA=12345.D0                                                    
      I1=J1T+J2T-J3T                                                    
      IF(BTEST(I1,0))GOTO 10                                            
      IF(I1.LT.0)GOTO 10                                                
      I2=J2T+J3T-J1T                                                    
      IF(I2.LT.0)GOTO 10                                                
      I3=J3T+J1T-J2T                                                    
      IF(I3.LT.0)GOTO 10                                                
      DELTA=DLOGF(I1)+DLOGF(I2)+DLOGF(I3)-DLOGF(J1T+J2T+J3T+2)          
   10 RETURN                                                            
      END                                                               
