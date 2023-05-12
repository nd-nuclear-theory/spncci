      SUBROUTINE IJFACT(LM,MU,I,J,IFIX,JFIX)                            
C     ------------------------------------------------------------------
C     SUBPROGRAM FOR FIXING THE (I,J) FACTORS IN THE SU3 -> R3 SCHEME   
C     ------------------------------------------------------------------
      I=IFIX                                                            
      J=JFIX                                                            
      IF(LM.GE.MU)RETURN                                                
      I=1-IFIX                                                          
      J=1-JFIX                                                          
      RETURN                                                            
      END                                                               
