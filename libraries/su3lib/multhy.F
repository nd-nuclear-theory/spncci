      FUNCTION MULTHY(L1,M1,L2,M2,L3,M3)                                
C     ------------------------------------------------------------------
C     MULTIPLICITY (THEORY) IN U3 COUPLING                              
C     ------------------------------------------------------------------
      DIMENSION IX(6)                                                   
      MULTHY=0                                                          
      IX(1)=L1+L2-L3+2*(M1+M2-M3)                                       
      IX(2)=M1+M2-M3+2*(L1+L2-L3)                                       
      IX(3)=2*L2+M2+M1-L1-M3+L3                                         
      IX(4)=2*M2+L2+L1-M1-L3+M3                                         
      IX(5)=L3+M2-L1+2*(M3+L2-M1)                                       
      IX(6)=M3+L2-M1+2*(L3+M2-L1)                                       
      IXMIN=1000                                                        
      DO 10 I=1,6                                                       
      IXDB3=IX(I)/3                                                     
      IF(3*IXDB3.LT.IX(I))RETURN                                        
      IF(IXDB3.LT.IXMIN)IXMIN=IXDB3                                     
   10 CONTINUE                                                          
      IF(IXMIN.LT.0)RETURN                                              
      MULTHY=MIN0(IXMIN,L2,M2)+1                                        
      RETURN                                                            
      END                                                               
