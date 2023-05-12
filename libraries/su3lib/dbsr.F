      SUBROUTINE DBSR(MA,DA,DB,DC,MD)                                   
C     ------------------------------------------------------------------
C     DOUBLE BACK SUBSTITUTION FOR SOLVING SIMULTANEOUS EQUATIONS       
C     ------------------------------------------------------------------
C     REFERENCES--COMPUTER SOLUTION OF LINEAR ALGEBRAIC SYSTEMS,        
C                 G.FORSYTHE AND C.MOLER, PRENTICE HALL                 
C     PARAMETERS--MA  ACTUAL NUMBER OF ROWS                             
C                 MD  DIMENSIONED NUMBER OF ROWS                        
C                 ND  DIMENSIONED NUMBER OF COLUMNS                     
C                 DA  MATRIX OF COEFFICIENTS PREPARED IN DLUT           
C                 DB  INPUT CONSTANTS REPLACED BY SOLUTIONS             
C                 DC  WORK VECTOR                                       
C     DIMENSIONS--DA(MD*ND),DB(ND),DC(ND)                               
C     ------------------------------------------------------------------
      IMPLICIT REAL*8(D)                                                
      DIMENSION DA(1),DB(1),DC(1)                                       
      DC(1)=DB(1)                                                       
      DO 15 I=2,MA                                                      
      NA=I-1                                                            
      DSUM=0.D0                                                         
      JQ=-MD                                                            
      DO 10 J=1,NA                                                      
      JQ=JQ+MD                                                          
      IJ=I+JQ                                                           
   10 DSUM=DSUM+DA(IJ)*DC(J)                                            
   15 DC(I)=DB(I)-DSUM                                                  
      MAMA=IJ+MD                                                        
      DB(MA)=DC(MA)/DA(MAMA)                                            
      IQ=NA*MD                                                          
      DO 25 IP=1,NA                                                     
      IQ=IQ-MD                                                          
      I=MA-IP                                                           
      II=I+IQ                                                           
      DSUM=0                                                            
      JQ=MA*MD                                                          
      DO 20 JP=1,IP                                                     
      JQ=JQ-MD                                                          
      J=MA-JP+1                                                         
      IJ=I+JQ                                                           
   20 DSUM=DSUM+DA(IJ)*DB(J)                                            
   25 DB(I)=(DC(I)-DSUM)/DA(II)                                         
      RETURN                                                            
      END                                                               
