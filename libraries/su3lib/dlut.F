      SUBROUTINE DLUT(MA,NA,IA,DA,MD)                                   
C     ------------------------------------------------------------------
C     DECOMPOSITION OF A REAL MATRIX INTO THE PRODUCT OF A LOWER AND    
C     AN UPPER TRIANGULAR MATRIX                                        
C     ------------------------------------------------------------------
C     REFERENCES--COMPUTER SOLUTION OF LINEAR ALGEBRAIC SYSTEMS,        
C                 G.FORSYTHE AND C.MOLER, PRENTICE HALL                 
C     PARAMETERS--MA  ACTUAL NUMBER OF ROWS                             
C                 NA  ACTUAL NUMBER OF COLUMNS                          
C                 MD  DIMENSIONED NUMBER OF ROWS                        
C                 ND  DIMENSIONED NUMBER OF COLUMNS                     
C                 DA  MATRIX TO BE DECOMPOSED                           
C                 IA  INPUT AS IA(I)=I RETURNED IN SORTED ORDER         
C     DIMENSIONS--DA(MD*ND),IA(MD)                                      
C     ------------------------------------------------------------------
      IMPLICIT REAL*8(D)                                                
      DIMENSION DA(1),IA(1)                                             
      LA=MIN0(MA,NA)                                                    
      LQ=-MD                                                            
      DO 35 L=1,LA                                                      
      IF(L.LT.LA)GO TO 10                                               
      IF(MA.LE.NA)RETURN                                                
   10 LQ=LQ+MD                                                          
      DBIG=0.D0                                                         
      DO 20 I=L,MA                                                      
      IL=I+LQ                                                           
      IF(DABS(DA(IL))-DABS(DBIG))20,20,15                               
   15 IBIG=I                                                            
      DBIG=DA(IL)                                                       
   20 CONTINUE                                                          
      ISAVE=IA(L)                                                       
      IA(L)=IA(IBIG)                                                    
      IA(IBIG)=ISAVE                                                    
      JQ=-MD                                                            
      DO 25 J=1,NA                                                      
      JQ=JQ+MD                                                          
      LJ=L+JQ                                                           
      IBIGJ=IBIG+JQ                                                     
      DSAVE=DA(LJ)                                                      
      DA(LJ)=DA(IBIGJ)                                                  
   25 DA(IBIGJ)=DSAVE                                                   
      K=L+1                                                             
      DO 30 I=K,MA                                                      
      IL=I+LQ                                                           
   30 DA(IL)=DA(IL)/DBIG                                                
      IF(L.EQ.LA)RETURN                                                 
      DO 35 I=K,MA                                                      
      IL=I+LQ                                                           
      JQ=LQ                                                             
      DO 35 J=K,NA                                                      
      JQ=JQ+MD                                                          
      IJ=I+JQ                                                           
      LJ=L+JQ                                                           
   35 DA(IJ)=DA(IJ)-DA(IL)*DA(LJ)                                       
      RETURN                                                            
      END                                                               
