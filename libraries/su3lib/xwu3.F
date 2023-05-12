      SUBROUTINE XWU3(LAM1,MU1,LAM2,MU2,LAM3,MU3,IE,JT,NEC,DEWU3,       
     1KR0MAX,INDMAX,DWU3,J1SMAX,J1TMAX,J2SMAX,J2TMAX,IESMAX,IE2MAX,     
     2INDMAT,N1,N2,N3,KIMAX2)                                           
C     ------------------------------------------------------------------
C     WIGNER COEFFICIENTS FOR U3 (X PREFIX FOR 6-81 VERSION)            
C     ------------------------------------------------------------------
C     UPDATE/MOD: (LSU,06-81)  J.P.DRAAYER        INDEXING OF DEWU3     
C                 (LSU,11-89)  J.P.DRAAYER        DWU3 ZERO-OUT RANGE   
C                                                                       
C     REFERENCES--J.P.DRAAYER AND Y.AKIYAMA, J.MATH.PHYS.14(1973)1904   
C                 K.T.HECHT, NUCL.PHYS.62(1965)1                        
C     PARAMETERS--SEE ALSO XEWU3                                        
C       EXTERNAL--N1=MAX(KR0MAX)                                        
C                 N2=MAX(IESMAX) SAFE TO SET N2=MAX(LAM2+MU2+1)         
C                 N3=MAX(DIM(LAM2,MU2))                                 
C                 NA=MAX(INDMAX)=NX*(NX+1)*(NX+2)/6, NX=MAX(LAM2+MU2+1) 
C                 NB=MAX(J2SMAX) SAFE TO SET NB=MAX(LAM2+MU2+1)         
C*                KIMAX1=MAX(KR0MAX*INDMAX) (SEE XEWU3)                 
C*                KIMAX2=MAX(KR0MAX*DIM(LAM2,MU2))                      
C       INTERNAL--X1=ABS(N1*N3)                                         
C                 X2=ABS(N2)                                            
C                 X3=ABS(N2*NB)                                         
C     EXTENSIONS--CHANGE EXTERNAL PARAMETERS IN CALL STATEMENT          
C                 ADJUST INTERNAL PARAMETERS BELOW                      
C*    DIMENSIONS--DEWU3(N1*NA),DWU3(N1*N3),J1SMAX(N2*NB),               
C                 J1TMAX(N2*NB),J2SMAX(N2),J2TMAX(N2),INDMAT(N2*NB),    
C                 DWU3P(X1),J2TMAP(X2),INDMAP(X3)                       
C*      COMMENTS--USE N1*NA->KIMAX1,N1*N3->KIMAX2                       
C                 ASSUME MAX N1=9,N2=42,N3=9030                         
C                        SET X1=27090,X2=42,X3=1764 (X1=3*N3,FIXED)     
C     ------------------------------------------------------------------
      IMPLICIT REAL*8(D),INTEGER(X)                                     
C     IMPLICIT INTEGER(X)                                               
      DIMENSION DEWU3(1),DWU3(1),J1SMAX(1),                             
     1          J1TMAX(1),J2SMAX(1),J2TMAX(1),INDMAT(1),                
     2          DWU3P(27090),J2TMAP(42),INDMAP(1764)                    
      INDEX(J1TD,LAM1,J1T,J2TD,LAM2,J2T)=1+J2TD*(J2TD+1)*(3*J1TD+J2TD+5)
     1/6+(J1TD+1)*(LAM2+J2TD-J2T)/2+(LAM1+J1TD-J1T)/2                   
      IDM(LAM,MU)=(LAM+1)*(MU+1)*(LAM+MU+2)/2                           
C     DIMENSION CHECKS (LSU,6-81)-START                                 
      IF(N1.GT.9)GO TO 130                                              
      IF(N2.GT.42)GO TO 140                                             
      IF(N3.GT.9030)GO TO 150                                           
      IDTEST=KR0MAX*IDM(LAM2,MU2)                                       
      IF(IDTEST.GT.KIMAX2.OR.IDTEST.GT.27090)GO TO 160                  
C     DIMENSION CHECKS (LSU,6-81)--STOP                                 
      LL1=LAM1+1                                                        
      MM1=MU1+1                                                         
      LL2=LAM2+1                                                        
      MM2=MU2+1                                                         
      LL3=LAM3+1                                                        
      MM3=MU3+1                                                         
      LM1=LAM1+MU1                                                      
      LM2=LAM2+MU2                                                      
      LLMM1=LL1+MM1                                                     
      LLMM2=LL2+MM2                                                     
      LLMM3=LL3+MM3                                                     
      JJTD=(IE+LAM3+2*MU3)/3+1                                          
      IP=(JJTD+JT-LL3)/2                                                
      NCC=NEC-1                                                         
      INC=1                                                             
      IQ3=0                                                             
      IP3=-1                                                            
      J3T=LAM3+IP3-IQ3                                                  
      DO 125 JJ3TD=1,JJTD                                               
C     DO 10 N=1,KIMAX2                                                  
      DO 10 N=1,IDTEST                                                  
   10 DWU3(N)=0.D0                                                      
      NCC=NCC+1                                                         
      IF(IP3.EQ.IP)INC=0                                                
      IF(INC.EQ.1)GO TO 15                                              
      IQ3=IQ3+1                                                         
      J3T=J3T-1                                                         
      NM=(LL3-IQ3)*IQ3*(LLMM3-IQ3)                                      
      GO TO 20                                                          
   15 IP3=IP3+1                                                         
      J3T=J3T+1                                                         
      NM=(MM3-IP3)*IP3*(LL3+IP3)                                        
   20 JJ2TDA=NCC-LM1                                                    
      IF(JJ2TDA.LT.0)JJ2TDA=0                                           
      JJ2TDA=JJ2TDA+1                                                   
      JJ2TDB=LM2                                                        
      IF(NCC.LT.JJ2TDB)JJ2TDB=NCC                                       
      JJ2TDB=JJ2TDB+1                                                   
      JJ2TDC=JJ2TDA                                                     
      IND=0                                                             
      IES=0                                                             
      DO 115 JJ2TD=JJ2TDA,JJ2TDB                                        
      J2TD=JJ2TD-1                                                      
      J1TD=NCC-J2TD                                                     
      IES=IES+1                                                         
      IIQ2A=J2TD-MU2                                                    
      IF(IIQ2A.LT.0)IIQ2A=0                                             
      IIQ2A=IIQ2A+1                                                     
      IIQ2B=J2TD                                                        
      IF(LAM2.LT.IIQ2B)IIQ2B=LAM2                                       
      IIQ2B=IIQ2B+1                                                     
      IIQ1A=J1TD-MU1                                                    
      IF(IIQ1A.LT.0)IIQ1A=0                                             
      IIQ1A=IIQ1A+1                                                     
      IIQ1B=J1TD                                                        
      IF(LAM1.LT.IIQ1B)IIQ1B=LAM1                                       
      IIQ1B=IIQ1B+1                                                     
      J2S=0                                                             
      DO 105 IIQ2=IIQ2A,IIQ2B                                           
      IQ2=IIQ2-1                                                        
      IP2=J2TD-IQ2                                                      
      J2T=LAM2+IP2-IQ2                                                  
      J23S=J2T+J3T                                                      
      J23D=J3T-J2T                                                      
      J23H=IABS(J23D)                                                   
      J1S=0                                                             
      DO 100 IIQ1=IIQ1A,IIQ1B                                           
      IQ1=IIQ1-1                                                        
      IP1=J1TD-IQ1                                                      
      J1T=LAM1+IP1-IQ1                                                  
      IF(J1T.LT.J23H)GO TO 100                                          
      IF(J1T.GT.J23S)GO TO 100                                          
      J1TS=J1T                                                          
      J2TS=J2T                                                          
      INDQ=IND*KR0MAX                                                   
      IND=IND+1                                                         
      J1S=J1S+1                                                         
      IF(JJ3TD.EQ.1)GO TO 90                                            
      JA=(J23S-J1T)/2                                                   
      JJA=JA+1                                                          
      JB=(J23D+J1T)/2                                                   
      JJB=JB+1                                                          
      JC=(J1T+J23S)/2+1                                                 
      JJC=JC+1                                                          
      JD=(J1T-J23D)/2+1                                                 
      JJD=JD-1                                                          
      IESP=J2TD-JJ2TDP                                                  
      DO 85 I=1,4                                                       
      IF(I.EQ.1)IESP=IESP+1                                             
      IF(I.EQ.3)IESP=IESP+1                                             
      IF(IESP.LT.1)GO TO 85                                             
      IF(IESP.GT.IESMAX)GO TO 85                                        
      GO TO (25,35,45,55),I                                             
   25 J2TP=J2T+1                                                        
      J1TP=J1T                                                          
      IF(J1TP.LT.IABS(J2TP-J3TP))GO TO 85                               
      IF(J1TP.GT.J2TP+J3TP)GO TO 85                                     
      M=IQ2                                                             
      IF(M.EQ.0)GO TO 85                                                
      N=LL2-M                                                           
      N=(LLMM2-M)*N                                                     
      J12TP=J2T+1                                                       
      IF(INC.EQ.1)GO TO 30                                              
      IAB=JJA                                                           
      ICD=JJC                                                           
      IPH=1                                                             
      GO TO 65                                                          
   30 IAB=JB                                                            
      ICD=JD                                                            
      IPH=-1                                                            
      GO TO 65                                                          
   35 J2TP=J2T-1                                                        
      J1TP=J1T                                                          
      IF(J1TP.LT.IABS(J2TP-J3TP))GO TO 85                               
      IF(J1TP.GT.J2TP+J3TP)GO TO 85                                     
      M=IP2                                                             
      IF(M.EQ.0)GO TO 85                                                
      N=MM2-M                                                           
      N=(LLMM2-N)*N                                                     
      J12TP=J2T                                                         
      IF(INC.EQ.1)GO TO 40                                              
      IAB=JJB                                                           
      ICD=JJD                                                           
      IPH=1                                                             
      GO TO 65                                                          
   40 IAB=JA                                                            
      ICD=JC                                                            
      IPH=1                                                             
      GO TO 65                                                          
   45 J2TP=J2T                                                          
      J1TP=J1T+1                                                        
      IF(J1TP.LT.IABS(J2TP-J3TP))GO TO 85                               
      IF(J1TP.GT.J2TP+J3TP)GO TO 85                                     
      M=IQ1                                                             
      IF(M.EQ.0)GO TO 85                                                
      N=LL1-M                                                           
      N=(LLMM1-M)*N                                                     
      J12TP=J1T+1                                                       
      IF(INC.EQ.1)GO TO 50                                              
      IAB=JJB                                                           
      ICD=JJC                                                           
      IPH=1                                                             
      GO TO 65                                                          
   50 IAB=JA                                                            
      ICD=JD                                                            
      IPH=1                                                             
      GO TO 65                                                          
   55 J2TP=J2T                                                          
      J1TP=J1T-1                                                        
      IF(J1TP.LT.IABS(J2TP-J3TP))GO TO 85                               
      IF(J1TP.GT.J2TP+J3TP)GO TO 85                                     
      M=IP1                                                             
      IF(M.EQ.0)GO TO 85                                                
      N=MM1-M                                                           
      N=(LLMM1-N)*N                                                     
      J12TP=J1T                                                         
      IF(INC.EQ.1)GO TO 60                                              
      IAB=JJA                                                           
      ICD=JJD                                                           
      IPH=-1                                                            
      GO TO 65                                                          
   60 IAB=JB                                                            
      ICD=JC                                                            
      IPH=1                                                             
   65 IF(J12TP.GT.0)GO TO 70                                            
      IF(INC.EQ.1)IAB=1                                                 
      IF(INC.EQ.0)ICD=1                                                 
      DC=1.D0                                                           
      GO TO 75                                                          
   70 DC=DFLOAT(J12TP*(J12TP+1))                                        
   75 J2SP=(J2TMAP(IESP)-J2TP+2)/2                                      
      INDP=(INDMAP(IESP+(J2SP-1)*N2)-J1TP)/2                            
      DC=DSQRT(DFLOAT(IAB*ICD*M*N)/(DFLOAT(NM)*DC))                     
      IF(IPH.LT.0)DC=-DC                                                
      INDPQ=(INDP-1)*KR0MAX                                             
      DO 80 KR0=1,KR0MAX                                                
      KI=KR0+INDQ                                                       
      KIP=KR0+INDPQ                                                     
   80 DWU3(KI)=DWU3(KI)+DC*DWU3P(KIP)                                   
   85 CONTINUE                                                          
      GO TO 100                                                         
   90 INDPQ=(INDEX(J1TD,LAM1,J1T,J2TD,LAM2,J2T)-1)*KR0MAX               
      DO 95 KR0=1,KR0MAX                                                
      KI=KR0+INDQ                                                       
      KIP=KR0+INDPQ                                                     
   95 DWU3(KI)=DEWU3(KIP)                                               
  100 CONTINUE                                                          
      IF(J1S.EQ.0)GO TO 105                                             
      IESJ2S=IES+J2S*N2                                                 
      J2S=J2S+1                                                         
      J1SMAX(IESJ2S)=J1S                                                
      J1TMAX(IESJ2S)=J1TS+2*(J1S-1)                                     
      INDMAT(IESJ2S)=2*IND+J1TS                                         
  105 CONTINUE                                                          
      IF(J2S.NE.0)GO TO 110                                             
      IES=IES-1                                                         
      IF(IES.EQ.0)JJ2TDC=JJ2TDC+1                                       
      GO TO 115                                                         
  110 J2SMAX(IES)=J2S                                                   
      J2TMAX(IES)=J2TS+2*(J2S-1)                                        
  115 CONTINUE                                                          
      IESMAX=IES                                                        
      IF(JJ3TD.EQ.JJTD)GO TO 125                                        
      J3TP=J3T                                                          
      JJ2TDP=JJ2TDC                                                     
      IND=0                                                             
      DO 120 IES=1,IESMAX                                               
      J2TMAP(IES)=J2TMAX(IES)                                           
      J2SB=J2SMAX(IES)                                                  
      J2SQ=-N2                                                          
      DO 120 J2S=1,J2SB                                                 
      J2SQ=J2SQ+N2                                                      
      IESJ2S=IES+J2SQ                                                   
      INDMAP(IESJ2S)=INDMAT(IESJ2S)                                     
      J1SB=J1SMAX(IESJ2S)                                               
      DO 120 J1S=1,J1SB                                                 
      INDQ=IND*KR0MAX                                                   
      IND=IND+1                                                         
      DO 120 KR0=1,KR0MAX                                               
      KI=KR0+INDQ                                                       
  120 DWU3P(KI)=DWU3(KI)                                                
  125 CONTINUE                                                          
      INDMAX=IND                                                        
      IE2MAX=-(LAM2+2*MU2)+3*(JJ2TDC-1)+3*(IES-1)                       
      RETURN                                                            
C     DIMENSION CHECKS (LSU,6-81)-START                                 
  130 WRITE(6,135)N1                                                    
  135 FORMAT(35H *****XWU3 DIMENSION OVERFLOW:  N1=,I10)                
      GO TO 170                                                         
  140 WRITE(6,145)N2                                                    
  145 FORMAT(35H *****XWU3 DIMENSION OVERFLOW:  N2=,I10)                
      GO TO 170                                                         
  150 WRITE(6,155)N3                                                    
  155 FORMAT(35H *****XWU3 DIMENSION OVERFLOW:  N3=,I10)                
      GO TO 170                                                         
  160 WRITE(6,165)IDTEST                                                
  165 FORMAT(39H *****XWU3 DIMENSION OVERFLOW:  IDTEST=,I10)            
  170 STOP                                                              
C     DIMENSION CHECKS (LSU,6-81)--STOP                                 
      END                                                               
