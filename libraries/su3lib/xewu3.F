      SUBROUTINE XEWU3(LAM1X,MU1X,LAM2X,MU2X,LAM3X,MU3X,I3,NEC,KR0MAX,  
     1INDMAX,DEWU3,J1TA,J2TA,IEA,N1,N2,KIMAX1)                          
C     ------------------------------------------------------------------
C     EXTREMAL WIGNER COEFFICIENTS FOR U3 (X PREFIX FOR 6-81 VERSION)   
C     ------------------------------------------------------------------
C     UPDATE/MOD: (LSU,05-80)  J.P.DRAAYER        LOG BINOMIALS         
C                 (LSU,06-81)  J.P.DRAAYER        INDEXING DEWU3        
C                 (LSU,03-83)  J.P.DRAAYER        SPACE SAVING MEASURE  
C                 (LSU,02-87)  J.P.DRAAYER        OVERFLOW CORRECTION   
C                 (LSU,10-89)  J.P.DRAAYER        ZERO OUT RELOCATED    
C                 (UoT,04-97)  C.BAHRI            NX=82
C                                                                       
C     REFERENCES--J.P.DRAAYER AND Y.AKIYAMA, J.MATH.PHYS.14(1973)1904   
C                 K.T.HECHT, NUCL.PHYS.62(1965)1                        
C     PARAMETERS--(I3) : (1)=GHW, (0)=GLW                               
C       EXTERNAL--N1=MAX(KR0MAX)                                        
C                 N2=MAX(INDMAX)=NX*(NX+1)*(NX+2)/6, NX=MAX(LAM2+MU2+1) 
C*                KIMAX1=MAX(KR0MAX*INDMAX)                             
C       INTERNAL--X1=ABS(N1*NX)                                         
C                 X2=ABS(NX)                                            
C     EXTENSIONS--CHANGE EXTERNAL PARAMETERS IN CALL                    
C                 ADJUST INTERNAL PARAMETERS BELOW                      
C*    DIMENSIONS--DEWU3(N1*N2->KIMAX1),J1TA(N2),J2TA(N2),IEA(N2),       
C                 DEWU3P(X1),DZ(X2),J1TAP(X2),IAB(X2),ICD(X2)           
CB      COMMENTS--ASSUME MAX N1=9,NX=42,N2=13244                        
CB                       SET X1=378,X2=42                               
C-                           N1=9,NX=82,N2=95284
C-                           X1=738,X2=82
C                 DZ ARRAY ADDED FOR CORRECTING THE OVERFLOW PROBLEM    
C     ------------------------------------------------------------------
      IMPLICIT REAL*8(D)                                                
      COMMON/BKDB/DBINO(8385),DBINV(8385),DTWOS(-128:128)               
      COMMON/BKDF/DLOGF(0:2000)                                         
      DIMENSION DEWU3(1),J1TA(1),J2TA(1),IEA(1),                        
     1          DEWU3P(738),DZ(82),J1TAP(82),IAB(82),ICD(82)            
CB   1          DEWU3P(378),DZ(42),J1TAP(42),IAB(42),ICD(42)            
      INDEX(J1TD,LAM1,J1T,J2TD,LAM2,J2T)=1+J2TD*(J2TD+1)*(3*J1TD+J2TD+5)
     1/6+(J1TD+1)*(LAM2+J2TD-J2T)/2+(LAM1+J1TD-J1T)/2                   
C     DIMENSION CHECKS (LSU,6-81)-START                                 
      IF(N1.GT.9)GO TO 200                                              
      NX=LAM2X+MU2X+1                                                   
CB    IF(NX.GT.42)GO TO 210                                             
      IF(NX.GT.82)GO TO 210                                             
C     DIMENSION CHECKS (LSU,6-81)-START                                 
      KR0MAX=MULTU3(LAM1X,MU1X,LAM2X,MU2X,LAM3X,MU3X)                   
      IF(KR0MAX.EQ.0)RETURN                                             
      IF(I3.EQ.1)GO TO 10                                               
      LAM1=MU1X                                                         
      LAM2=MU2X                                                         
      LAM3=MU3X                                                         
      MU1=LAM1X                                                         
      MU2=LAM2X                                                         
      MU3=LAM3X                                                         
      GO TO 15                                                          
   10 LAM1=LAM1X                                                        
      LAM2=LAM2X                                                        
      LAM3=LAM3X                                                        
      MU1=MU1X                                                          
      MU2=MU2X                                                          
      MU3=MU3X                                                          
   15 NEC=(LAM1+LAM2-LAM3+2*(MU1+MU2-MU3))/3                            
      IAH=(LAM2+LAM3-LAM1-NEC)/2                                        
      IBH=(LAM3+LAM1-LAM2+NEC+2)/2                                      
      ICH=(LAM1+LAM2-LAM3-NEC)/2                                        
      IDH=(LAM1+LAM2+LAM3-NEC+2)/2                                      
      DO 20 I=1,NEC                                                     
      IAB(I)=(IAH+I)*(IBH-I)                                            
   20 ICD(I)=(ICH+I)*(IDH+I)                                            
      NCDMAX=MULTHY(LAM1,MU1,LAM2,MU2,LAM3,MU3)                         
      NEC=NEC-NCDMAX                                                    
      LAM2=LAM2-NCDMAX                                                  
      MU2=MU2-NCDMAX                                                    
      NCDMIN=1                                                          
   25 IF(NCDMIN.EQ.NCDMAX)GO TO 30                                      
      IF(MULTU3(LAM1,MU1,LAM2+1,MU2+1,LAM3,MU3).GT.0)GO TO 30           
      NEC=NEC+1                                                         
      LAM2=LAM2+1                                                       
      MU2=MU2+1                                                         
      NCDMIN=NCDMIN+1                                                   
      GO TO 25                                                          
C     DIMENSION MODIFICATION (LSU,6-81)-START                           
   30 NNCMAX=NEC+NCDMAX-NCDMIN+2                                        
      KITEST=KR0MAX*(NNCMAX)*(NNCMAX+1)*(NNCMAX+2)/6                    
      IF(KITEST.GT.KIMAX1)GO TO 220                                     
C     DIMENSION MODIFICATION (LSU,6-81)--STOP                           
      DO I=1,KITEST                                                     
         DEWU3(I)=0.D0                                                  
      ENDDO                                                             
      LL1=LAM1+1                                                        
      MM1=MU1+1                                                         
      LL2=LAM2+1                                                        
      MM2=MU2+1                                                         
      IA1=2*LAM1+4*MU1                                                  
      IB1=4*LAM1+2*MU1                                                  
      IC1=IB1-IA1                                                       
      IA2=2*LAM2+4*MU2                                                  
      IB2=4*LAM2+2*MU2                                                  
      IC2=IB2-IA2                                                       
      IS1=LL1+MM1                                                       
      IS2=LL2+MM2                                                       
      ISS=MM1+LAM2+MU2-NEC                                              
      IE3=-(LAM3+2*MU3)                                                 
      IEH=-(LAM2+2*MU2+3)                                               
      KR0CNT=0                                                          
      DO 135 NCD=NCDMIN,NCDMAX                                          
      NEC=NEC+1                                                         
      LAM2=LAM2+1                                                       
      MU2=MU2+1                                                         
      NNC=NEC+1                                                         
      INDMAX=NNC*(NNC+1)*(NNC+2)/6                                      
      IA2=IA2+6                                                         
      IB2=IB2+6                                                         
      IS2=IS2+2                                                         
      ISS=ISS+1                                                         
      IEH=IEH-3                                                         
      LL2=LAM2+1                                                        
      MM2=MU2+1                                                         
      LN1=LAM1+NEC                                                      
      LN2=LAM2+NEC                                                      
      INN=NEC*NNC/2                                                     
      IF(NCD.EQ.NCDMIN)GO TO 40                                         
      DO 35 I=1,KITEST                                                  
   35 DEWU3(I)=0.D0                                                     
   40 DO 45 IND=1,INDMAX                                                
      IEA(IND)=-1000                                                    
      J2TA(IND)=-1000                                                   
   45 J1TA(IND)=-1000                                                   
      IE2=IEH                                                           
      I=1000                                                            
      DO 55 IIE=1,NNC                                                   
      IE2=IE2+3                                                         
      IE1=IE3-IE2                                                       
      J2TD=IIE-1                                                        
      J1TD=NNC-IIE                                                      
      JJ2TA=IA2-IE2                                                     
      JJ2TB=IB2+IE2                                                     
      IF(JJ2TB.LT.JJ2TA)JJ2TA=JJ2TB                                     
      JJ2TA=JJ2TA/3+1                                                   
      JJ2TB=JJ2TA-IABS(IC2-IE2)/3                                       
      JJ1TA=IA1-IE1                                                     
      JJ1TB=IB1+IE1                                                     
      IF(JJ1TB.LT.JJ1TA)JJ1TA=JJ1TB                                     
      JJ1TA=JJ1TA/3+1                                                   
      JJ1TB=JJ1TA-IABS(IC1-IE1)/3                                       
      J=0                                                               
      DO 50 JJ2T=1,JJ2TB,2                                              
      J2T=JJ2TA-JJ2T                                                    
      L=IABS(J2T-LAM3)                                                  
      M=J2T+LAM3                                                        
      DO 50 JJ1T=1,JJ1TB,2                                              
      J1T=JJ1TA-JJ1T                                                    
      IF(J1T.LT.L)GO TO 50                                              
      IF(J1T.GT.M)GO TO 50                                              
      IND=INDEX(J1TD,LAM1,J1T,J2TD,LAM2,J2T)                            
      IEA(IND)=IE2                                                      
      J2TA(IND)=J2T                                                     
      J1TA(IND)=J1T                                                     
      J=J+1                                                             
   50 CONTINUE                                                          
      IF(J.LT.I)I=J                                                     
   55 CONTINUE                                                          
      IF(I.EQ.0)GO TO 135                                               
      IF(KR0CNT.EQ.0)GO TO 80                                           
C     GENERATE <(LAM1,MU1)????;(LAM2,MU2)HIGH::KR0(LAM3,MU3)HIGH>       
C     FROM <(LAM1,MU1)????;(LAM2-1,MU2-1)HIGH::KR0(LAM3,MU3)HIGH>       
      INDQ=-KR0MAX                                                      
      DO 75 IND=1,NNC                                                   
      INDQ=INDQ+KR0MAX                                                  
      J1T=J1TA(IND)                                                     
      IF(J1T.LT.0)GO TO 75                                              
      IQ1=(LN1-J1T)/2                                                   
      IF(IQ1.LE.0)GO TO 65                                              
      J1TP=J1T+1                                                        
      INDP=(LN1-J1TP-1)/2+1                                             
      IF(J1TAP(INDP).LT.0)GO TO 65                                      
      I=IAB(IQ1)*IQ1*(LL1-IQ1)*(IS1-IQ1)                                
      DC=-DSQRT(DFLOAT(I)/DFLOAT((J1T+2)*J1TP))                         
      INDPQ=(INDP-1)*KR0MAX                                             
      DO 60 KR0=1,KR0CNT                                                
      KI=KR0+INDQ                                                       
      KIP=KR0+INDPQ                                                     
   60 DEWU3(KI)=DC*DEWU3P(KIP)                                          
   65 IP1=NEC-IQ1                                                       
      IF(IP1.LE.0)GO TO 75                                              
      J1TP=J1T-1                                                        
      INDP=(LN1-J1TP-1)/2+1                                             
      IF(J1TAP(INDP).LT.0)GO TO 75                                      
      I=ICD(NNC-IND)*IP1*(MM1-IP1)*(LL1+IP1)                            
      DC=DSQRT(DFLOAT(I)/DFLOAT((J1TP+2)*J1T))                          
      INDPQ=(INDP-1)*KR0MAX                                             
      DO 70 KR0=1,KR0CNT                                                
      KI=KR0+INDQ                                                       
      KIP=KR0+INDPQ                                                     
   70 DEWU3(KI)=DEWU3(KI)+DC*DEWU3P(KIP)                                
   75 CONTINUE                                                          
      INC=0                                                             
      IF(KR0CNT.EQ.KR0MAX)GO TO 125                                     
C     EVALUATE <(LAM1,MU1)HIGH;(LAM2,MU2)????::KR0(LAM3,MU3)HIGH>       
C     WITH (LAM2,MU2) A MINIMUM FOR KR0=KR0CNT                          
   80 KR0CNT=KR0CNT+1                                                   
      I=0                                                               
      IND=INDEX(0,LAM1,LAM1,NEC,LAM2,LN2)-1                             
      INDQ=-KR0MAX                                                      
      DO 85 IIQ2=1,NNC                                                  
      INDQ=INDQ+KR0MAX                                                  
      IND=IND+1                                                         
      KI=KR0CNT+INDQ                                                    
      DEWU3P(KI)=0.D0                                                   
      IF(J1TA(IND).LT.0)GO TO 85                                        
      I=I+1                                                             
      IIQ2B=IIQ2                                                        
   85 CONTINUE                                                          
C                                                                       
C     *****MODIFIED TO AVOID OVERFLOW (LSU,5-80)-START                  
C      ****DLOGB CHANGED TO SAVE SPACE (LSU,3-83)****                   
C      ****FURTHER OVERFLOW CORRECTION (LSU,2-87)****                   
C                                                                       
      IIQ2A=IIQ2B-I+1                                                   
      IQ2B=IIQ2B-1                                                      
      INDQ=(IIQ2A-2)*KR0MAX                                             
      IZ=0                                                              
      DO 115 IIQ2=IIQ2A,IIQ2B                                           
      IZ=IZ+1                                                           
      DZ(IZ)=1.D0                                                       
      INDQ=INDQ+KR0MAX                                                  
      L=LL2-IIQ2                                                        
C                                                                       
C --> START NUMERATOR PRODUCT LOOP                                      
C                                                                       
      IX=L-IIQ2+NNC+1                                                   
      IF(IX.EQ.0)GO TO 120                                              
      IY=IABS(IX)                                                       
      IN=IX/IY                                                          
      DN=DLOG(DFLOAT(IY))                                               
      IF(IIQ2A.EQ.IIQ2B)GO TO 95                                        
      DO 90 I=IIQ2A,IQ2B                                                
      J=NNC-I                                                           
      IF(I.LT.IIQ2)THEN                                                 
      K=IAB(J)*ICD(J)*(IS2-I)                                           
      ELSE                                                              
      K=MM2-J                                                           
      ENDIF                                                             
      IF(K.EQ.0)GO TO 120                                               
      IF(K.LT.0)IN=-IN                                                  
   90 DN=DLOG(DFLOAT(IABS(K)))+DN                                       
   95 DN=DN+DLOG(DBINO(INN+IIQ2))                                       
C                                                                       
C --> END NUMERATOR PRODUCT LOOP & START DENOMINATOR PRODUCT LOOP       
C                                                                       
      ID=1                                                              
      DD=0.D0                                                           
      DO 100 I=1,NNC                                                    
      IX=I+L                                                            
      IF(IX.LT.0)ID=-ID                                                 
  100 DD=DLOG(DFLOAT(I+L))+DD                                           
C                                                                       
C --> END DENOMINATOR PRODUCT LOOP & START INNER PRODUCT/SUM LOOP       
C                                                                       
      IP2=NNC-IIQ2                                                      
C                                                                       
C     MULTIPLY BY SMALL NUMBER --> DEXP(-172) LIMIT FOR IBM SYSTEMS     
C                                                                       
      DZ(IZ)=DEXP(-DMIN1(DLOGF(2*IP2),172.D0))                          
      IIP2=IP2+1                                                        
      M=IP2*IIP2/2                                                      
      DS=0.D0                                                           
      DO 110 I=1,IIP2                                                   
      DC=DZ(IZ)*DBINO(I+M)                                              
      IF(IIP2.EQ.1)GO TO 110                                            
      DO 105 J=1,IP2                                                    
      IF(J.LT.I)THEN                                                    
         K=(J+L)*(ISS+J)                                                
      ELSE                                                              
         K=IAB(J)                                                       
      ENDIF                                                             
  105 DC=DFLOAT(K)*DC                                                   
  110 DS=DS+DC                                                          
C                                                                       
C --> END INNER PRODUCT/SUM LOOP & ASSIGN unnormalized DEWU3P VALUE     
C                                                                       
      IF(2*(IP2/2).NE.IP2)DS=-DS                                        
      KI=KR0CNT+INDQ                                                    
  115 DEWU3P(KI)=DFLOAT(IN*ID)*DS*DEXP((DN-DD)/2.D0)                    
C                                                                       
C --> START renormalization PROCEDURE                                   
C                                                                       
      DMIN=1.D0                                                         
      IZ=0                                                              
      DO 117 IIQ2=IIQ2A,IIQ2B                                           
      IZ=IZ+1                                                           
  117 DMIN=DMIN1(DMIN,DZ(IZ))                                           
      INDQ=(IIQ2A-2)*KR0MAX                                             
      IZ=0                                                              
      DO 118 IIQ2=IIQ2A,IIQ2B                                           
      IZ=IZ+1                                                           
      INDQ=INDQ+KR0MAX                                                  
      KI=KR0CNT+INDQ                                                    
  118 DEWU3P(KI)=(DMIN/DZ(IZ))*DEWU3P(KI)                               
C                                                                       
C     *****MODIFIED TO AVOID OVERFLOW (LSU,5-80)--STOP                  
C                                                                       
      KR0A=KR0CNT                                                       
      KR0B=KR0CNT                                                       
C     GENERATE <(LAM1,MU1)????;(LAM2,MU2)????::KR0(LAM3,MU3)HIGH>       
C     FROM <(LAM1,MU1)HIGH;(LAM2,MU2)????::KR0(LAM3,MU3)HIGH>           
      CALL XEWU3S(1,LAM1,MU1,LAM2,MU2,NEC,NNC,KR0A,KR0B,DEWU3P,J1TA,    
     1IAB,ICD,INDMAX,DEWU3,KR0MAX)                                      
      INC=1                                                             
      GO TO 125                                                         
  120 KR0CNT=KR0CNT-1                                                   
      WRITE(1,195)LAM1X,MU1X,LAM2X,MU2X,LAM3X,MU3X,KR0MAX,KR0CNT        
  125 IF(KR0CNT.EQ.0)GO TO 135                                          
      INDQ=-KR0MAX                                                      
      DO 130 IND=1,NNC                                                  
      INDQ=INDQ+KR0MAX                                                  
      J1TAP(IND)=J1TA(IND)                                              
      DO 130 KR0=1,KR0CNT                                               
      KI=KR0+INDQ                                                       
  130 DEWU3P(KI)=DEWU3(KI)                                              
  135 CONTINUE                                                          
      IF(KR0CNT.EQ.0)RETURN                                             
      KR0A=1                                                            
      KR0B=KR0CNT-INC                                                   
      IF(KR0B.EQ.0)GO TO 140                                            
C     GENERATE <(LAM1,MU1)????;(LAM2,MU2)????::KR0(LAM3,MU3)HIGH>       
C     FROM <(LAM1,MU1)????;(LAM2,MU2)HIGH::KR0(LAM3,MU3)HIGH>           
      CALL XEWU3S(0,LAM1,MU1,LAM2,MU2,NEC,NNC,KR0A,KR0B,DEWU3P,J1TA,    
     1IAB,ICD,INDMAX,DEWU3,KR0MAX)                                      
  140 CONTINUE                                                          
C     RENORMALIZE VIA LARGEST ELEMENT TO AVOID OVERFLOW (LSU,5-80)-START
      DO 142 KR0=1,KR0CNT                                               
      DC=1.D0                                                           
      INDQ=-KR0MAX                                                      
      DO 141 IND=1,INDMAX                                               
      INDQ=INDQ+KR0MAX                                                  
      KI=KR0+INDQ                                                       
  141 DC=DMAX1(DC,DABS(DEWU3(KI)))                                      
      INDQ=-KR0MAX                                                      
      DO 142 IND=1,INDMAX                                               
      INDQ=INDQ+KR0MAX                                                  
      KI=KR0+INDQ                                                       
  142 DEWU3(KI)=DEWU3(KI)/DC                                            
C     RENORMALIZE VIA LARGEST ELEMENT TO AVOID OVERFLOW (LSU,5-80)--STOP
C     ORTHONORMALIZATION OF SOLUTIONS                                   
      DO 165 KR0=1,KR0CNT                                               
      KR0PA=1                                                           
      IF(INC.EQ.1.AND.KR0.EQ.KR0CNT)KR0PA=KR0CNT                        
      DO 165 KR0P=KR0PA,KR0                                             
      DN=0.D0                                                           
      INDQ=-KR0MAX                                                      
      DO 145 IND=1,INDMAX                                               
      INDQ=INDQ+KR0MAX                                                  
      KI=KR0+INDQ                                                       
      KIP=KR0P+INDQ                                                     
  145 DN=DN+DEWU3(KI)*DEWU3(KIP)                                        
      IF(KR0P.EQ.KR0)GO TO 155                                          
      IF(DABS(DN).LT.1.D-12)GO TO 165                                   
      INDQ=-KR0MAX                                                      
      DO 150 IND=1,INDMAX                                               
      INDQ=INDQ+KR0MAX                                                  
      KI=KR0+INDQ                                                       
      KIP=KR0P+INDQ                                                     
  150 DEWU3(KI)=DEWU3(KI)-DN*DEWU3(KIP)                                 
      GO TO 165                                                         
  155 DN=1.D0/DSQRT(DN)                                                 
      INDQ=-KR0MAX                                                      
      DO 160 IND=1,INDMAX                                               
      INDQ=INDQ+KR0MAX                                                  
      KI=KR0+INDQ                                                       
  160 DEWU3(KI)=DN*DEWU3(KI)                                            
  165 CONTINUE                                                          
C     SET PHASE CONVENTION (K.T.HECHT, NUCL.PHYS.62(1965)1)             
      IE2=IE3+(LAM1+2*MU1)                                              
      IPH=2*(LAM1+LAM2-LAM3+MU1+MU2-MU3+KR0MAX)                         
      INDQ=-KR0MAX                                                      
      DO 180 IND=1,INDMAX                                               
      INDQ=INDQ+KR0MAX                                                  
      IF(IEA(IND).NE.IE2)GO TO 180                                      
      I=IPH+J1TA(IND)+J2TA(IND)-LAM3                                    
      DO 175 KR0=1,KR0MAX                                               
      I=I-2                                                             
      KI=KR0+INDQ                                                       
      J=I                                                               
      IF(DEWU3(KI).LT.0.D0)J=J-2                                        
      IF(4*(J/4).EQ.J)GO TO 175                                         
      INDPQ=-KR0MAX                                                     
      DO 170 INDP=1,INDMAX                                              
      INDPQ=INDPQ+KR0MAX                                                
      KIP=KR0+INDPQ                                                     
  170 DEWU3(KIP)=-DEWU3(KIP)                                            
  175 CONTINUE                                                          
      GO TO 185                                                         
  180 CONTINUE                                                          
  185 IF(I3.EQ.1)RETURN                                                 
      INDQ=-KR0MAX                                                      
      DO 190 IND=1,INDMAX                                               
      INDQ=INDQ+KR0MAX                                                  
      IEA(IND)=-IEA(IND)                                                
      I=IPH+J1TA(IND)+J2TA(IND)-LAM3                                    
      DO 190 KR0=1,KR0MAX                                               
      I=I-2                                                             
      KI=KR0+INDQ                                                       
      IF(4*(I/4).NE.I)DEWU3(KI)=-DEWU3(KI)                              
  190 CONTINUE                                                          
  195 FORMAT(28H *****U3 COUPLING ERROR*****,3X,3(4X,2I3),3X,           
     115HRO(ABSOLUTE) = ,I2,3X,15HRO(RELATIVE) = ,I2,4X,                
     226H*****REPORT TO AUTHOR*****)                                    
      RETURN                                                            
C     DIMENSION CHECKS (LSU,6-81)-START                                 
  200 WRITE(6,205)N1                                                    
  205 FORMAT(36H ***** XEWU3 DIMENSION OVERFLOW: N1=,I10)               
      GO TO 230                                                         
  210 WRITE(6,215)NX                                                    
  215 FORMAT(36H ***** XEWU3 DIMENSION OVERFLOW: NX=,I10)               
      GO TO 230                                                         
  220 WRITE(6,225)KITEST,KIMAX1                                         
  225 FORMAT(40H ***** XEWU3 DIMENSION OVERFLOW: KITEST=,I10,5X,        
     17HKIMAX1=,I10)                                                    
  230 STOP                                                              
C     DIMENSION CHECKS (LSU,6-81)--STOP                                 
      END                                                               
