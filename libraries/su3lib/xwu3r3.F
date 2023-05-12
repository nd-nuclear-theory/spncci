      SUBROUTINE XWU3R3(I1,J1,LAM1,MU1,KA1MAX,L1,DON1,N1,I2,J2,LAM2,MU2,
     1KA2MAX,L2,DON2,N2,I3,J3,LAM3,MU3,KA3MAX,L3,DON3,N3,KR0MAX,INDMAX, 
     2DEWU3,J1TA,J2TA,IEA,DWU3R3,N0)                                    
C     ------------------------------------------------------------------
C     WIGNER COEFFICIENTS FOR U3R3 REDUCTION (X PREFIX FOR 6-81 VERSION)
C     ------------------------------------------------------------------
C     UPDATE/MOD: (LSU,06-81)  J.P.DRAAYER        INDEXING DEWU3        
C                                                                       
C     REFERENCES--J.P.DRAAYER AND Y.AKIYAMA, J.MATH.PHYS.14(1973)1904   
C                 J.P.DRAAYER, NUCL.PHYS.A129(1969)647-665              
C     PARAMETERS--SEE ALSO XEWU3 AND DTU3R3 AND/OR CONMAT               
C       EXTERNAL--N0=MAX(KR0MAX)                                        
C                 N1=MAX(KA1MAX)                                        
C                 N2=MAX(KA2MAX)                                        
C                 N3=MAX(KA3MAX)                                        
C                 NA=MAX(INDMAX)=NX*(NX+1)*(NX+2)/6, NX=MAX(LAM2+MU2+1) 
C       INTERNAL--X1=ABS(N1)                                            
C                 X2=ABS(N2)                                            
C                 X3=ABS(N1*N2*N3)                                      
C                 X4=ABS(N1*N2)                                         
C     EXTENSIONS--CHANGE EXTERNAL PARAMETERS IN CALL STATEMENT AND      
C                 ADJUST INTERNAL PARAMETERS IN STATEMENT 28            
C*    DIMENSIONS--DON1(N1*N1),DON2(N2*N2),DON3(N3*N3),DEWU3(N0*NA),     
C                 J1TA(NA),J2TA(NA),IEA(NA),DWU3R3(N0*N1*N2*N3),        
C                 DT1A(X1),DT2A(X2),DS1A(X3),DS2A(X4)                   
C       COMMENTS--ASSUME MAX N0=9,N1=9,N2=9,N3=9                        
C                        SET X1=9,X2=9,X3=729,X4=81                     
C     ------------------------------------------------------------------
      IMPLICIT REAL*8(D),INTEGER(X)                                     
C     IMPLICIT INTEGER(X)                                               
      DIMENSION DON1(1),DON2(1),DON3(1),DEWU3(1),                       
     1          J1TA(1),J2TA(1),IEA(1),DWU3R3(1),                       
     2          DT1A(9),DT2A(9),DS1A(729),DS2A(81)                      
      KSTART(LAM,MU,L)=MOD(LAM,2)+2*(MAX0(0,(L-MU)/2)+MOD(LAM+1,2)*     
     1MOD(IABS(L-MU),2))                                                
C     DIMENSION CHECKS (LSU,6-81)-START                                 
C     icount = 0
      IF(N1.GT.9)GO TO 145                                              
      IF(N2.GT.9)GO TO 155                                              
      IF(N3.GT.9)GO TO 165                                              
C     DIMENSION CHECKS (LSU,6-81)--STOP                                 
      N01=N0*N1                                                         
      N012=N01*N2                                                       
      N0123=N012*N3                                                     
      N12=N1*N2                                                         
      N123=N12*N3                                                       
      DO 10 N=1,N0123                                                   
   10 DWU3R3(N)=0.D0                                                    
      IF(KR0MAX.EQ.0)RETURN                                             
      IF(I3.EQ.1)GO TO 15                                               
      K3S=KSTART(MU3,LAM3,L3)-2                                         
      IE3=2*LAM3+MU3                                                    
      J3T=MU3                                                           
      M3T=-MU3                                                          
      GO TO 20                                                          
   15 K3S=KSTART(LAM3,MU3,L3)-2                                         
      IE3=-(LAM3+2*MU3)                                                 
      J3T=LAM3                                                          
      M3T=LAM3                                                          
   20 IF(I3.NE.J3)M3T=-M3T                                              
      I2X=LAM2+MU2+L2                                                   
      IF(I2.EQ.1)GO TO 25                                               
      K2S=KSTART(MU2,LAM2,L2)-2                                         
      IF(I2.EQ.J2)I2X=I2X-MU2                                           
      I2Y=6*(LAM2+MU2)                                                  
      GO TO 30                                                          
   25 K2S=KSTART(LAM2,MU2,L2)-2                                         
      IF(I2.EQ.J2)I2X=I2X-LAM2                                          
      I2Y=0                                                             
   30 I1X=LAM1+MU1+L1                                                   
      IF(I1.EQ.1)GO TO 35                                               
      K1S=KSTART(MU1,LAM1,L1)-2                                         
      IF(I1.EQ.J1)I1X=I1X-MU1                                           
      I1Y=6*(LAM1+MU1)                                                  
      GO TO 40                                                          
   35 K1S=KSTART(LAM1,MU1,L1)-2                                         
      IF(I1.EQ.J1)I1X=I1X-LAM1                                          
      I1Y=0                                                             
   40 IS=I2X                                                            
      I2X=1                                                             
      IF(2*(IS/2).EQ.IS)I2X=0                                           
      IS=I1X                                                            
      I1X=1                                                             
      IF(2*(IS/2).EQ.IS)I1X=0                                           
      I2Y=2*LAM2+4*MU2+6*(L2+K2S)-I2Y                                   
      I1Y=2*LAM1+4*MU1+6*(L1+K1S)-I1Y                                   
      INDP=-KR0MAX                                                      
      DO 130 IND=1,INDMAX                                               
      INDP=INDP+KR0MAX                                                  
      IF(J1TA(IND).LT.0)GO TO 130                                       
      IE2=IEA(IND)                                                      
      IE1=IE3-IE2                                                       
      I2Z=(I2Y-IE2)/3                                                   
      I1Z=(I1Y-IE1)/3                                                   
      J2T=J2TA(IND)                                                     
      J1T=J1TA(IND)                                                     
      MM2TA=J1T+M3T                                                     
      IF(J2T.LT.MM2TA)MM2TA=J2T                                         
      MM2TA=MM2TA+1                                                     
      MM2TB=J1T-M3T                                                     
      IF(J2T.LT.MM2TB)MM2TB=J2T                                         
      MM2TB=MM2TA+MM2TB                                                 
      DO 45 N=1,N123                                                    
   45 DS1A(N)=0.D0                                                      
      K3=K3S                                                            
      KA3Q=-N12                                                         
      DO 120 KA3=1,KA3MAX                                               
      KA3Q=KA3Q+N12                                                     
      K3=K3+2                                                           
      MM2A=L1+K3                                                        
      IF(L2.LT.MM2A)MM2A=L2                                             
      MM2A=MM2A+1                                                       
      MM2B=L1-K3                                                        
      IF(L2.LT.MM2B)MM2B=L2                                             
      MM2B=MM2A+MM2B                                                    
      DO 120 MM2=1,MM2B                                                 
      M2=MM2A-MM2                                                       
      IF(IABS(M2).GT.J2T)GO TO 120                                      
      I2S=J2T+M2                                                        
      IF(2*(I2S/2).NE.I2S)GO TO 120                                     
      M1=K3-M2                                                          
      IF(IABS(M1).GT.J1T)GO TO 120                                      
      I1S=J1T+M1                                                        
      IF(2*(I1S/2).NE.I1S)GO TO 120                                     
      DO 50 N=1,N12                                                     
   50 DS2A(N)=0.D0                                                      
      DO 110 MM2T=1,MM2TB,2                                             
      M2T=MM2TA-MM2T                                                    
      IF(M2T.NE.0)GO TO 55                                              
      IF(4*(I2S/4).NE.I2S)GO TO 110                                     
   55 IF(M2.NE.0)GO TO 60                                               
      IS=I2Z+M2T                                                        
      IF(4*(IS/4).NE.IS)GO TO 110                                       
      IS=(J2T+M2T)/2                                                    
      IF(2*(IS/2).NE.IS)GO TO 110                                       
   60 M1T=M3T-M2T                                                       
      IF(M1T.NE.0)GO TO 65                                              
      IF(4*(I1S/4).NE.I1S)GO TO 110                                     
   65 IF(M1.NE.0)GO TO 70                                               
      IS=I1Z+M1T                                                        
      IF(4*(IS/4).NE.IS)GO TO 110                                       
      IS=(J1T+M1T)/2                                                    
      IF(2*(IS/2).NE.IS)GO TO 110                                       
   70 K1=K1S                                                            
      KA1Q=-N1                                                          
      DO 85 KA1=1,KA1MAX                                                
      KA1Q=KA1Q+N1                                                      
      K1=K1+2                                                           
      IF(K1.NE.0)GO TO 75                                               
      DT1A(KA1)=0.D0                                                    
      IF(I1X.NE.0)GO TO 85                                              
   75 DT1A(KA1)=DTU3R3(I1,J1,LAM1,MU1,IE1,J1T,M1T,K1,L1,M1)             
C     icount = icount + 1
CW    write(6,9000) ie1,j1t,m1t,k1,l1,m1,icount, 1, dt1a(ka1)
C9000 format(' Calling for DTU3R3 <',3i4,'|',3i4,'>', i10,'(',i1,')',
C    1   f12.5 )
      DS=0.D0                                                           
      DO 80 KAS=1,KA1                                                   
      KASKA1=KAS+KA1Q                                                   
   80 DS=DS+DON1(KASKA1)*DT1A(KAS)                                      
      DT1A(KA1)=DS                                                      
   85 CONTINUE                                                          
      K2=K2S                                                            
      KA2Q=-N2                                                          
      DO 100 KA2=1,KA2MAX                                               
      KA2Q=KA2Q+N2                                                      
      K2=K2+2                                                           
      IF(K2.NE.0)GO TO 90                                               
      DT2A(KA2)=0.D0                                                    
      IF(I2X.NE.0)GO TO 100                                             
   90 DT2A(KA2)=DTU3R3(I2,J2,LAM2,MU2,IE2,J2T,M2T,K2,L2,M2)             
C     icount = icount + 1
C     write(6,9000) ie2,j2t,m2t,k2,l2,m2,icount, 2, dt2a(ka2)
      DS=0.D0                                                           
      DO 95 KAS=1,KA2                                                   
      KASKA2=KAS+KA2Q                                                   
   95 DS=DS+DON2(KASKA2)*DT2A(KAS)                                      
      DT2A(KA2)=DS                                                      
  100 CONTINUE                                                          
      DC=DWR3(J1T,J2T,J3T,-M1T,-M2T,-M3T)                               
      KA2Q=-N1                                                          
      DO 105 KA2=1,KA2MAX                                               
      KA2Q=KA2Q+N1                                                      
      DO 105 KA1=1,KA1MAX                                               
      KA1KA2=KA1+KA2Q                                                   
  105 DS2A(KA1KA2)=DS2A(KA1KA2)+DC*DT1A(KA1)*DT2A(KA2)                  
  110 CONTINUE                                                          
      DC=DWR3(2*L1,2*L2,2*L3,2*M1,2*M2,2*K3)                            
      KA2Q=-N1                                                          
      DO 115 KA2=1,KA2MAX                                               
      KA2Q=KA2Q+N1                                                      
      KA23Q=KA2Q+KA3Q                                                   
      DO 115 KA1=1,KA1MAX                                               
      KA1KA2=KA1+KA2Q                                                   
      KA123=KA1+KA23Q                                                   
  115 DS1A(KA123)=DS1A(KA123)+DC*DS2A(KA1KA2)                           
  120 CONTINUE                                                          
      KA3P=-N012                                                        
      KA3Q=-N12                                                         
      DO 125 KA3=1,KA3MAX                                               
      KA3Q=KA3Q+N12                                                     
      KA3P=KA3P+N012                                                    
      KA2P=-N01                                                         
      KA2Q=-N1                                                          
      DO 125 KA2=1,KA2MAX                                               
      KA2Q=KA2Q+N1                                                      
      KA2P=KA2P+N01                                                     
      KA23P=KA2P+KA3P                                                   
      KA23Q=KA2Q+KA3Q                                                   
      KA1P=-N0                                                          
      DO 125 KA1=1,KA1MAX                                               
      KA1P=KA1P+N0                                                      
      KA123P=KA1P+KA23P                                                 
      KA123=KA1+KA23Q                                                   
      DO 125 KR0=1,KR0MAX                                               
      KR0123=KR0+KA123P                                                 
      KR0IND=KR0+INDP                                                   
  125 DWU3R3(KR0123)=DWU3R3(KR0123)+DEWU3(KR0IND)*DS1A(KA123)
  130 CONTINUE                                                          
      KA3P=-N012                                                        
      KA3Q=-N3                                                          
      DO 140 KA3=1,KA3MAX                                               
      KA3Q=KA3Q+N3                                                      
      KA3P=KA3P+N012                                                    
      KA2P=-N01                                                         
      DO 140 KA2=1,KA2MAX                                               
      KA2P=KA2P+N01                                                     
      KA23P=KA2P+KA3P                                                   
      KA1P=-N0                                                          
      DO 140 KA1=1,KA1MAX                                               
      KA1P=KA1P+N0                                                      
      KA123P=KA1P+KA23P                                                 
      KA12P=KA1P+KA2P                                                   
      DO 140 KR0=1,KR0MAX                                               
      KR0123=KR0+KA123P                                                 
      KR012=KR0+KA12P                                                   
      DS=0.D0                                                           
      KASP=-N012                                                        
      DO 135 KAS=1,KA3                                                  
      KASP=KASP+N012                                                    
      KASKA3=KAS+KA3Q                                                   
      KR012S=KR012+KASP                                                 
  135 DS=DS+DON3(KASKA3)*DWU3R3(KR012S)
  140 DWU3R3(KR0123)=DS                                                 
C     write(6,*) ' icount ', icount
      RETURN                                                            
C     DIMENSION CHECKS (LSU,6-81)-START                                 
  145 WRITE(6,150)N1                                                    
  150 FORMAT(37H *****CWU3R3 DIMENSION OVERFLOW:  N1=,I10)              
      GO TO 175                                                         
  155 WRITE(6,160)N2                                                    
  160 FORMAT(37H *****CWU3R3 DIMENSION OVERFLOW:  N2=,I10)              
      GO TO 175                                                         
  165 WRITE(6,170)N3                                                    
  170 FORMAT(37H *****CWU3R3 DIMENSION OVERFLOW:  N3=,I10)              
C     DIMENSION CHECKS (LSU,6-81)--STOP                                 
  175 STOP                                                              
      END                                                               
