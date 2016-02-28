      SUBROUTINE WZU3OPTIMIZED(LAM1,MU1,LAM2,MU2,LAM,MU,LAM3,MU3,LAM12,
     1MU12,LAM23,MU23,KR0A,KR0B,KR0C,KR0D,DZU3,NABCD)
C     ------------------------------------------------------------------
C     TRANSFORMATION FROM (AB)E,D;C  TO  (AD)F,B;C                      
C     ------------------------------------------------------------------
C     AUTHOR: MILLENER VIA RESKE VIA DRAAYER                            
C     UPDATE: 05/80  REWRITE (DRAAYER) OF THE MILLENER/RESKE CODE       
C             08/81  ADOPTED RUNNING INDEX FOR DEWU3 ... CZU3 --> XZU3  
C             02/88  INCLUDED BINARY TREE TECHNOLOGY ... XZU3 --> YZU3  
C             11/89  UPDATED ... BINARY TREE PACKAGE                    
C             03/94  UPDATED ... WST AS FOR PUBLICATION  YZU3 --> WZU3
C                                ELIMINATE EJR CHECK (UNNECESSARY)
C             07/95  COMMON BLOCK FOR WST ARRAYS (C. BAHRI)
C             07/11  WST & MULTU3 OVERHEAD REMOVED (T. DYTRYCH)
C     ------------------------------------------------------------------
C     PARAMERTERS--SEE ALSO XEWU3 AND XWU3                              
C       EXTERNAL--NABCD=MAX(KR0A*KR0B*KR0C*KR0D)
C       INTERNAL--X1=ABS(NA)                                            
C                 X2=ABS(NB)                                            
C                 X3=ABS(NC)                                            
C                 X4=ABS(ND)                                            
C                 X5=NX*(NX+1)*(NX+2)/6, NX=ABS(MAX(LAM1+MU1+1))        
C                 X6=NX*(NX+1)*(NX+2)/6, NX=ABS(MAX(LAM3+MU3+1))        
C                 X7=NX*(NX+1)*(NX+2)/6, NX=ABS(MAX(LAM3+MU3+1))        
C                 X8=NX*(NX+1)*(NX+2)/6, NX=ABS(MAX(LAM1+MU1+1))        
C*                X9=X1*X5   ->KIMAX1                                   
C*                X10=X2*X6  ->KIMAX1                                   
C*                X11=X3*X7  ->KIMAX1                                   
C*                X12=X4*X8  ->KIMAX1  (SEE XEWU3)                      
C                 X13=ABS(MAX(LAM1+MU1+1))                              
C                 X14=X13*X4                                            
C                 X15=ABS(MAX(DIM(LAM1,MU1)))                           
C*                X16=X1*X15 ->KIMAX2  (SEE XWU3)                       
C                 X17=ABS(MAX(IESMAX))                                  
C                     SAFE TO SET X17=ABS(MAX(LAM1+MU1+1))              
C                 X18=X17*ABS(MAX(J1SMAX))                              
C                     SAFE TO SET X18=X17*ABS(MAX(LAM1+MU1+1))          
C     EXTENSIONS--CHANGE EXTERNAL PARAMERTERS IN CALL                   
C                 ADJUST INTERNAL PARAMETERS BELOW                      
C*    DIMENSIONS--DZU3(NA*NB*NC*ND),                                    
C                 JXTA(X5),JYTA(X5),IEA(X5),JXTB(X6),JYTB(X6),IEB(X6),  
C                 JXTC(X7),JYTC(X7),IEC(X7),JXTD(X8),JYTD(X8),IED(X8),  
C*                DEWU3A(X9),DEWU3B(X10),DEWU3C(X11),DEWU3D(X12),       
C*                IA(X13),DA(X14),DB(X4),DT(X4),DWU3(X16),J2SMAX(X18),  
C                 J2TMAX(X18),J1SMAX(X17),J1TMAX(X17),INDMAT(X18)       
C       COMMENTS--ASSUME MAX NA=9,NB=9,NC=9,ND=9                        
C                        SET X'S AS IN DATA STATEMENT                   
C     ------------------------------------------------------------------
      IMPLICIT REAL*8(D),INTEGER(X)                                     
      EQUIVALENCE(JXTA(1),JXTC(1),JXTD(1))                              
      EQUIVALENCE(JYTA(1),JYTC(1),JYTD(1))                              
      EQUIVALENCE(IEA(1),IEC(1),IED(1))                                 
      DIMENSION DZU3(1),                                                
     1JXTA(13244),JYTA(13244),IEA(13244),JXTB(13244),JYTB(13244),       
     2IEB(13244),                                                       
     3JXTC(13244),JYTC(13244),IEC(13244),JXTD(13244),JYTD(13244),       
     4IED(13244)                                                        
      COMMON/BKSAVE/                                                    
     1DEWU3A(39732),DEWU3B(39732),DEWU3C(39732),DEWU3D(39732)           
      DIMENSION                                                         
     1IA(42),DA(378),DB(9),DT(9),DWU3(27090),J2SMAX(1764),              
     2J2TMAX(1764),J1SMAX(42),J1TMAX(42),INDMAT(1764)                   
      DATA X1,X2,X3,X4,X5,X6,X7,X8,X9,X10                               
     1/9,9,9,9,13244,13244,13244,13244,39732,39732/                     
      DATA X11,X12,X13,X14,X15,X16,X17,X18                              
     1/39732,39732,42,378,9030,27090,42,1764/                           
C     OTHER FUNCTIONS                                                   
      INDEX(J1TD,LAM1,J1T,J2TD,LAM2,J2T)=1+J2TD*(J2TD+1)*(3*J1TD+J2TD+5)
     1/6+(J1TD+1)*(LAM2+J2TD-J2T)/2+(LAM1+J1TD-J1T)/2                   

      NA=KR0A                                                           
      NAB=NA*KR0B                                                       
      NABC=NAB*KR0C                                                     
      KIMAX1=X9                                                         
      KIMAX2=X16                                                        
      CALL XEWU3(LAM23,MU23,LAM1,MU1,LAM,MU,1,NECD,KR0D,INDMAX,         
     1DEWU3D,JXTD,JYTD,IED,X4,X8,KIMAX1)                                
      NNCD=NECD+1                                                       
      IDQ=(INDMAX-NNCD-1)*KR0D                                          
      DO 15 JD=1,NNCD                                                   
      IDQ=IDQ+KR0D                                                      
      IA(JD)=JD                                                         
      KDQ=-X13                                                          
      DO 15 KD=1,KR0D                                                   
      KDQ=KDQ+X13                                                       
      JDKD=JD+KDQ                                                       
      KDID=KD+IDQ                                                       
   15 DA(JDKD)=DEWU3D(KDID)                                             
      CALL DLUT(NNCD,KR0D,IA,DA,X13)                                    
      CALL XEWU3(LAM2,MU2,LAM1,MU1,LAM12,MU12,1,NECA,KR0A,INDMAX,       
     1DEWU3A,JXTA,JYTA,IEA,X1,X5,KIMAX1)                                
      CALL XEWU3(LAM2,MU2,LAM3,MU3,LAM23,MU23,1,NECC,KR0C,INDMAX,       
     1DEWU3C,JXTC,JYTC,IEC,X3,X7,KIMAX1)                                
      CALL XEWU3(LAM12,MU12,LAM3,MU3,LAM,MU,1,NECB,KR0B,INDMAX,         
     1DEWU3B,JXTB,JYTB,IEB,X2,X6,KIMAX1)                                
      I1=LAM+2*MU                                                       
      I2=LAM12+2*MU12                                                   
      I3=4*LAM12+2*MU12                                                 
      I4=2*I2                                                           
      I5=2*(LAM12-MU12)                                                 
      I6=LAM3+2*MU3                                                     
      I7=2*(LAM3-MU3)                                                   
      IE1=LAM23+2*MU23-I1                                               
      J1TS=LAM1+NECD+2                                                  
      J3TDMA=MIN0(NECB,NECC)+1                                          
      DO 35 J3S=1,J3TDMA                                                
      J3TD=J3S-1                                                        
      IE3=-I6+3*J3TD                                                    
      JJ3TA=LAM3+J3TD                                                   
      IS=I6-J3TD                                                        
      IF(IS.LT.JJ3TA)JJ3TA=IS                                           
      JJ3TA=JJ3TA+1                                                     
      IS=(I7-IE3)/3                                                     
      JJ3TB=JJ3TA-IABS(IS)                                              
      J12TD=NECB-J3TD                                                   
      J2TD=NECC-J3TD                                                    
      IE12=-I2+3*J12TD                                                  
      JJ12TA=I3+IE12                                                    
      IS=I4-IE12                                                        
      IF(IS.LT.JJ12TA)JJ12TA=IS                                         
      JJ12TA=JJ12TA/3+1                                                 
      IS=(I5-IE12)/3                                                    
      JJ12TB=JJ12TA-IABS(IS)                                            
      DO 35 JJ12T=1,JJ12TB,2                                            
      J12T=JJ12TA-JJ12T                                                 
      CALL XWU3(LAM2,MU2,LAM1,MU1,LAM12,MU12,IE12,J12T,NECA,DEWU3A,     
     1KR0A,INDMAX,DWU3,J2SMAX,J2TMAX,J1SMAX,J1TMAX,IESMAX,IE1MAX,INDMAT,
     2X1,X17,X15,KIMAX2)                                                
      IE1TES=IE1MAX-IE1                                                 
      IF(IE1TES.LT.0)GOTO 35                                            
      IF(IE1TES.GT.3*(IESMAX-1))GOTO 35                                 
      IES=(IE1-IE1MAX+3*IESMAX)/3                                       
      KDQ=-NABC                                                         
      DO 30 KD=1,KR0D                                                   
      KDQ=KDQ+NABC                                                      
      J1T=J1TS-2*IA(KD)                                                 
      J1S=(J1TMAX(IES)+2-J1T)/2                                         
      IF(J1S.LT.1)GOTO 30                                               
      IF(J1S.GT.J1SMAX(IES))GOTO 30                                     
      J1SQ=X17*(J1S-1)                                                  
      IESJ1S=IES+J1SQ                                                   
      J2T=J2TMAX(IESJ1S)+2                                              
      INDA=(INDMAT(IESJ1S)-J2T)/2                                       
      J2SB=J2SMAX(IESJ1S)                                               
      DO 25 J2S=1,J2SB                                                  
      J2T=J2T-2                                                         
      IAQ=INDA*KR0A                                                     
      INDA=INDA+1                                                       
      DO 25 JJ3T=1,JJ3TB,2                                              
      J3T=JJ3TA-JJ3T                                                    
      INDC=INDEX(J2TD,LAM2,J2T,J3TD,LAM3,J3T)                           
      IF(JXTC(INDC).LT.0)GOTO 25                                        
      INDB=INDEX(J12TD,LAM12,J12T,J3TD,LAM3,J3T)                        
      IF(JXTB(INDB).LT.0)GOTO 25                                        
      IPH=J2T+LAM+3*(J12T+LAM23)                                        
      D2=DFLOAT(1-MOD(IABS(IPH),4))                                     
      DC=D2*DRR3(J1T,J2T,LAM,J3T,J12T,LAM23)                            
      DC=DC*DSQRT(DFLOAT((J12T+1)*(LAM23+1)))                           
      ICQ=(INDC-1)*KR0C                                                 
      IBQ=(INDB-1)*KR0B                                                 
      KCQ=-NAB                                                          
      DO 20 KC=1,KR0C                                                   
      KCQ=KCQ+NAB                                                       
      KCDQ=KCQ+KDQ                                                      
      KCIC=KC+ICQ                                                       
      KBQ=-NA                                                           
      DO 20 KB=1,KR0B                                                   
      KBQ=KBQ+NA                                                        
      KBCDQ=KBQ+KCDQ                                                    
      KBIB=KB+IBQ                                                       
      DO 20 KA=1,KR0A                                                   
      KABCD=KA+KBCDQ                                                    
      KAIA=KA+IAQ                                                       
   20 DZU3(KABCD)=DZU3(KABCD)+DC*DEWU3C(KCIC)*DEWU3B(KBIB)*DWU3(KAIA)   
   25 CONTINUE                                                          
   30 CONTINUE                                                          
   35 CONTINUE                                                          
      KCQ=-NAB                                                          
      DO 60 KC=1,KR0C                                                   
      KCQ=KCQ+NAB                                                       
      KBQ=-NA                                                           
      DO 60 KB=1,KR0B                                                   
      KBQ=KBQ+NA                                                        
      KBCQ=KBQ+KCQ                                                      
      DO 60 KA=1,KR0A                                                   
      KABC=KA+KBCQ                                                      
      KDQ=-NABC                                                         
      DO 40 KD=1,KR0D                                                   
      KDQ=KDQ+NABC                                                      
      KABCD=KABC+KDQ                                                    
   40 DB(KD)=DZU3(KABCD)                                                
      IF(KR0D.GT.1)GOTO 45                                              
      DB(1)=DB(1)/DA(1)                                                 
      GOTO 50                                                           
   45 CALL DBSR(KR0D,DA,DB,DT,X13)                                      
   50 KDQ=-NABC                                                         
      DO 55 KD=1,KR0D                                                   
      KDQ=KDQ+NABC                                                      
      KABCD=KABC+KDQ                                                    
   55 DZU3(KABCD)=DB(KD)                                                
   60 CONTINUE                                                          
   75 RETURN                                                            
      END                                                               
