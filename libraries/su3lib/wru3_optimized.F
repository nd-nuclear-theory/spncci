      SUBROUTINE WRU3OPTIMIZED(LAM1,MU1,LAM2,MU2,LAM,MU,LAM3,MU3,LAM12,
     1MU12,LAM23,MU23,KR0A,KR0B,KR0C,KR0D,DRU3,NABCD)                  
C     ------------------------------------------------------------------
C     RACAH COEFFICIENTS FOR U3 (X PREFIX FOR 6-81 VERSION)             
C     ------------------------------------------------------------------
C     AUTHOR: J.P.DRAAYER                                               
C     UPDATE: 06/81  ADOPTED RUNNING INDEX FOR DEWU3 ... CRU3 --> XRU3  
C             02/88  INCLUDED BINARY TREE TECHNOLOGY ... XRU3 --> YRU3  
C             11/89  UPDATED ... BINARY TREE PACKAGE                    
C             07/91  UPDATED ... WST & DRU3 INDEXING ... YRU3 --> WRU3  
C             03/94  UPDATED ... WST AS FOR PUBLICATION (C. BAHRI)
C             07/95  COMMON BLOCK FOR WST ARRAYS (C. BAHRI)
C             07/11  WST & MULTU3 OVERHEAD REMOVED (T. DYTRYCH)
C             04/16  MADE COMMON BLOCKS THREAD PRIVATE (A. MCCOY)
C                                                                       
C     REFERENCES--J.P.DRAAYER AND Y.AKIYAMA, J.MATH.PHYS.14(1973)1904   
C                 K.T.HECHT, NUCL.PHYS.62(1965)1                        
C     PARAMETERS--SEE ALSO XEWU3 AND XWU3                               
C       EXTERNAL--NABCD=MAX(KR0A*KR0B*KR0C*KR0D)                        
C       INTERNAL--X1=ABS(NA)                                            
C                 X2=ABS(NB)                                            
C                 X3=ABS(NC)                                            
C                 X4=ABS(ND)                                            
C                 X5=NX*(NX+1)*(NX+2)/6, NX=ABS(MAX(LAM2+MU2+1))        
C                 X6=NX*(NX+1)*(NX+2)/6, NX=ABS(MAX(LAM3+MU3+1))        
C                 X7=NX*(NX+1)*(NX+2)/6, NX=ABS(MAX(LAM3+MU3+1))        
C                 X8=NX*(NX+1)*(NX+2)/6, NX=ABS(MAX(LAM23+MU23+1))      
C*                X9=X1*X5    ->KIMAX1                                  
C*                X10=X2*X6   ->KIMAX1                                  
C*                X11=X3*X7   ->KIMAX1                                  
C*                X12=X4*X8   ->KIMAX1 (SEE XEWU3)                      
C                 X13=ABS(MAX(LAM23+MU23+1))                            
C                 X14=X13*X4                                            
C                 X15=ABS(MAX(DIM(LAM3,MU3)))                           
C*                X16=X3*X15  ->KIMAX2 (SEE XWU3)                       
C                 X17=ABS(MAX(IESMAX))                                  
C                     SAFE TO SET X17=ABS(MAX(LAM3+MU3+1))              
C                 X18=X17*ABS(MAX(J3SMAX))                              
C                     SAFE TO SET X18=X17*ABS(MAX(LAM3+MU3+1))          
C     EXTENSIONS--CHANGE EXTERNAL PARAMETERS IN CALL                    
C                 ADJUST INTERNAL PARAMETERS BELOW                      
C     DIMENSIONS--DRU3(NABCD),                                          
C                 JXTA(X5),JYTA(X5),IEA(X5),JXTB(X6),JYTB(X6),IEB(X6),  
C                 JXTC(X7),JYTC(X7),IEC(X7),JXTD(X8),JYTD(X8),IED(X8),  
C*                DEWU3A(X9),DEWU3B(X10),DEWU3C(X11),DEWU3D(X12),       
C*                IA(X13),DA(X14),DB(X4),DT(X4),DWU3(X16),J2SMAX(X18),  
C                 J2TMAX(X18),J3SMAX(X17),J3TMAX(X17),INDMAT(X18)       
C       COMMENTS--SET X'S AS IN DATA STATEMENT                          
C     ------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER LAM1,MU1,LAM2,MU2,LAM,MU,LAM3,MU3,LAM12,MU12,LAM23,MU23
      INTEGER KR0A,KR0B,KR0C,KR0D
      INTEGER J1TD, J1T, J2TD, J2T
      INTEGER I1, I2, I3, I4, I5, I6, IA, IAQ, IBQ, ICQ, IE12, IE23
      INTEGER INDA, INDB, INDC, INDMAT, INDMAX
      INTEGER IS, J12T, J12TD, J23T, J23TS, J2S, J2SB, J2SMAX, J2TMAX
      INTEGER J3S, J3SB, J3SMAX, J3SQ, J3T, J3TD, J3TMAX, JJ12T, JJ12TA
      INTEGER JJ12TB, JTDMAX, JTDMIN, JXTA, JXTB, JXTC, JXTD
      INTEGER IE3MAX, IEA, IEB, IEC, IED, IES, IESJ3S, IESMAX
      INTEGER NABCD
      INTEGER NA, NAB, NABC, NNCD, IDQ, JD, KD, JDKD, KDID
      INTEGER X1,X2,X3,X4,X5,X6,X7,X8,X9,X10
      INTEGER X11,X12,X13,X14,X15,X16,X17,X18
      INTEGER JYTA, JYTB, JYTC, JYTD, KA, KABC, KABCD, KAIA
      INTEGER KB, KBCDQ, KBCQ, KBIB, KBQ, KC, KCDQ, KCIC
      INTEGER KCQ, KDQ, KIMAX1, KIMAX2, NECA, NECB, NECC, NECD
      INTEGER IDM, INDEX
      REAL*8 DRU3, D1, D2, DA, DB, DC, DT
      REAL*8 DWU3, DEWU3A, DEWU3B, DEWU3C, DEWU3D
      REAL*8 DRR3

      EQUIVALENCE(JXTB(1),JXTC(1),JXTD(1))                              
      EQUIVALENCE(JYTB(1),JYTC(1),JYTD(1))                              
      EQUIVALENCE(IEB(1),IEC(1),IED(1))                                 
      DIMENSION DRU3(1),                                                
     1JXTA(13244),JYTA(13244),IEA(13244),JXTB(13244),JYTB(13244),       
     2IEB(13244),                                                       
     3JXTC(13244),JYTC(13244),IEC(13244),JXTD(13244),JYTD(13244),       
     4IED(13244)                                                        
      COMMON/BKSAVE/                                                    
     1DEWU3A(39732),DEWU3B(39732),DEWU3C(39732),DEWU3D(39732)
!$omp threadprivate(/BKSAVE/)           
      DIMENSION                                                         
     1IA(42),DA(378),DB(9),DT(9),DWU3(27090),J2SMAX(1764),              
     2J2TMAX(1764),J3SMAX(42),J3TMAX(42),INDMAT(1764)                   
      DATA X1,X2,X3,X4,X5,X6,X7,X8,X9,X10                               
     1/9,9,9,9,13244,13244,13244,13244,39732,39732/                     
      DATA X11,X12,X13,X14,X15,X16,X17,X18                              
     1/39732,39732,42,378,9030,27090,42,1764/                           
C     OTHER FUNCTIONS                                                   
      IDM(LAM,MU)=(LAM+1)*(MU+1)*(LAM+MU+2)/2                           
      INDEX(J1TD,LAM1,J1T,J2TD,LAM2,J2T)=1+J2TD*(J2TD+1)*(3*J1TD+J2TD+5)
     1/6+(J1TD+1)*(LAM2+J2TD-J2T)/2+(LAM1+J1TD-J1T)/2                   

C For debugging purposes ...
C      write(*,*) LAM1,MU1,LAM2,MU2,LAM,MU,LAM3,MU3,LAM12,MU12,LAM23,MU23
C      write(*,*) KR0A,KR0B,KR0C,KR0D
C      write(*,*) J1TD, J1T, J2TD, J2T
C      write(*,*) I1, I2, I3, I4, I5, I6, IA, IAQ, IBQ, ICQ, IE12, IE23
C      write(*,*) INDA, INDB, INDC, INDMAT, INDMAX
C      write(*,*) IS, J12T, J12TD, J23T, J23TS, J2S, J2SB, J2SMAX, J2TMAX
C      write(*,*) J3S, J3SB, J3SMAX, J3SQ, J3T, J3TD, J3TMAX, JJ12T
C      write(*,*) JJ12TB, JTDMAX, JTDMIN, JXTA, JXTB, JXTC, JXTD
C      write(*,*) IE3MAX, IEA, IEB, IEC, IED, IES, IESJ3S, IESMAX
C      write(*,*) NABCD
C      write(*,*) NA, NAB, NABC, NNCD, IDQ, JD, KD, JDKD, KDID
C      write(*,*) X1,X2,X3,X4,X5,X6,X7,X8,X9,X10
C      write(*,*) X11,X12,X13,X14,X15,X16,X17,X18
c      write(*,*) JYTA, JYTB, JYTC, JYTD, KA, KABC, KABCD, KAIA
C      write(*,*) KB, KBCDQ, KBCQ, KBIB, KBQ, KC, KCDQ, KCIC
C      write(*,*) KCQ, KDQ, KIMAX1, KIMAX2, NECA, NECB, NECC, NECD
C      write(*,*) DRU3, D1, D2, DA, DB, DC, DT
C      write(*,*) DWU3, DEWU3A, DEWU3B, DEWU3C, DEWU3D

      NA=KR0A                                                           
      NAB=NA*KR0B                                                       
      NABC=NAB*KR0C                                                     
      KIMAX1=X9                                                         
      KIMAX2=X16                                                        
      CALL XEWU3(LAM1,MU1,LAM23,MU23,LAM,MU,1,NECD,KR0D,INDMAX,DEWU3D,  
     1JXTD,JYTD,IED,X4,X8,KIMAX1)                                       
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
      CALL XEWU3(LAM2,MU2,LAM3,MU3,LAM23,MU23,1,NECC,KR0C,INDMAX,       
     1DEWU3C,JXTC,JYTC,IEC,X3,X7,KIMAX1)                                
      CALL XEWU3(LAM12,MU12,LAM3,MU3,LAM,MU,1,NECB,KR0B,INDMAX,         
     1DEWU3B,JXTB,JYTB,IEB,X2,X6,KIMAX1)                                
      CALL XEWU3(LAM12,MU12,MU2,LAM2,LAM1,MU1,1,NECA,KR0A,INDMAX,       
     1DEWU3A,JXTA,JYTA,IEA,X1,X5,KIMAX1)                                
      I1=LAM+2*MU                                                       
      I2=LAM12+2*MU12                                                   
      I3=4*LAM12+2*MU12                                                 
      I4=2*I2                                                           
      I5=2*(LAM12-MU12)                                                 
      JTDMIN=MAX0(0,NECA-LAM2-MU2,NECB-LAM3-MU3)                        
      JTDMAX=MIN0(NECA,NECB,LAM12+MU12)                                 
      D1=DFLOAT((LAM1+1)*IDM(LAM12,MU12))/DFLOAT(IDM(LAM1,MU1))         
      IE23=LAM1+2*MU1-I1                                                
      J23TS=LAM23+NECD+2                                                
      KDQ=-NABC                                                         
      DO 30 KD=1,KR0D                                                   
      KDQ=KDQ+NABC                                                      
      J23T=J23TS-2*IA(KD)                                               
      D2=DSQRT(D1*DFLOAT(J23T+1))                                       
      CALL XWU3(LAM2,MU2,LAM3,MU3,LAM23,MU23,IE23,J23T,NECC,DEWU3C,     
     1KR0C,INDMAX,DWU3,J2SMAX,J2TMAX,J3SMAX,J3TMAX,IESMAX,IE3MAX,INDMAT,
     2X3,X17,X15,KIMAX2)                                                
      IE12=3*IESMAX-IE3MAX-I1                                           
      J12TD=(IE12+I2)/3                                                 
      DO 30 IES=1,IESMAX                                                
      IE12=IE12-3                                                       
      J12TD=J12TD-1                                                     
      IF(J12TD.LT.JTDMIN)GOTO 30                                        
      IF(J12TD.GT.JTDMAX)GOTO 30                                        
      JJ12TA=I3+IE12                                                    
      IS=I4-IE12                                                        
      IF(IS.LT.JJ12TA)JJ12TA=IS                                         
      JJ12TA=JJ12TA/3+1                                                 
      IS=(I5-IE12)/3                                                    
      JJ12TB=JJ12TA-IABS(IS)                                            
      I6=2*LAM1+IS                                                      
      J2TD=NECA-J12TD                                                   
      J3TD=NECB-J12TD                                                   
      J3T=J3TMAX(IES)+2                                                 
      J3SB=J3SMAX(IES)                                                  
      J3SQ=-X17                                                         
      DO 25 J3S=1,J3SB                                                  
      J3SQ=J3SQ+X17                                                     
      IESJ3S=IES+J3SQ                                                   
      J3T=J3T-2                                                         
      J2T=J2TMAX(IESJ3S)+2                                              
      INDC=(INDMAT(IESJ3S)-J2T)/2                                       
      J2SB=J2SMAX(IESJ3S)                                               
      DO 25 J2S=1,J2SB                                                  
      J2T=J2T-2                                                         
      ICQ=INDC*KR0C                                                     
      INDC=INDC+1                                                       
      DO 25 JJ12T=1,JJ12TB,2                                            
      J12T=JJ12TA-JJ12T                                                 
      INDA=INDEX(J12TD,LAM12,J12T,J2TD,MU2,J2T)                         
      IF(JXTA(INDA).LT.0)GOTO 25                                        
      INDB=INDEX(J12TD,LAM12,J12T,J3TD,LAM3,J3T)                        
      IF(JXTB(INDB).LT.0)GOTO 25                                        
      DC=D2*DRR3(LAM1,J2T,LAM,J3T,J12T,J23T)                            
      IS=J12T+I6                                                        
      IF(4*(IS/4).NE.IS)DC=-DC                                          
      IAQ=(INDA-1)*KR0A                                                 
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
   20 DRU3(KABCD)=DRU3(KABCD)+DC*DEWU3A(KAIA)*DEWU3B(KBIB)*DWU3(KCIC)   
   25 CONTINUE                                                          
   30 CONTINUE                                                          
      KCQ=-NAB                                                          
      DO 55 KC=1,KR0C                                                   
      KCQ=KCQ+NAB                                                       
      KBQ=-NA                                                           
      DO 55 KB=1,KR0B                                                   
      KBQ=KBQ+NA                                                        
      KBCQ=KBQ+KCQ                                                      
      DO 55 KA=1,KR0A                                                   
      KABC=KA+KBCQ                                                      
      KDQ=-NABC                                                         
      DO 35 KD=1,KR0D                                                   
      KDQ=KDQ+NABC                                                      
      KABCD=KABC+KDQ                                                    
   35 DB(KD)=DRU3(KABCD)                                                
      IF(KR0D.GT.1)GOTO 40                                              
      DB(1)=DB(1)/DA(1)                                                 
      GOTO 45                                                           
   40 CALL DBSR(KR0D,DA,DB,DT,X13)                                      
   45 KDQ=-NABC                                                         
      DO 50 KD=1,KR0D                                                   
      KDQ=KDQ+NABC                                                      
      KABCD=KABC+KDQ                                                    
   50 DRU3(KABCD)=DB(KD)                                                
   55 CONTINUE                                                          
   70 RETURN                                                            
      END                                                               
