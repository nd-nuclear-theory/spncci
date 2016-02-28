      SUBROUTINE WU39LM(LA,MA,LB,MB,LC,MC,LD,MD,LE,ME,LF,               
     1                  MF,LG,MG,LH,MH,LI,MI,S,N9LMAX,*)                
C     ------------------------------------------------------------------
C     AUTHOR: MILLENER VIA RESKE VIA DRAAYER                            
C     UPDATE: 05/80  ORIGINAL CODE FOR THE MILLENER VIA RESKE 9LM LOGIC 
C             08/81  DRAAYER REWRITE OF THE FIRST VERSION, NAMED SU39LM 
C             03/83  GENERATED FULL DOCUMENTATION FOR EASE OF REFERENCE 
C             02/88  ADDED BINARY TREE TECHNOLOGY ... SU39LM --> YU39LM 
C             10/13  INTRODUCED ALTERNATE RETURN FOR UNALLOWED COUPLING 
C             11/89  UPDATED ... NEW BINARY TREE PACKAGE                
C             03/94  UPDATED ... WST (C. BAHRI) ...   YU39LM --> WU39LM 
C             07/95  COMMON BLOCK FOR WST ARRAYS (C. BAHRI)
C             03/11  MEMORY OUT-OF-BOUNDARY WRITE FIXED (T. DYTRYCH)
C     ------------------------------------------------------------------
C                                                                       
C                | (XA) (XB) (XC) R1 |   |  (X1)  (X2)   (X12)    R12  |
C     S(INDEX) = | (XD) (XE) (XF) R2 | = |  (X3)  (X4)   (X34)    R34  |
C                | (XG) (XH) (XI) R6 |   | (X13) (X24)   (LM0)  R13,24 |
C                |  R4   R5   R3     |   |  R13   R24   R12,34         |
C                                                                       
C                   INDEX = R1              INDEX = R12                 
C                   + NU*(R2-1)             + N1*(R34-1)                
C                   + NV*(R3-1)             + N1*N2*(R12,34-1)          
C                   + NW*(R4-1)             + N1*N2*N3*(R13-1)          
C                   + NX*(R5-1)             + N1*N2*N3*N4*(R34-1)       
C                   + NY*(R6-1)             + N1*N2*N3*N4*N5*(R13,24-1) 
C                                                                       
C                     WHERE --- N1=MULTU3(LA,MA,LB,MB,LC,MC)            
C                               N2=MULTU3(LD,MD,LE,ME,LF,MF)            
C                               N3=MULTU3(LC,MC,LF,MF,LI,MI)            
C                               N4=MULTU3(LA,MA,LD,MD,LG,MG)            
C                               N5=MULTU3(LB,MB,LE,ME,LH,MH)            
C                               N6=MULTU3(LG,MG,LH,MH,LI,MI)            
C                                                                       
C     ------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)                                          
      PARAMETER (NTAL=6561)
      DIMENSION S(1),URU3(NTAL),VRU3(NTAL),TRU3(NTAL)                   

      IMAX = N9LMAX;
C     INITIALIZE SUMMATION ARRAY                                        
      DO I=1,IMAX                                                       
         S(I)=0.D0                                                      
      ENDDO                                                             
CTD   INDEX MUST BE INITIALIZED. IF NOT OUT-OF-BOUND MEMORY WRITES OCCUR
      INDEX = IMAX
CTD   INDEX = # 9-(l m) COEFFICIENTS 

      LMB=LB+MB                                                         
      IXM=MIN0(MG,LB)                                                   
      IXMP1=IXM+1                                                       
      DO 25 IXP1=1,IXMP1                                                
         IYM=MIN0(LG,LMB-(IXP1-1))                                      
         IYMP1=IYM+1                                                    
      DO 25 IYP1=1,IYMP1                                                
         IZMI=MAX0(0,MB-MG+(IXP1-1)-(IYP1-1))                           
         IZMA=MIN0(MB,LMB-(IXP1-1)-(IYP1-1))                            
         IZMIP1=IZMI+1                                                  
         IZMAP1=IZMA+1                                                  
      DO 25 IZP1=IZMIP1,IZMAP1                                          
         LJ=LG+LMB-2*(IYP1-1)-(IXP1-1)-(IZP1-1)                         
         MJ=MG-MB-(IXP1-1)+(IYP1-1)+2*(IZP1-1)                          
         IP=(2*(MJ-MG)+LJ-LG-LB+MB)/3                                   
         IN=(LMB+LB-2*(LJ-LG)-MJ+MG)/3                                  
         MIM=MAX0(0,-IP,IN-LG)                                          
         IRHO=(IXP1-1)-MIM+1                                            
         IF(IRHO.NE.1)GOTO 25                                           
         KJE=MULTU3(LJ,MJ,LE,ME,LI,MI)                                  
         IF(KJE.EQ.0)GOTO 25                                            
         KCD=MULTU3(LC,MC,LD,MD,LJ,MJ)                                  
         IF(KCD.EQ.0)GOTO 25                                            
C                                                                       
         CALL WRU3OPTIMIZED(LG,MG,LB,MB,LI,MI,LE,ME,LJ,MJ,LH,MH,
     1             I1,I2,I3,I4,URU3,NTAL)
C                                                                       
         CALL WZU3OPTIMIZED(LB,MB,LA,MA,LJ,MJ,LD,MD,LC,MC,LG,MG,
     1             J1,J2,J3,J4,VRU3,NTAL)
C                                                                       
         CALL WRU3OPTIMIZED(LC,MC,LD,MD,LI,MI,LE,ME,LJ,MJ,LF,MF,
     1             K1,K2,K3,K4,TRU3,NTAL)
C                                                                       
         MM1=J1                                                         
         MM2=MM1*K3                                                     
         MM3=MM2*K4                                                     
         MM4=MM3*J3                                                     
         MM5=MM4*I3                                                     
         MM6=MM5*I4                                                     
C
         NWI=I1
         NXI=NWI*I2
         NYI=NXI*I3
         NWN=J1
         NXN=NWN*J2
         NYN=NXN*J3
         NWJ=K1
         NXJ=NWJ*K2
         NYJ=NXJ*K3
C
         DO 20 IR=1,J1                                                  
            JAQ=-NXJ                                                    
         DO 20 JR=1,K3                                                  
            JAQ=JAQ+NXJ                                                 
            JBQ=-NYJ                                                    
         DO 20 KR=1,K4                                                  
            JBQ=JBQ+NYJ                                                 
            JB=JBQ+JAQ                                                  
            NAQ=-NXN                                                    
         DO 20 LR=1,J3                                                  
            NAQ=NAQ+NXN                                                 
            NA=NAQ+IR                                                   
            IAQ=-NXI                                                    
         DO 20 MR=1,I3                                                  
            IAQ=IAQ+NXI                                                 
            IBQ=-NYI                                                    
         DO 20 NR=1,I4                                                  
            INDEX=IR+MM1*(JR-1)+MM2*(KR-1)+MM3*(LR-1)                   
     #                         +MM4*(MR-1)+MM5*(NR-1)                   
            IF(INDEX.GT.N9LMAX)GOTO 50                                  
            IBQ=IBQ+NYI                                                 
            IB=IBQ+IAQ                                                  
            TEMP=0.D0                                                   
            NBQ=-NYN                                                    
            DO 15 JS=1,I1                                               
               NBQ=NBQ+NYN                                              
               NB=NBQ+NA                                                
               IA=IB+JS                                                 
               JCQ=-NWJ                                                 
               ICQ=-NWI                                                 
            DO 15 KS=1,K2                                               
               JCQ=JCQ+NWJ                                              
               JC=JCQ+JB                                                
               ICQ=ICQ+NWI                                              
               IC=ICQ+IA                                                
               A=URU3(IC)                                               
               NCQ=-NWN                                                 
            DO 15 LS=1,K1                                               
               NCQ=NCQ+NWN                                              
               NC=NCQ+NB                                                
               JD=JC+LS                                                 
               TEMP=TEMP+A*VRU3(NC)*TRU3(JD)                            
   15       CONTINUE                                                    
            S(INDEX)=S(INDEX)+TEMP                                      
   20    CONTINUE                                                       
   25 CONTINUE                                                          
C
   40 RETURN                                                            
C                                                                       
C     OVERFLOW CHECK/STOP                                               
C                                                                       
   50 WRITE(6,55)INDEX,N9LMAX                                           
   55 FORMAT('0SU39LM DIMENSION OVERFLOW: INDEX=',I10,5X,'N9LMAX=',I10) 
      STOP                                                              
      END                                                               
