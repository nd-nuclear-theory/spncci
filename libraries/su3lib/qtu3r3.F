      FUNCTION QTU3R3(I,J,LAMX,MUX,IEBX,JBT,MBTX,K,L,M)                 
C     ------------------------------------------------------------------
C     EXTREMAL U3R3 TRANSFORMATION--CHECKS OUT                          
C     ------------------------------------------------------------------
C     UPDATE/MOD: (LSU,11-78)  J.P.DRAAYER        EXPANDED RANGE        
C                 (LSU,05-80)  J.P.DRAAYER        INTEGER->REAL         
C                 (LSU,06-81)  J.P.DRAAYER        4->4.D0 IN DTS        
C                 (LSU,08-81)  J.P.DRAAYER        LOG BINOMIALS         
C                 (LSU,01-83)  J.P.DRAAYER        EXTENDED PRECISION    
C                 (LSU,03-83)  J.P.DRAAYER        SPACE SAVING MEASURE  
C                 (LSU,10-89)  J.P.DRAAYER        VS FORTRAN UPGRADE    
C                 (ND,11-11)   M.A.CAPRIO         DOUBLE/QUAD SWITCHABLE
C                                                                       
C     REFERENCES--J.P.DRAAYER AND Y.AKIYAMA, J.MATH.PHYS.14(1973)1904   
C                 J.P.DRAAYER, NUCL.PHYS.A129(1969)647-665              
C     PARAMETERS--(I,J) : (1,1)=GHW, (1,0)=GHW', (0,0)=GLW, (0,1)=GLW'  
C     ------------------------------------------------------------------
      IMPLICIT REAL*8(D)
#if defined(SU3DBL)
      IMPLICIT REAL*8(Q)
#elif defined(SU3QUAD)
      IMPLICIT REAL*16(Q)                                     
#endif
C                                                                       
      COMMON/BKQB/QBINO(18721),QBINV(18721),QTWOS(-192:192)             
C                                                                       
      ID(I)=I*(I+1)/2+1                                                 
      ND(I,J)=I*(I+1)/2+J+1                                             
C                                                                       
C     SET HIGHEST WEIGHT/LOWEST WEIGHT                                  
C                                                                       
      IF(I.EQ.0)THEN                                                    
         LAM=MUX                                                        
         MU=LAMX                                                        
         IEB=-IEBX                                                      
         MBT=-MBTX                                                      
      ELSE                                                              
         LAM=LAMX                                                       
         MU=MUX                                                         
         IEB=IEBX                                                       
         MBT=MBTX                                                       
      ENDIF                                                             
C                                                                       
C     OVERALL PHASE                                                     
C                                                                       
      IF(I.EQ.J)THEN                                                    
         IPH=0                                                          
      ELSE                                                              
         IPH=(LAM+K)/2                                                  
      ENDIF                                                             
C                                                                       
C     INITIALIZATIONS                                                   
C                                                                       
      QTU3R3=0.D0                                                       
      KBAR=IABS(K)                                                      
      MBAR=IABS(M)                                                      
      IAH=(LAM+KBAR)/2                                                  
      IBH=(JBT+MBAR)/2                                                  
      IPB=(2*(LAM-MU)+3*JBT-IEB)/3                                      
      QTS=QTWOS(-IPB)                                                   
      IPB=IPB/2                                                         
      IQB=MU+IPB-JBT                                                    
      IRB=(JBT+MBT)/2                                                   
      IF(K.LT.0)IPH=IPH+IPB+IQB+L+M                                     
      IF(M.LT.0)IPH=IPH+IQB+IRB+L+K                                     
      IW=L+KBAR                                                         
      IX=L-KBAR                                                         
      IY=L-MBAR                                                         
      IS=IW+IX                                                          
C                                                                       
C     OVERALL FACTOR                                                    
C                                                                       
      QTS=QTS*DFLOAT(IS+1)*                                             
#if defined(SU3DBL)
     1    DSQRT(QBINO(ND(LAM,IPB))*QBINO(ND(MU,IQB))*                   
#else defined(SU3QUAD)
     1    QSQRT(QBINO(ND(LAM,IPB))*QBINO(ND(MU,IQB))*                   
#endif
     2    QBINO(ND(LAM+MU+1,IQB))*QBINV(ND(IPB+MU+1,IQB))*              
     3    QBINV(ND(JBT,IRB))*QBINO(ND(IS,IX))*QBINV(ND(IS,IY)))         
      IF(BTEST(L-IPB+IPH,0))QTS=-QTS                                    
C                                                                       
C     SETUP CONSTANTS                                                   
C                                                                       
      NTB=ID(IPB)                                                       
      MA2=LAM                                                           
      MA1=0                                                             
      NA1=1                                                             
      NA2=ND(MA2,IAH)                                                   
      IB1=IRB-MBT                                                       
      MB2=IPB                                                           
      MB1=MU-IQB                                                        
      NB1=ID(MB1)                                                       
      NB2=NTB+IB1                                                       
      IB1=IB1+1                                                         
      IB2=IPB-IB1-1                                                     
      MC2=LAM+MU+L                                                      
      MC1=IPB+IQB                                                       
      NC1=ID(MC1)                                                       
      NC2=ND(MC2,(IX+IY-L+MC2-MC1)/2)                                   
      NC3=ID(IX)                                                        
      NC4=ND(IW,IY)                                                     
C                                                                       
C     OUTER LOOP                                                        
C                                                                       
      IZMAX=IY                                                          
      IZMIN=IZMAX-IW                                                    
      IF(IX.LT.IZMAX)IZMAX=IX                                           
      IF(IZMIN.LT.0)IZMIN=0                                             
      DO IW=0,IPB                                                       
         IF(K.EQ.0)THEN                                                 
            IF(BTEST(IW,0))GOTO 10                                      
         ENDIF                                                          
         IF(M.EQ.0)THEN                                                 
            IF(BTEST(IB1-MB2+1,0))GOTO 10                               
         ENDIF                                                          
C                                                                       
C        INNER LOOP "A"                                                 
C                                                                       
         QAS=0.D0                                                       
         IXMAX=IAH                                                      
         IXMIN=IXMAX-MA2                                                
         IF(MA1.LT.IXMAX)IXMAX=MA1                                      
         IF(IXMIN.LT.0)IXMIN=0                                          
         DO IX=IXMIN,IXMAX                                              
               QAX=QBINO(NA1+IX)*QBINO(NA2-IX)                          
            IF(BTEST(IX,0))THEN                                         
               QAS=QAS-QAX                                              
            ELSE                                                        
               QAS=QAS+QAX                                              
            ENDIF                                                       
         ENDDO                                                          
C                                                                       
C        INNER LOOP "B"                                                 
C                                                                       
         QBS=0.D0                                                       
         IXMAX=IB1-1                                                    
         IXMIN=IXMAX-MB2                                                
         IF(MB1.LT.IXMAX)IXMAX=MB1                                      
         IF(IXMIN.LT.0)IXMIN=0                                          
         MS2=2*(IXMIN+1)+IB2-IW                                         
         MS1=JBT-MS2                                                    
         NS1=ID(MS1)                                                    
         NS2=ND(MS2,IBH)                                                
         DO IX=IXMIN,IXMAX                                              
            IYMAX=IBH                                                   
            IYMIN=IYMAX-MS2                                             
            IF(MS1.LT.IYMAX)IYMAX=MS1                                   
            IF(IYMIN.LT.0)IYMIN=0                                       
            QS=0.D0                                                     
            DO IY=IYMIN,IYMAX                                           
                  QX=QBINO(NS1+IY)*QBINO(NS2-IY)                        
               IF(BTEST(IY,0))THEN                                      
                  QS=QS-QX                                              
               ELSE                                                     
                  QS=QS+QX                                              
               ENDIF                                                    
            ENDDO                                                       
            IF(IX.NE.IXMAX)THEN                                         
               NS1=NS1-2*MS1+1                                          
               MS1=MS1-2                                                
               MS2=MS2+2                                                
               NS2=NS2+2*MS2-1                                          
            ENDIF                                                       
            QBS=QBS+QBINO(NB1+IX)*QBINO(NB2-IX)*QS                      
         ENDDO                                                          
C                                                                       
C        "C" LOOP                                                       
C                                                                       
         QCS=0.D0                                                       
         IYMAX=MC1                                                      
         DO IX=IZMIN,IZMAX                                              
            QS=0.D0                                                     
            DO IY=0,IYMAX                                               
                  QX=QBINO(NC1+IY)*QBINV(NC2+IY-IX)                     
               IF(BTEST(IY,0))THEN                                      
                  QS=QS-QX                                              
               ELSE                                                     
                  QS=QS+QX                                              
               ENDIF                                                    
            ENDDO                                                       
               QCX=QBINO(NC3+IX)*QBINO(NC4-IX)*QS                       
            IF(BTEST(IX,0))THEN                                         
               QCS=QCS-QCX                                              
            ELSE                                                        
               QCS=QCS+QCX                                              
            ENDIF                                                       
         ENDDO                                                          
         QTU3R3=QTU3R3+QAS*QBS*QCS*QBINO(NTB+IW)/DFLOAT(MC2+1)          
   10 IF(IW.NE.IPB)THEN                                                 
         MA1=MA1+1                                                      
         NA1=NA1+MA1                                                    
         NA2=NA2-MA2                                                    
         MA2=MA2-1                                                      
         MB1=MB1+1                                                      
         NB1=NB1+MB1                                                    
         NB2=NB2-MB2                                                    
         MB2=MB2-1                                                      
         NC1=NC1-MC1                                                    
         MC1=MC1-1                                                      
         NC2=NC2-MC2                                                    
         MC2=MC2-1                                                      
      ENDIF                                                             
      ENDDO                                                             
      QTU3R3=QTS*QTU3R3                                                 
      RETURN                                                            
      END                                                               
