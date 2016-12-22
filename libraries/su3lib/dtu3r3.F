      FUNCTION DTU3R3(I,J,LAMX,MUX,IEBX,JBT,MBTX,K,L,M)                 
C     ------------------------------------------------------------------
C     EXTREMAL U3R3 TRANSFORMATION--CHECKS OUT                          
C     ------------------------------------------------------------------
C     UPDATE/MOD: (LSU,11-78)  J.P.DRAAYER        EXPANDED RANGE        
C                 (LSU,05-80)  J.P.DRAAYER        INTEGER->REAL         
C                 (LSU,06-81)  J.P.DRAAYER        4->4.D0 IN DTS        
C                 (LSU,08-81)  J.P.DRAAYER        LOG BINOMIALS         
C                 (LSU,03-83)  J.P.DRAAYER        SPACE AND PRECISION   
C                 (LSU,10-84)  J.P.DRAAYER        BINARY TREE TECHNOLOGY
C                 (LSU,01-88)  J.P.DRAAYER        SOON PARK BT CODE     
C                 (LSU,10-89)  J.P.DRAAYER        VS FORTRAN UPGRADE    
C                 (LSU,03-90)  J.P.DRAAYER        REVISED BT LABELING   
C                 (LSU,07-91)  J.P.DRAAYER        CONVERTED BT TO WST   
C                 (UoT,07-95)  C.BAHRI            COMMON BLOCK FOR WST
C                 (ND,11-11)   M.A.CAPRIO         DOUBLE/QUAD SWITCHABLE
C                 (ND,11-22)   A.E.MCCOY          ADDED FLAG SU3QUAD_GNU
C                                                                       
C     REFERENCES--J.P.DRAAYER AND Y.AKIYAMA, J.MATH.PHYS.14(1973)1904   
C                 J.P.DRAAYER, NUCL.PHYS.A129(1969)647-665              
C     PARAMETERS--(I,J) : (1,1)=GHW, (1,0)=GHW', (0,0)=GLW, (0,1)=GLW'  
C     ------------------------------------------------------------------
      IMPLICIT REAL*8(D)
#if defined(SU3DBL)
      IMPLICIT REAL*8(Q)
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
      IMPLICIT REAL*16(Q)
#endif
C                                                                       
      COMMON/BKDB/DBINO(8385),DBINV(8385),DTWOS(-128:128)               
C                                                                       
C     OTHER FUNCTIONS                                                   
C                                                                       
      ID(I)=I*(I+1)/2+1                                                 
      ND(I,J)=I*(I+1)/2+J+1                                             
C                                                                       
C     SET HIGHEST WEIGHT/LOWEST WEIGHT                                  
C                                                                       
      IF(I.NE.1)THEN                                                    
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
C     EXTENDED PRECISION OPTION INVOKED IF (LAM+MU+L) > 17              
C                                                                       
      IF(LAM+MU+L.GT.17)GOTO 40                                         
C                                                                       
C     INITIALIZATIONS                                                   
C                                                                       
      IPH=0                                                             
      DTU3R3=0.D0                                                       
      KBAR=IABS(K)                                                      
      MBAR=IABS(M)                                                      
      IAH=(LAM+KBAR)/2                                                  
      IBH=(JBT+MBAR)/2                                                  
      IPB=(2*(LAM-MU)+3*JBT-IEB)/3                                      
      DTS=DTWOS(-IPB)                                                   
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
      DTS=DTS*DFLOAT(IS+1)*                                             
     1    DSQRT(DBINO(ND(LAM,IPB))*DBINO(ND(MU,IQB))*                   
     2    DBINO(ND(LAM+MU+1,IQB))*DBINV(ND(IPB+MU+1,IQB))*              
     3    DBINV(ND(JBT,IRB))*DBINO(ND(IS,IX))*DBINV(ND(IS,IY)))         
      IF(BTEST(L-IPB+IPH,0))DTS=-DTS                                    
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
            IF(BTEST(IW,0))GOTO 30                                      
         ENDIF                                                          
         IF(M.EQ.0)THEN                                                 
            IF(BTEST(IB1-MB2+1,0))GOTO 30                               
         ENDIF                                                          
C                                                                       
C        INNER LOOP "A"                                                 
C                                                                       
         DAS=0.D0                                                       
         IXMAX=IAH                                                      
         IXMIN=IXMAX-MA2                                                
         IF(MA1.LT.IXMAX)IXMAX=MA1                                      
         IF(IXMIN.LT.0)IXMIN=0                                          
         DO IX=IXMIN,IXMAX                                              
               DAX=DBINO(NA1+IX)*DBINO(NA2-IX)                          
            IF(BTEST(IX,0))THEN                                         
               DAS=DAS-DAX                                              
            ELSE                                                        
               DAS=DAS+DAX                                              
            ENDIF                                                       
         ENDDO                                                          
C                                                                       
C        INNER LOOP "B"                                                 
C                                                                       
         DBS=0.D0                                                       
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
            DS=0.D0                                                     
            DO IY=IYMIN,IYMAX                                           
                  DX=DBINO(NS1+IY)*DBINO(NS2-IY)                        
               IF(BTEST(IY,0))THEN                                      
                  DS=DS-DX                                              
               ELSE                                                     
                  DS=DS+DX                                              
               ENDIF                                                    
            ENDDO                                                       
            IF(IX.NE.IXMAX)THEN                                         
               NS1=NS1-2*MS1+1                                          
               MS1=MS1-2                                                
               MS2=MS2+2                                                
               NS2=NS2+2*MS2-1                                          
            ENDIF                                                       
            DBS=DBS+DBINO(NB1+IX)*DBINO(NB2-IX)*DS                      
         ENDDO                                                          
C                                                                       
C        "C" LOOP                                                       
C                                                                       
         DCS=0.D0                                                       
         IYMAX=MC1                                                      
         DO IX=IZMIN,IZMAX                                              
            DS=0.D0                                                     
            DO IY=0,IYMAX                                               
                  DX=DBINO(NC1+IY)*DBINV(NC2+IY-IX)                     
               IF(BTEST(IY,0))THEN                                      
                  DS=DS-DX                                              
               ELSE                                                     
                  DS=DS+DX                                              
               ENDIF                                                    
            ENDDO                                                       
               DCX=DBINO(NC3+IX)*DBINO(NC4-IX)*DS                       
            IF(BTEST(IX,0))THEN                                         
               DCS=DCS-DCX                                              
            ELSE                                                        
               DCS=DCS+DCX                                              
            ENDIF                                                       
         ENDDO                                                          
         DTU3R3=DTU3R3+DAS*DBS*DCS*DBINO(NTB+IW)/DFLOAT(MC2+1)          
   30 IF(IW.NE.IPB)THEN                                                 
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
      DTU3R3=DTS*DTU3R3                                                 
      GOTO 70
   40 DTU3R3=QTU3R3(1,1,LAM,MU,IEB,JBT,MBT,K,L,M)                       
C                                                                       
C     ADD OVERALL PHASE FOR (I.NE.J) CASE                               
C                                                                       
   70 IF(I.NE.J.AND.BTEST((LAM+K)/2,0))THEN                             
         DTU3R3=-DTU3R3                                                 
      ENDIF                                                             
C                                                                       
C     EXIT CALCULATION                                                  
C                                                                       
      RETURN                                                            
      END                                                               
