      SUBROUTINE WU3R3W(LM1,MU1,LM2,MU2,LM3,MU3,L1,L2,L3,               
     X                  K0MAX,K1MAX,K2MAX,K3MAX,DVALU,*)                
C     ------------------------------------------------------------------
C     AUTHOR: J.P.DRAAYER                                               
C     UPDATE: 05/83  ORIGINAL VERSION OF THE CODE CALLED XU3R3X         
C             02/88  ADDED BINARY TREE TECHNOLOGY AND RENAMED YU3R3W    
C             10/89  INTRODUCED ALTERNATE RETURN FOR UNALLOWED COUPLING 
C             11/89  UPDATED ... NEW BINARY TREE PACKAGE                
C             03/90  FUNCTION MULTU3 REPLACED BY SUBROUTINE U3MULT      
C             12/91  UPDATED ... WST (C.BAHRI) ... YU3R3W --> WU3R3W    
C             07/95  COMMON BLOCK FOR WST ARRAYS (C. BAHRI)
C             08/11  REMOVE BINARY TREE (T. DYTRYCH)
C             03/16  ADDED FLAG SU3LM82 FOR CHANGING LIMITS ON LAMBDA+MU<42 
C                    TO <82 (A.E.MCCOY)
C     ------------------------------------------------------------------
C     PROGRAM TO GENERATE U3R3 WIGNER COEFFICIENTS                      
C          NOTE:  IF THE ROUTINES IN THE SU3 PACKAGE ARE MODIFIED,      
C                 NCE1,NCE2,NCEX & NCW0,NCW1,NCW2,NCW3 AS WELL AS THE   
C                 DIMENSIONS USED IN THIS SAMPLE CONTROL ROUTINE MUST   
C                 ALSO BE ADJUSTED.  THE APPROPRIATE CHANGES FOLLOW     
C                 FROM COMMENT STATEMENTS DEFINING PARAMETERS IN THE    
C                 CALLING SEQUENCES OF XEWU3 (NCE'S) AND XWU3R3 (NCW'S).
C          NOTE:  SAME AS U3R3WG BUT WITH I,J SELECTION AUTOMATICALLY   
C                 MADE USING THE IJFACT ROUTINE--SEE U3R3WIG SUBROUTINE.
C          NOTE:  THIS VERSION USES THE "X" VERSIONS OF THE SU3 PACKAGE.
C     ------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(D)                                      
      DIMENSION DON1(81),DON2(81),DON3(81),DEWU3(39732),                
     1DWU3R3(6561),J1TA(13244),J2TA(13244),IEA(13244),DVALU(9,9,9,9)    
      DATA KIMAX1/39732/                                                
C SUPPORT FOR LARGE LAMBDA+MU REQUIRES C. BAHRI MOD IN XEWU3.F ON 4/97
#if defined(SU3LM82)
      DATA NCE1,NCE2,NCEX,NCW0,NCW1,NCW2,NCW3/9,95284,82,9,9,9,9/
#else
      DATA NCE1,NCE2,NCEX,NCW0,NCW1,NCW2,NCW3/9,13244,42,9,9,9,9/       
#endif

C     COMBINATIONS OF THE ABOVE: NCW01=NCW0*NCW1,NCW012=NCW01*NCW02,... 
      DATA NCW01,NCW012,NXALL/81,729,6561/                              
C     MULTIPLICITY FUNCTION FOR L-VALUES IN (LM,MU) SU3 IRREP           
      MULT(LM,MU,L)=MAX0(0,(LM+MU+2-L)/2)-MAX0(0,(LM+1-L)/2)-           
     1 MAX0(0,(MU+1-L)/2)                                               
C     SET EVERYTHING TO ZERO INITIALLY --> SAFETY MEASURE               
      ITYPE=0                                                           
C     ZERO OUT FEATURE CANCELLED (10/89)                                
C     K0MAX=0                                                           
C     K1MAX=0                                                           
C     K2MAX=0                                                           
C     K3MAX=0                                                           
C     DO 10 K3=1,9                                                      
C     DO 10 K2=1,9                                                      
C     DO 10 K1=1,9                                                      
C     DO 10 K0=1,9                                                      
C  10 DVALU(K0,K1,K2,K3)=0.E0                                           
C     CHECK IF THIS IS AN ALLOWED CASE                                  
C     K0MAX=MULTU3(LM1,MU1,LM2,MU2,LM3,MU3)                             
C     IF(K0MAX.EQ.0)GOTO 60                                             
      CALL U3MULT(LM1,MU1,LM2,MU2,LM3,MU3,K0MAX,*60)
C     CHECK ANGULAR MOMENTUM CONTENT                                    
      K1MAX=MULT(LM1,MU1,L1)                                            
      IF(K1MAX.EQ.0)GOTO 60
      K2MAX=MULT(LM2,MU2,L2)                                            
      IF(K2MAX.EQ.0)GOTO 60
      K3MAX=MULT(LM3,MU3,L3)                                            
      IF(K3MAX.EQ.0)GOTO 60
C     CHECK ANGULAR MOMENTUM COUPLING                                   
      IF(L1+L2.LT.L3)GOTO 60
      IF(IABS(L1-L2).GT.L3)GOTO 60
      KXALL=K0MAX*K1MAX*K2MAX*K3MAX                                     
      IF(K0MAX.GT.NCE1)GOTO 30                                          
      ITYPE=1                                                           
      IF(KXALL.GT.NXALL)GOTO 30                                         
      ITYPE=4                                                           
      IF(LM2+MU2+1.GT.NCEX)GOTO 30                                      
      ITYPE=7                                                           
      CALL IJFACT(LM1,MU1,I1,J1,0,1)                                    
      CALL IJFACT(LM2,MU2,I2,J2,0,1)                                    
      CALL IJFACT(LM3,MU3,I3,J3,0,1)                                    
      CALL XEWU3(LM1,MU1,LM2,MU2,LM3,MU3,I3,NEC,K0MAX,INDMAX,           
     1DEWU3,J1TA,J2TA,IEA,NCE1,NCE2,KIMAX1)                             
C     USING NEEDED STORAGE ONLY (ZERO OUT ENHANCEMENT)                  
      CALL CONMAT(I1,J1,LM1,MU1,K1MAX,L1,DON1,K1MAX)                    
      CALL CONMAT(I2,J2,LM2,MU2,K2MAX,L2,DON2,K2MAX)                    
      CALL CONMAT(I3,J3,LM3,MU3,K3MAX,L3,DON3,K3MAX)
C     USING NEEDED STORAGE ONLY (ZERO OUT ENHANCEMENT)                  
      CALL XWU3R3(I1,J1,LM1,MU1,K1MAX,L1,DON1,K1MAX,I2,J2,LM2,MU2,      
     1K2MAX,L2,DON2,K2MAX,I3,J3,LM3,MU3,K3MAX,L3,DON3,K3MAX,K0MAX,      
     2INDMAX,DEWU3,J1TA,J2TA,IEA,DWU3R3,K0MAX)
19    KOVER=0                                                           
      DO 20 K3=1,K3MAX                                                  
      DO 20 K2=1,K2MAX                                                  
      DO 20 K1=1,K1MAX                                                  
      DO 20 K0=1,K0MAX                                                  
         KOVER=KOVER+1 
         DVALU(K0,K1,K2,K3)=DWU3R3(KOVER)
   20    CONTINUE                                                       
      GOTO 50                                                           
   30 WRITE(6,40)LM1,MU1,LM2,MU2,LM3,MU3,L1,L2,L3,                      
     1K0MAX,K1MAX,K2MAX,K3MAX,ITYPE                                     
   40 FORMAT(1H0,13I4,2X,17H*ERROR* PAR-CHECK,2X,6HTYPE =,I2)           
      STOP                                                              
   50 RETURN                                                            
   60 RETURN 1                                                          
      END                                                               
