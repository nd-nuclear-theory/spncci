/****************************************************************
  lsushell_wru3.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include <algorithm>


#include "sp3rlib/lsushell_wru3.h"

namespace u3
{
  namespace lsu
  {
  
    extern "C" { 
      //    #ifndef AIX
      extern void xewu3_(int&, int&, int&, int&, int&, int&, int&, int&, int&, int&, double[], int[], int[], int[], int&, int&, int&);
      extern void xwu3_(int&, int&, int&, int&, int&, int&, int&, int&, int&, double[], int&, int&, double[], int[], int[], int[], int[], int&, int&, int[], int&, int&, int&, int&);
      extern void dlut_(int&, int&, int[], double[], int&);
      extern void dbsr_(int&, double[], double[], double[], int&);
      extern double drr3_(int&, int&, int&, int&, int&, int&); 
      // #else
      // extern void xewu3(int&, int&, int&, int&, int&, int&, int&, int&, int&, int&, double[], int[], int[], int[], int&, int&, int&);
      // extern void xwu3(int&, int&, int&, int&, int&, int&, int&, int&, int&, double[], int&, int&, double[], int[], int[], int[], int[], int&, int&, int[], int&, int&, int&, int&);
      // extern void dlut(int&, int&, int[], double[], int&);
      // extern void dbsr(int&, double[], double[], double[], int&);
      // extern double drr3(int&, int&, int&, int&, int&, int&); 
      // #endif
    }

    inline int IDM(int LAM, int MU) 
    {
      return (LAM+1)*(MU+1)*(LAM+MU+2)/2;
    }

    inline int INDEX(int J1TD, int LAM1, int J1T, int J2TD, int LAM2, int J2T)
    {
      return 1 + J2TD *(J2TD + 1)*(3*J1TD + J2TD + 5)/6 + (J1TD + 1)*(LAM2 + J2TD - J2T)/2 + (LAM1 + J1TD - J1T)/2;
    }

    int X1(9), X2(9), X3(9), X4(9), X5(13244), X6(13244), X7(13244), X8(13244), X9(39732), X10(39732);
    int X11(39732), X12(39732), X13(42), X14(378), X15(9030), X16(27090), X17(42), X18(1764);
    int HW(1);


    void wru3(int LAM1, int MU1, int LAM2, int MU2, int LAM, int MU, int LAM3, int MU3, int LAM12, int MU12, int LAM23, int MU23, int KR0A, int KR0B, int KR0C, int KR0D, double* DRU3, int NABCD)
    {
      int J1TD, J1T, J2TD, J2T;
      int I1, I2, I3, I4, I5, I6,  IAQ, IBQ, ICQ, IE12, IE23;
      int INDA, INDB, INDC,  INDMAX;
      int IS, J12T, J12TD, J23T, J23TS, J2S, J2SB;
      int J3S, J3SB, J3SQ, J3T, J3TD,  JJ12T, JJ12TA;
      int JJ12TB, JTDMAX, JTDMIN;
      int IE3MAX, IES, IESJ3S, IESMAX;
      int NA, NAB, NABC, NNCD, IDQ, JD, KD, JDKD, KDID;
      int KA, KABC, KABCD, KAIA;
      int KB, KBCDQ, KBCQ, KBIB, KBQ, KC, KCDQ, KCIC;
      int KCQ, KDQ, KIMAX1, KIMAX2, NECA, NECB, NECC, NECD;

      int IA[X13]; 
      int INDMAT[X18], J2SMAX[X18], J2TMAX[X18];
      int J3SMAX[X17], J3TMAX[X17];
	
      int JXTA[X5], JYTA[X5], IEA[X5]; 
      int JXTB[X6], JYTB[X6], IEB[X6];  
      int JXTC[X7], JYTC[X7], IEC[X7]; 
      int JXTD[X8], JYTD[X8], IED[X8]; 
	
      double DEWU3A[X9], DEWU3B[X10], DEWU3C[X11], DEWU3D[X12]; 
      double DWU3[X16]; 
      double DA[X14], DB[X4], DT[X4]; 
	
      double D1, D2, DC; 
	
      NA = KR0A;
      NAB = NA*KR0B;
      NABC = NAB*KR0C;
      KIMAX1 = X9;
      KIMAX2 = X16;
      xewu3_(LAM1, MU1, LAM23, MU23, LAM, MU, HW, NECD, KR0D, INDMAX, DEWU3D, JXTD, JYTD, IED, X4, X8, KIMAX1);
      NNCD = NECD + 1;
      IDQ = (INDMAX - NNCD - 1)*KR0D;
      for (JD = 1; JD <= NNCD; JD++) // DO 15 JD=1,NNCD
        { 
          IDQ = IDQ + KR0D; 
          IA[JD-1] = JD; 
          KDQ = -X13; 
          for (KD = 1; KD <= KR0D; ++KD) // DO 15 KD=1,KR0D
            { 
              KDQ = KDQ + X13; 
              JDKD = JD + KDQ; 
              KDID = KD + IDQ; 
              DA[JDKD-1] = DEWU3D[KDID -1];
            }
        } 
      dlut_(NNCD, KR0D, IA, DA, X13);

      xewu3_(LAM2, MU2, LAM3, MU3, LAM23, MU23, HW, NECC, KR0C, INDMAX, DEWU3C, JXTC, JYTC, IEC, X3, X7, KIMAX1);
      xewu3_(LAM12, MU12, LAM3, MU3, LAM, MU, HW, NECB, KR0B, INDMAX, DEWU3B, JXTB, JYTB, IEB, X2, X6, KIMAX1);
      xewu3_(LAM12, MU12, MU2, LAM2, LAM1, MU1, HW, NECA, KR0A, INDMAX, DEWU3A, JXTA, JYTA, IEA, X1, X5, KIMAX1); 
      I1 = LAM + 2*MU;
      I2 = LAM12 + 2*MU12;
      I3 = 4*LAM12 + 2*MU12;
      I4 = 2*I2;
      I5 = 2*(LAM12 - MU12);
      JTDMIN = std::max(0, std::max(NECA - LAM2 - MU2, NECB - LAM3 - MU3));
      JTDMAX = std::min(NECA, std::min(NECB, LAM12 + MU12));
      D1 = (double)((LAM1 + 1)*IDM(LAM12,MU12))/(double)(IDM(LAM1,MU1));
      IE23 = LAM1 + 2*MU1 - I1;
      J23TS = LAM23 + NECD + 2;
      KDQ = -NABC;
      for (KD = 1; KD <= KR0D; ++KD) // DO 30 KD=1,KR0D
        { 
          KDQ = KDQ + NABC; 
          J23T = J23TS - 2*IA[KD-1];                                               
          D2 = sqrt(D1*(J23T + 1));                                       
          xwu3_(LAM2, MU2, LAM3, MU3, LAM23, MU23, IE23, J23T, NECC, DEWU3C, KR0C, INDMAX, DWU3, J2SMAX, J2TMAX, J3SMAX, J3TMAX, IESMAX, IE3MAX, INDMAT, X3, X17, X15, KIMAX2); 
          IE12 = 3*IESMAX - IE3MAX - I1; 
          J12TD = (IE12 + I2)/3;  
          for (IES = 1; IES <= IESMAX; ++IES) //  DO 30 IES=1,IESMAX                                                
            { 
              IE12 = IE12 - 3;
              J12TD = J12TD - 1; 
              if (J12TD < JTDMIN || J12TD > JTDMAX) // IF(J12TD.LT.JTDMIN)GOTO 30 IF(J12TD.GT.JTDMAX)GOTO 30
                {
                  continue;
                }
              JJ12TA = I3 + IE12; 
              IS = I4 - IE12; 
              if (IS < JJ12TA) // IF(IS.LT.JJ12TA)
                {
                  JJ12TA = IS;
                }
              JJ12TA = JJ12TA/3 + 1;
              IS = (I5 - IE12)/3;
              JJ12TB = JJ12TA - abs(IS); 
              I6 = 2*LAM1 + IS; 
              J2TD = NECA - J12TD; 
              J3TD = NECB - J12TD; 
              J3T = J3TMAX[IES-1] + 2;                                                 
              J3SB = J3SMAX[IES-1]; 
              J3SQ = -X17;                          
              for (J3S = 1; J3S <= J3SB; ++J3S) // DO 25 J3S=1,J3SB
                { 
                  J3SQ = J3SQ + X17;
                  IESJ3S = IES + J3SQ;
                  J3T = J3T - 2;
                  J2T = J2TMAX[IESJ3S-1]+2; 
                  INDC = (INDMAT[IESJ3S-1] - J2T)/2; 
                  J2SB = J2SMAX[IESJ3S-1]; 
                  for (J2S=1; J2S <= J2SB; ++J2S) // DO 25 J2S=1,J2SB
                    { 
                      J2T = J2T - 2;
                      ICQ = INDC*KR0C; 
                      INDC = INDC + 1;
                      for (JJ12T=1; JJ12T <= JJ12TB; JJ12T += 2) // DO 25 JJ12T=1,JJ12TB,2                                            
                        { 
                          J12T = JJ12TA - JJ12T; 

                          INDA = INDEX(J12TD, LAM12, J12T, J2TD, MU2, J2T); 
                          if (JXTA[INDA-1] < 0) 
                            {
                              continue; //GOTO 25
                            } 

                          INDB = INDEX(J12TD, LAM12, J12T, J3TD, LAM3, J3T);
                          if (JXTB[INDB-1] < 0)
                            {
                              continue;
                            } 
                          DC = D2*drr3_(LAM1, J2T, LAM, J3T, J12T, J23T); 
                          IS = J12T + I6; 
                          if (4*(IS/4) != IS)
                            {
                              DC = -DC;
                            } 
                          IAQ = (INDA-1)*KR0A; 
                          IBQ = (INDB-1)*KR0B; 
                          KCQ = -NAB;
                          for (KC = 1; KC <= KR0C; ++KC) // DO 20 KC=1,KR0C
                            { 
                              KCQ = KCQ + NAB; 
                              KCDQ = KCQ + KDQ; 
                              KCIC = KC + ICQ; 
                              KBQ = -NA;                                                           
                              for (KB = 1; KB <= KR0B; ++KB) // DO 20 KB=1,KR0B
                                { 
                                  KBQ = KBQ + NA; 
                                  KBCDQ = KBQ + KCDQ; 
                                  KBIB = KB + IBQ; 
                                  for (KA = 1; KA <= KR0A; ++KA) // DO 20 KA=1,KR0A
                                    { 
                                      KABCD = KA + KBCDQ; 
                                      KAIA = KA + IAQ; 
                                      DRU3[KABCD-1] += DC*DEWU3A[KAIA-1]*DEWU3B[KBIB-1]*DWU3[KCIC-1];
                                    }
                                }
                            }
                        } // JJ12T: 25 CONTINUE
                    } // J2S: 25 CONTINUE
                } // J3S: 25 CONTINUE
            } // IES: 30 CONTINUE
        } // KD: 30 CONTINUE 
	
      KCQ = -NAB; 
      for (KC = 1; KC <= KR0C; ++KC) // DO 55 KC=1,KR0C
        { 
          KCQ = KCQ + NAB;                                                       
          KBQ = -NA;                                                           
          for (KB = 1; KB <= KR0B; ++KB) // DO 55 KB=1,KR0B
            { 
              KBQ = KBQ + NA;
              KBCQ = KBQ + KCQ; 
              for (KA = 1; KA <= KR0A; ++KA) // DO 55 KA=1,KR0A                                                   
                { 
                  KABC = KA + KBCQ; 
                  KDQ = -NABC; 
                  for (KD = 1; KD <= KR0D; ++KD) // DO 35 KD=1,KR0D                                                   
                    { 
                      KDQ = KDQ + NABC; 
                      KABCD = KABC + KDQ; 
                      DB[KD-1] = DRU3[KABCD-1]; // 35: KD
                    }	

                  /*	
                        IF(KR0D.GT.1)GOTO 40
                        DB(1)=DB(1)/DA(1) 
                        GOTO 45    
                        40 CALL DBSR(KR0D,DA,DB,DT,X13)
                        45 KDQ=-NABC
                  */
                  if (KR0D > 1) 
                    { 
                      dbsr_(KR0D, DA, DB, DT, X13);
                    }
                  else
                    { 
                      DB[0] = DB[0]/DA[0];
                    }
                  KDQ = -NABC;
                  for (KD = 1; KD <= KR0D; ++KD) // DO 50 KD=1,KR0D
                    { 
                      KDQ = KDQ + NABC; 
                      KABCD = KABC + KDQ; 
                      DRU3[KABCD-1] = DB[KD-1]; // KD: 50
                    } 
                } // KA: 55 continue
            } // KB: 55 continue
        } // KC: 55 CONTINUE                      
    }

    void wzu3(int LAM1, int MU1, int LAM2, int MU2, int LAM, int MU, int LAM3, int MU3, int LAM12, int MU12, int LAM23, int MU23, int KR0A, int KR0B, int KR0C, int KR0D, double* DZU3, int NABCD)
    {
      int J1TD, J1T, J2TD, J2T;
      int I1, I2, I3, I4, I5, I6, I7,  IAQ, IBQ, ICQ, IE12, IE23;
      int IE1, IE3, J1TS, J3TDMA, JJ3TA, JJ3TB, IE1MAX, IE1TES, J1S, J1SQ, IESJ1S, JJ3T, IPH;
      int INDA, INDB, INDC,  INDMAX;
      int IS, J12T, J12TD, J23T, J23TS, J2S, J2SB;
      int J3S, J3SB, J3SQ, J3T, J3TD,  JJ12T, JJ12TA;
      int JJ12TB, JTDMAX, JTDMIN;
      int IE3MAX, IES, IESJ3S, IESMAX;
      int NA, NAB, NABC, NNCD, IDQ, JD, KD, JDKD, KDID;
      int KA, KABC, KABCD, KAIA;
      int KB, KBCDQ, KBCQ, KBIB, KBQ, KC, KCDQ, KCIC;
      int KCQ, KDQ, KIMAX1, KIMAX2, NECA, NECB, NECC, NECD;

      int IA[X13]; 
      int INDMAT[X18], J2SMAX[X18], J2TMAX[X18];
      int J3SMAX[X17], J3TMAX[X17];
	
      int JXTA[X5], JYTA[X5], IEA[X5]; 
      int JXTB[X6], JYTB[X6], IEB[X6];  
      int JXTC[X7], JYTC[X7], IEC[X7]; 
      int JXTD[X8], JYTD[X8], IED[X8]; 
      int J1SMAX[X17], J1TMAX[X17];
	
      double DEWU3A[X9], DEWU3B[X10], DEWU3C[X11], DEWU3D[X12]; 
      double DWU3[X16]; 
      double DA[X14], DB[X4], DT[X4]; 
	
      double D1, D2, DC; 
	
      NA = KR0A;
      NAB = NA*KR0B; 
      NABC = NAB*KR0C; 
      KIMAX1 = X9; 
      KIMAX2 = X16; 
      xewu3_(LAM23, MU23, LAM1, MU1, LAM, MU, HW, NECD, KR0D, INDMAX, DEWU3D, JXTD, JYTD, IED, X4, X8, KIMAX1);
      NNCD = NECD + 1; 
      IDQ = (INDMAX - NNCD - 1)*KR0D; 
      for (JD = 1; JD <= NNCD; ++JD) // DO 15 JD=1,NNCD                                                   
        { 
          IDQ = IDQ + KR0D; 
          IA[JD-1] = JD; 
          KDQ = -X13;
          for (KD = 1; KD <= KR0D; ++KD) // DO 15 KD=1,KR0D                                                   
            { 
              KDQ = KDQ + X13; 
              JDKD = JD + KDQ; 
              KDID = KD + IDQ; 
              DA[JDKD-1]=DEWU3D[KDID-1];
            }
        } 
      dlut_(NNCD, KR0D, IA, DA, X13);

      xewu3_(LAM2, MU2, LAM1, MU1, LAM12, MU12, HW, NECA, KR0A, INDMAX, DEWU3A, JXTA, JYTA, IEA, X1, X5, KIMAX1); 
      xewu3_(LAM2, MU2, LAM3, MU3, LAM23, MU23, HW, NECC, KR0C, INDMAX, DEWU3C, JXTC, JYTC, IEC, X3, X7, KIMAX1); 
      xewu3_(LAM12, MU12, LAM3, MU3, LAM, MU, HW, NECB, KR0B, INDMAX, DEWU3B, JXTB, JYTB, IEB, X2, X6, KIMAX1); 
      I1 = LAM + 2*MU; 
      I2 = LAM12 + 2*MU12; 
      I3 = 4*LAM12 + 2*MU12; 
      I4 = 2*I2; 
      I5 = 2*(LAM12 - MU12); 
      I6 = LAM3 + 2*MU3; 
      I7 = 2*(LAM3 - MU3); 
      IE1 = LAM23 + 2*MU23 - I1;                                               
      J1TS = LAM1 + NECD + 2;                                                  
      J3TDMA = std::min(NECB, NECC) + 1;                                          
      for (J3S = 1; J3S <= J3TDMA; ++J3S) // DO 35 J3S=1,J3TDMA
        { 
          J3TD = J3S - 1;
          IE3 = -I6 + 3*J3TD; 
          JJ3TA = LAM3 + J3TD; 
          IS = I6 - J3TD; 
          if (IS < JJ3TA)
            {
              JJ3TA=IS;
            } 
          JJ3TA = JJ3TA + 1; 
          IS = (I7 - IE3)/3; 
          JJ3TB = JJ3TA - abs(IS);
          J12TD = NECB - J3TD; 
          J2TD = NECC - J3TD; 
          IE12 = -I2 + 3*J12TD; 
          JJ12TA = I3 + IE12; 
          IS = I4 - IE12; 
          if (IS < JJ12TA)
            {
              JJ12TA=IS;
            }
          JJ12TA = JJ12TA/3 + 1;
          IS = (I5-IE12)/3;
          JJ12TB = JJ12TA - abs(IS);
          for (JJ12T = 1; JJ12T <= JJ12TB; JJ12T += 2) // DO 35 JJ12T=1,JJ12TB,2
            {
              J12T = JJ12TA - JJ12T;
              xwu3_(LAM2, MU2, LAM1, MU1, LAM12, MU12, IE12, J12T, NECA, DEWU3A, KR0A, INDMAX, DWU3, J2SMAX, J2TMAX, J1SMAX, J1TMAX, IESMAX, IE1MAX, INDMAT, X1, X17, X15, KIMAX2);
              IE1TES = IE1MAX - IE1;
              if (IE1TES < 0)
                {
                  continue;
                } 
              if (IE1TES > 3*(IESMAX-1))
                {
                  continue;
                }
              IES = (IE1 - IE1MAX + 3*IESMAX)/3; 
              KDQ = -NABC;
              for (KD = 1; KD <= KR0D; ++KD) // DO 30 KD=1,KR0D
                {
                  KDQ = KDQ + NABC;
                  J1T = J1TS - 2*IA[KD-1];
                  J1S = (J1TMAX[IES-1] + 2 - J1T)/2;
                  if (J1S < 1)
                    {
                      continue;
                    }
                  if (J1S > J1SMAX[IES-1])
                    {
                      continue;
                    }
                  J1SQ = X17*(J1S - 1);
                  IESJ1S = IES + J1SQ;
                  J2T = J2TMAX[IESJ1S-1] + 2;
                  INDA = (INDMAT[IESJ1S-1] - J2T)/2;
                  J2SB = J2SMAX[IESJ1S-1];
                  for (J2S = 1; J2S <= J2SB; ++J2S) // DO 25 J2S=1,J2SB                                                  
                    { 
                      J2T = J2T - 2; 
                      IAQ = INDA*KR0A; 
                      INDA = INDA + 1;
                      for (JJ3T = 1; JJ3T <= JJ3TB; JJ3T += 2) // DO 25 JJ3T=1,JJ3TB,2                                              
                        { 
                          J3T = JJ3TA - JJ3T; 
                          INDC = INDEX(J2TD, LAM2, J2T, J3TD, LAM3, J3T); 
                          if (JXTC[INDC-1] < 0) 
                            {
                              continue;
                            } 
                          INDB = INDEX(J12TD, LAM12, J12T, J3TD, LAM3, J3T); 
                          if (JXTB[INDB-1] < 0)
                            {
                              continue;
                            }
                          IPH = J2T + LAM + 3*(J12T + LAM23);
                          D2 = (double)(1 - abs(IPH)%4);
                          DC = D2*drr3_(J1T, J2T, LAM, J3T, J12T, LAM23);
                          DC = DC*sqrt((double)((J12T + 1)*(LAM23 + 1))); 
                          ICQ = (INDC-1)*KR0C; 
                          IBQ = (INDB-1)*KR0B;                                                 
                          KCQ = -NAB;
                          for (KC = 1; KC <= KR0C; ++KC) // DO 20 KC=1,KR0C
                            { 
                              KCQ = KCQ + NAB; 
                              KCDQ = KCQ + KDQ; 
                              KCIC = KC + ICQ; 
                              KBQ = -NA; 
                              for (KB = 1; KB <= KR0B; ++KB) // DO 20 KB=1,KR0B
                                { 
                                  KBQ = KBQ + NA; 
                                  KBCDQ = KBQ + KCDQ; 
                                  KBIB = KB + IBQ; 
                                  for (KA = 1; KA <= KR0A; ++KA)// DO 20 KA=1,KR0A
                                    {
                                      KABCD = KA + KBCDQ;
                                      KAIA = KA + IAQ;
                                      DZU3[KABCD-1] = DZU3[KABCD-1]+DC*DEWU3C[KCIC-1]*DEWU3B[KBIB-1]*DWU3[KAIA-1];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } 
	
      KCQ = -NAB;
      for (KC = 1; KC <= KR0C; ++KC) // DO 60 KC=1,KR0C                                                   
        { 
          KCQ = KCQ + NAB; 
          KBQ = -NA;
          for (KB = 1; KB <= KR0B; ++KB) // DO 60 KB=1,KR0B                                                   
            { 
              KBQ = KBQ + NA; 
              KBCQ = KBQ + KCQ; 
              for (KA = 1; KA <= KR0A; ++KA) // DO 60 KA=1,KR0A                                                   
                { 
                  KABC = KA + KBCQ; 
                  KDQ = -NABC;
                  for (KD = 1; KD <= KR0D; ++KD) // DO 40 KD=1,KR0D                                                   
                    { 
                      KDQ = KDQ + NABC; 
                      KABCD = KABC + KDQ; 
                      DB[KD-1] = DZU3[KABCD-1];
                    }
                  /*
                    IF(KR0D.GT.1)GOTO 45                                              
                    DB(1)=DB(1)/DA(1)                                                 
                    GOTO 50                                                           
                    45 CALL DBSR(KR0D,DA,DB,DT,X13)                                      
                  */
                  if (KR0D > 1)
                    {
                      dbsr_(KR0D, DA, DB, DT, X13);
                    }
                  else
                    {
                      DB[0] = DB[0]/DA[0];
                    }
                  KDQ = -NABC;
                  for (KD = 1; KD <= KR0D; ++KD) // DO 55 KD=1,KR0D                                                   
                    { 
                      KDQ = KDQ + NABC; 
                      KABCD = KABC + KDQ; 
                      DZU3[KABCD-1] = DB[KD-1];
                    }
                }
            }
        }
    }


  } // namespace 
} // namespace 
