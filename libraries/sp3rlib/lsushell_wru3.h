/****************************************************************
  lsushell_wru3.h

  SU(3) coupling coefficient port from Akiyama and Draayer su3lib,
  from LSU3shell.

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  4/22/16 (aem,mac): Created from code in T. Dytrych
  CWig6lmLookUpTable.cpp.

****************************************************************/

#ifndef LSUSHELL_WRU3_H_
#define LSUSHELL_WRU3_H_

namespace u3
{
  namespace lsu
  {

    void wru3(int LAM1, int MU1, int LAM2, int MU2, int LAM, int MU, int LAM3, int MU3, int LAM12, int MU12, int LAM23, int MU23, int KR0A, int KR0B, int KR0C, int KR0D, double* DRU3, int NABCD);

  } //namespace 
} //namespace 






#endif
