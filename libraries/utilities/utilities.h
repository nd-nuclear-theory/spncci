/****************************************************************
  utilities.h                       

  Define arithmetic shorthands.

  Mark A. Caprio
  University of Notre Dame

  Created by Mark A. Caprio on 2/17/11.   
  Drawing upon libmc/mcutils C code.
  2/23/11 (mac): Renamed from mc_arithmetic to arithmetic.
  3/9/16 (mac): Imported into spncci project as utilities.h.

****************************************************************/

#ifndef UTILITIES_H_
#define UTILITIES_H_
#include "gsl/gsl_sf.h"

// ONLYIF(cond,x) evaluates and returns x only if cond is true
#define ONLYIF(cond,x) ( (cond) ? (x) : 0)    

// sqr(x) returns the arithmetic square of x by self-multiplication
//   Note: Use of inline template avoids double evaluation of x which
//   would occur in a macro implementation.

template <typename T>
inline
T sqr(const T& x) 
{
  return x*x;
}

template <typename T>
inline
int KroneckerDelta(const T& x, const T& y)
// Return Kronecker delta of variables x and y.
//
// That is, returns 1 if x==y, 0 otherwise.
{
  return int(x==y);  // Is int(true) guaranteed to be 1?
}



inline int Choose(int x, int y)
{
  int choose=0;
  if ((x>=y)&&(y>=0))
    {
      gsl_sf_result result;
      gsl_sf_choose_e(x,y,&result);
      choose=result.val;
    }
  return choose;
}

inline int Factorial(int x)
{
  return gsl_sf_fact(x);
}

#endif
