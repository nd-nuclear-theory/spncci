/****************************************************************
  su3lib_test.cpp

  Simple testbed program for C++ linkage from SU3LIB.
 
  Compilation and linkage:
    g++ -c su3lib_test.cpp
    g++ -o su3lib_test *.o -lgfortran


  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT

  3/9/16 (aem,mac): Created.
****************************************************************/



#include <cstdlib>
#include <iostream>



const size_t MAX_K = 9;

extern "C" 
{ 
  extern void wu3r3w_(const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, double[MAX_K][MAX_K][MAX_K][MAX_K]);
  extern void wru3optimized_(const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, double[], const int&);
  extern void wzu3optimized_(const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, double[], const int&);
  //extern void xewu3(const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, double[], int[], int[], int[], const int&, const int&, const int&);
  //extern void xwu3(const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, double[], const int&, const int&, double[], int[], int[], int[], int[], const int&, const int&, int[], const int&, const int&, const int&, const int&);
  extern void wu39lm_(const int&, const int& , const int&, const int&, const int& , const int& , const int& , const int&, const int&, const int&, const int&, const int&, const int& , const int& , const int& , const int&, const int&, const int&, double[], const int&);
  extern void blocks_(void);
}

void W_test()
{
    
  int x1lambda = 2;
  int x1mu = 0;
  int x2lambda = 1;
  int x2mu = 0;
  int x3lambda = 3;
  int x3mu = 0;

  int k1 = 1;
  int L1 = 0;
  int k2 = 1;
  int L2 = 1;
  int k3 = 1;
  int L3 = 1;

  int k1max = 1; // not used
  int k2max = 1;
  int k3max = 1;

  int r0 = 1;

  double w_array[MAX_K][MAX_K][MAX_K][MAX_K];

  wu3r3w_(
	  x1lambda, x1mu, x2lambda, x2mu, x3lambda, x3mu, 
	  L1, L2, L3, 
	  r0, 1, 1, 1, 
	  w_array
	  );
  double c = w_array[r0-1][k1-1][k2-1][k3-1];

  std::cout << c << std::endl;


}

int main(int argc, char **argv)
{
  blocks_();
  W_test();
}

