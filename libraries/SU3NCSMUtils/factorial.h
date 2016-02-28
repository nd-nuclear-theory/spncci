#ifndef factorial_h
#define factorial_h 
#include <cmath>

const int MAX_LOGFACT = 128; // when decomposing 2LaLb up to nmax = 8 ... I need at most LogFact(35) ==> 128 should be enough
const int MAX_SQRT = 1024; 

double logFact (int number);
double SQRT(int number); 
//inline double SQRT(int number) {return sqrt(number);}
void InitSqrtLogFactTables();
#endif 
