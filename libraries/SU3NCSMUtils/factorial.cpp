#include <SU3NCSMUtils/factorial.h>

#include <iostream>

using namespace std;

static double lnFact_table[MAX_LOGFACT];
static double sqrt_table[MAX_SQRT];

double logFact (int number) {return (number < MAX_LOGFACT) ?  lnFact_table[number] : logFact(number-1) + log(number);}

/*	
double logFact (int number) {return (number < MAX_LOGFACT) ?  lnFact_table[number] : logFact(number-1) + log(number);}
{
	static int maxnumber = 0;
	if (number > maxnumber)
	{
		maxnumber = number;
		cout << "maxnumber = " << maxnumber << endl;
	}
	return (number < MAX_LOGFACT) ?  lnFact_table[number] : logFact(number-1) + log(number);
}
*/	

double SQRT(int number) { return (number < MAX_SQRT) ? sqrt_table[number] : sqrt(number); }

void InitSqrtLogFactTables()
{
	lnFact_table[0] = lnFact_table[1] = 0;
	for (int i = 2; i < MAX_LOGFACT; i++) 
	{
		lnFact_table[i] = log(i) + lnFact_table[i-1];
	}

	for (int i = 0; i < MAX_SQRT; ++i)
	{
		sqrt_table[i] = sqrt(i);
	}
}
