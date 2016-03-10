#include <cstdlib>



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

int main(int argc, char **argv)
{
	blocks_();
}

//g++ -c su3lib_test.cpp
//g++ -o su3lib_test blocks.o su3lib_test.o -lgfortran