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

double W(const u3::SU3& x1, int k1, int L1, const u3::SU3& x2, int k2, int L2, const u3::SU3& x3, int k3, int L3, int r0)
  {
    

    //int r_max=u3::OuterMultiplicity(x1,x2,x3);
    //int k1max=BranchingMultiplicity(x1, L1);
    //int k2max=BranchingMultiplicity(x2, L2);
    //int k3max=BranchingMultiplicity(x3, L3);
    //int dummy=1;
    // Dimension of array that is returned contained coefficients is 9. 
    double w_array[su3lib::MAX_K][su3lib::MAX_K][su3lib::MAX_K][su3lib::MAX_K];
    su3lib::wu3r3w_(x1.lambda, x1.mu, x2.lambda, x2.mu, x3.lambda, x3.mu, L1 , L2, L3, r0,1,1,1, w_array);
    return w_array[r0-1][k1-1][k2-1][k3-1];
}
//g++ -c su3lib_test.cpp
//g++ -o su3lib_test blocks.o su3lib_test.o -lgfortran