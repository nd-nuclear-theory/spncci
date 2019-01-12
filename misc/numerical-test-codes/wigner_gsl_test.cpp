/****************************************************************
 wigner_gsl_test.cpp

 Checks orthogonality of wigner coefficients
****************************************************************/
#include "fmt/format.h"
#include "am/wigner_gsl.h"

int main(int argc, char **argv)
{
  HalfInt S_min=0, S_max=1, L_max=2;
  
  HalfInt L1,S1,L2,S2,L,S,J1,J2,J;
  for(HalfInt L1=0; L1<=L_max; ++L1)
    for(HalfInt L2=0; L2<=L_max; ++L2)
      for(HalfInt S1=S_min; S1<=S_max; ++S1)
        for(HalfInt S2=S_min; S2<=S_max; ++S2)
          for(HalfInt J1=abs(L1-S1); J1<=(L1+S1); ++J1)
            for(HalfInt J2=abs(L2-S2); J2<=(L2+S2); ++J2)
              for(HalfInt J=abs(J1-J2); J<=(J1+J2); ++J)
              {
                double coef1=0;
                double coef2=0;
                for(HalfInt L=abs(L1-L2); L<=(L1+L2); ++L)
                  for(HalfInt S=abs(S1-S2); S<=(S1+S2); ++S)
                    {
                      coef1+=am::Unitary9J(L1,S1,J1,L2,S2,J2,L,S,J)
                              *am::Unitary9J(L1,S1,J1,L2,S2,J2,L,S,J);
                      coef2+=am::Unitary9J(L1,L2,L,S1,S2,S,J1,J2,J)
                              *am::Unitary9J(L1,L2,L,S1,S2,S,J1,J2,J);
                    }
                if(fabs(coef1-1)>10e-10)
                  std::cout<<fmt::format("([{} {}]{} [{} {}]{}|[{} {}] {})   {}",L1,S1,J1,L2,S2,J2,L,S,J, coef1)<<std::endl;
                if(fabs(coef2-1)>10e-10)
                  std::cout<<fmt::format("([{} {}]{} [{} {}]{}|[{} {}] {})   {}",L1,L2,L,S1,S2,S,J1,J2,J, coef2)<<std::endl;

              }

}

