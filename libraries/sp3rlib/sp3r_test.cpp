/****************************************************************
  sp3r_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/10/16 (aem,mac): Created.

****************************************************************/

#include "sp3rlib/sp3r.h"
#include "cppformat/format.h"
#include <iostream>
#include "sp3rlib/u3coef.h"
int main(int argc, char **argv)
{
  u3::U3CoefInit();

  ////////////////////////////////////////////////////////////////
  // Sp(3,R) irrep construction test
  ////////////////////////////////////////////////////////////////

  // 16(2,1)  Nn_max=6

  u3::U3 sigma = u3::U3(16,u3::SU3(2,1));
  int Nn_max = 6;

  // examine raising polynomials
  //
  // We generate and print them here for testing purposes only.
  // Normally sp3r::RaisingPolynomialLabels is called directly by the
  // Sp3RSpace constructor.

  std::vector<u3::U3> polynomial_labels = sp3r::RaisingPolynomialLabels(Nn_max);
  for (auto it = polynomial_labels.begin(); it != polynomial_labels.end(); ++it)
    std::cout << " label " << (*it).Str() << std::endl;
    
  // examine Sp3RSpace object
  sp3r::Sp3RSpace irrep(sigma,Nn_max);
  std::cout << irrep.DebugStr();
  std::cout<<irrep.size()<<std::endl;

  std::cout<<"Restricted irrep"<<std::endl;
  HalfInt Nsigma(9,2);
  u3::U3 sigma1 = u3::U3(Nsigma,u3::SU3(0,0));
  bool restrict_sp3r_to_u3_branching=true;
  sp3r::Sp3RSpace irrep1(sigma1,Nn_max,restrict_sp3r_to_u3_branching);
  std::cout << irrep1.DebugStr();
  std::cout<<irrep1.size()<<std::endl;


  std::cout<<"Bcoef cache check"<<std::endl;
  Nn_max=8;
  sp3r::BCoefCache cache;
  sp3r::GenerateBCoefCache(cache, Nn_max);
  u3::U3 n1,n2,n3;
  int rho;
  
  for(auto it=cache.begin(); it!=cache.end(); ++it)
    {
      std::tie(n1,n2,n3,rho)=it->first;
      double coef_flip=cache[sp3r::BCoefLabels(n2,n1,n3,rho)];
      double coef=it->second;
      if((fabs(coef_flip-coef)>1e-8)||(fabs(coef)<1e-8))
        std::cout<<fmt::format("B({} {} {} {})  {}  {}  {}",n1.Str(), n2.Str(),n3.Str(),rho,coef,coef_flip,fabs(coef_flip-coef))<<std::endl;
    }

} //main
