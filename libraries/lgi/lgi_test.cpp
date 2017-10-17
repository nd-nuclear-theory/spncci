/****************************************************************
  lgi_test.cpp
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/7/16 (aem,mac): Created.
****************************************************************/
#include "lgi/lgi.h"

#include "cppformat/format.h"

int main(int argc, char **argv)
{
	u3::U3 omega(16,u3::SU3(4,0));
	HalfInt S=0, Sp=0, Sn=0;
	u3::U3S u3s(omega,S);
	u3shell::U3SPN u3spn(u3s,Sp,Sn);
	u3shell::U3SPN u3spn2(omega,S,Sp,Sn);
	lgi::LGI lgi(u3spn,0);
	std::cout<<lgi.Str()<<std::endl;
	std::string filename="../../data/lgi_set/lgi_test.dat";
	lgi::MultiplicityTaggedLGIVector lgi_vector;
	std::cout<<"lgi's from lgi_test.dat"<<std::endl;
	ReadLGISet(lgi_vector, filename);
	for(auto lgi : lgi_vector)
		std::cout<<lgi.irrep.Str()<<"  "<<lgi.tag<<std::endl;
	


  std::cout << "Nsigma0ForNuclide" << std::endl;
  std::cout
    << fmt::format(
        "3He {} 6Li {}",
        lgi::Nsigma0ForNuclide({2,1}),
        lgi::Nsigma0ForNuclide({3,3})
      )
    << std::endl;

}
