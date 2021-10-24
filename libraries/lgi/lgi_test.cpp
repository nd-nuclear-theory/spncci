/****************************************************************
  lgi_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT

  3/7/16 (aem,mac): Created.
  2/15/18 (aem): Update tests for Nsigma0ForNuclide and ReadLGISet
****************************************************************/
#include "lgi/lgi.h"

#include "am/halfint.h"
#include "am/halfint_fmt.h"
#include "fmt/format.h"
#include "utilities/nuclide.h"

int main(int argc, char **argv)
{
	u3::U3 omega(16,u3::SU3(4,0));
	HalfInt S=0, Sp=0, Sn=0;
	u3::U3S u3s(omega,S);
	u3shell::U3SPN u3spn(u3s,Sp,Sn);
	u3shell::U3SPN u3spn2(omega,S,Sp,Sn);
	lgi::LGI lgi(u3spn,0);
	std::cout<<lgi.Str()<<std::endl;
	std::string filename="data/lgi_set/lgi_test_full.dat";
	lgi::MultiplicityTaggedLGIVector lgi_vector;
	std::cout<<"lgi's from lgi_test.dat"<<std::endl;
	HalfInt Nsigma0=nuclide::Nsigma0ForNuclide({3,3});
	std::cout<<"Nsigma0 "<<Nsigma0<<std::endl;
	lgi::ReadLGISet(filename, Nsigma0,lgi_vector);
	// for(auto lgi : lgi_vector)
	// 	std::cout<<lgi.irrep.Str()<<"  "<<lgi.tag<<std::endl;

	//// Testing construction from lsu3shell
	std::cout<<"List of lgi's generated using lsu3shell basis constructors 6Li"<<std::endl;
	nuclide::NuclideType nuclide({3,3});
	int Nmax=2;
	lgi::MultiplicityTaggedLGIVector lgi_vector2 = lgi::get_lgi_vector(nuclide, Nsigma0,Nmax);
	// for(const auto& lgi_tagged : lgi_vector2)
	// 	std::cout<<lgi_tagged.irrep.Str()<<"  "<<lgi_tagged.tag<<std::endl;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//To compare lgi input file generated with old ordering, need to change ordering, 
	// so turn lgi_vector into map which does sorting and then write back to vector in new order.
	std::map<lgi::LGI,int> temp_map;
	for(const auto& lgi : lgi_vector)
		temp_map[lgi.irrep]=lgi.tag;

	lgi_vector.resize(0);
	for(auto lgi : temp_map)
		lgi_vector.emplace_back(lgi.first,lgi.second);

	for(int i=0; i<lgi_vector.size(); ++i)
		{
			if (not ((lgi_vector[i].irrep==lgi_vector2[i].irrep)&(lgi_vector[i].tag==lgi_vector2[i].tag)))
				std::cout<<"ERROR"<<std::endl;
			std::cout<<fmt::format("{}  {:6d}   {}  {:6d}", 
				lgi_vector[i].irrep.Str(),lgi_vector[i].tag,
				lgi_vector2[i].irrep.Str(),lgi_vector2[i].tag
			)<<std::endl;
		}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
