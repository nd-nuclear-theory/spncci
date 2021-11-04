/****************************************************************
  lgi_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT

  3/7/16 (aem,mac): Created.
  2/15/18 (aem): Update tests for Nsigma0ForNuclide and ReadLGISet
  10/24/21 (aem): Add test for lgi vector construction using lsu3shell
****************************************************************/
#include "lgi/lgi.h"

#include <cstdlib>

#include "am/halfint.h"
#include "am/halfint_fmt.h"
#include "fmt/format.h"
#include "utilities/nuclide.h"
#include "utilities/utilities.h"


int main(int argc, char **argv)
{
	
	u3::U3 omega(16,u3::SU3(4,0));
	HalfInt S=0, Sp=0, Sn=0;
	u3::U3S u3s(omega,S);
	u3shell::U3SPN u3spn(u3s,Sp,Sn);
	u3shell::U3SPN u3spn2(omega,S,Sp,Sn);
	lgi::LGI lgi(u3spn,0);
	std::cout<<fmt::format("{}",lgi)<<std::endl;
	std::cout<<lgi.Str()<<std::endl;

	// Read in lgi vector from file
	std::string spncci_root_dir=get_spncci_project_root_dir();
	std::string filename=fmt::format("{}/spncci/data/lgi_set/lgi_test_full.dat",spncci_root_dir);
	lgi::MultiplicityTaggedLGIVector lgi_vector;
	std::cout<<"lgi's from lgi_test.dat"<<std::endl;
	HalfInt Nsigma0=nuclide::Nsigma0ForNuclide({3,3});
	std::cout<<"Nsigma0 "<<Nsigma0<<std::endl;
	lgi::ReadLGISet(filename, Nsigma0,lgi_vector);
	for(const auto& [lgi,multiplicity] : lgi_vector)
		std::cout<<fmt::format("{}  {}",lgi,multiplicity)<<std::endl;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
