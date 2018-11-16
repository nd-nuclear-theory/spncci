/****************************************************************
  branching2_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/branching2.h"


void GetBabySpncciBasis(spncci::BabySpNCCISpace& baby_spncci_space)
{
	// read in LGIs
  std::string filename 
  ="/afs/crc.nd.edu/user/a/amccoy/research/codes/spncci/data/lgi_set/lgi_test.dat";
  // = "~/research/codes/spncci/data/lgi_set/lgi_test.dat";  // test file in data/lgi_set/lgi_test.dat
  lgi::MultiplicityTaggedLGIVector multiplicity_tagged_lgi_vector;
  HalfInt Nsigma0=lgi::Nsigma0ForNuclide({3,3});
  lgi::ReadLGISet(filename, Nsigma0, multiplicity_tagged_lgi_vector);

  // generate SpNCCI space from LGIs
  // HalfInt Nsigma_0 = HalfInt(11,1);
  int Nmax = 2;
  spncci::SpNCCISpace spncci_space;
  spncci::SigmaIrrepMap sigma_irrep_map;  // dictionary from sigma to branching
  spncci::NmaxTruncator truncator(Nsigma0,Nmax);
  spncci::GenerateSpNCCISpace(multiplicity_tagged_lgi_vector,truncator,spncci_space,sigma_irrep_map);

  ////////////////////////////////////////////////////////////////
  // construct flattened baby SpNCCI space
  ////////////////////////////////////////////////////////////////

  // put SpNCCI space into standard linearized container
  baby_spncci_space=spncci::BabySpNCCISpace(spncci_space);

  // diagnostic
  std::cout << "baby_spncci_space" << std::endl;
  for (int subspace_index=0; subspace_index<baby_spncci_space.size(); ++subspace_index)
    std::cout << baby_spncci_space.GetSubspace(subspace_index).DebugStr()
              << std::endl;
  std::cout << std::endl;


}

int main(int argc, char **argv)
{
	HalfInt J(1,1); 
	spncci::BabySpNCCISpace baby_spncci_space;
	GetBabySpncciBasis(baby_spncci_space);
	spncci::SpaceSpBasis spbasis(baby_spncci_space,J);
  std::cout<<"full dimension "<<spbasis.FullDimension()<<std::endl;
  std::cout<<"size "<<spbasis.size()<<std::endl;
	std::cout<<spbasis.DebugStr(true)<<std::endl;


}