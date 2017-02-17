/****************************************************************
  spncci_basis_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/
#include "spncci/spncci_basis.h"

#include "cppformat/format.h"
#include "lgi/lgi.h"

#include "u3shell/relative_operator.h"
#include "spncci/spncci_branching_u3s.h"
#include "spncci/spncci_branching_u3lsj.h"

int main(int argc, char **argv)
{

  ////////////////////////////////////////////////////////////////
  // set up SpNCCI space
  ////////////////////////////////////////////////////////////////

  // read in LGIs
  std::string filename = "lgi_test.dat";  // test file in data/lgi_set/lgi_test.dat
  lgi::MultiplicityTaggedLGIVector multiplicity_tagged_lgi_vector;
  lgi::ReadLGISet(multiplicity_tagged_lgi_vector,filename);

  // diagnostic -- inspect LGI listing
  std::cout << "LGI set" << std::endl;
  for (int i=0; i<multiplicity_tagged_lgi_vector.size(); ++i)
    std::cout << i << " " << multiplicity_tagged_lgi_vector[i].Str() << std::endl;
  std::cout << "********************************" << std::endl;

  // generate SpNCCI space from LGIs
  HalfInt Nsigma_0 = HalfInt(11,1);
  int Nmax = 2;
  spncci::SpNCCISpace spncci_space;
  spncci::SigmaIrrepMap sigma_irrep_map;  // dictionary from sigma to branching
  spncci::NmaxTruncator truncator(Nsigma_0,Nmax);
  spncci::GenerateSpNCCISpace(multiplicity_tagged_lgi_vector,truncator,spncci_space,sigma_irrep_map);

  // diagnostic -- inspect irrep families
  if(true)
  {
    std::cout << "SpNCCI space" << std::endl;
    for (const spncci::SpNCCIIrrepFamily& spncci_irrep_family : spncci_space)
      {
        std::cout << "Irrep family" << std::endl;
        std::cout << "  " << spncci_irrep_family.Str() << std::endl;
        std::cout << "  " << fmt::format("Multiplicity gamma_max {}",spncci_irrep_family.gamma_max()) << std::endl;
        std::cout << "  " << "Irrep contents" << std::endl;
        std::cout << spncci_irrep_family.Sp3RSpace().DebugStr();
        std::cout << std::endl;
      }
    std::cout << "********************************" << std::endl;
  }

  ////////////////////////////////////////////////////////////////
  // count dimensions
  ////////////////////////////////////////////////////////////////

  if(true)
  {
    std::cout << fmt::format("  Irrep families {}",spncci_space.size()) << std::endl;
    std::cout << fmt::format("  TotalU3Subspaces {}",spncci::TotalU3Subspaces(spncci_space)) << std::endl;
    std::cout << fmt::format("  TotalDimensionU3 {}",spncci::TotalDimensionU3S(spncci_space)) << std::endl;
    std::cout << fmt::format("  TotalDimensionU3LS {}",spncci::TotalDimensionU3LS(spncci_space)) << std::endl;
    std::cout << "TotalDimensionU3LSJConstrained ";
    for (HalfInt J=0; J<10; ++J)
      std::cout << J << " " << spncci::TotalDimensionU3LSJConstrained(spncci_space,J) << "    ";
    std::cout << std::endl;
    std::cout << fmt::format("  TotalDimensionU3LSJAll {}",spncci::TotalDimensionU3LSJAll(spncci_space)) << std::endl;
  }

} //main
