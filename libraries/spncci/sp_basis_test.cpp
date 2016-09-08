/****************************************************************
  sp_basis_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/sp_basis.h"
#include "lgi/lgi.h"

int main(int argc, char **argv)
{

  ////////////////////////////////////////////////////////////////
  // SpIrrep list
  ////////////////////////////////////////////////////////////////

  HalfInt Nsigma_0 = HalfInt(11,1);
  std::string filename = "libraries/spncci/lgi-3-3-2-fql-mini.dat";
  lgi::LGIVector lgi_vector;
  spncci::SpIrrepVector sp_irrep_vector;
  lgi::ReadLGISet(lgi_vector, filename);
  // lgi::ReadLGISet(lgi_vector,filename);

  // std::cout << "LGI vector" << std::endl;
  // for (int i=0; i<lgi_vector.size(); ++i)
  //   std::cout << i << " " << lgi_vector[i].Str() << std::endl;
  // std::cout << "********************************" << std::endl;

  ////////////////////////////////////////////////////////////////
  // build subspaces
  ////////////////////////////////////////////////////////////////

  int Nmax = 4;

  spncci::SigmaIrrepMap sigma_irrep_map;
  spncci::NmaxTruncator truncator(Nsigma_0,Nmax);
  spncci::GenerateSp3RIrreps(lgi_vector,truncator,sp_irrep_vector,sigma_irrep_map);

  std::cout << "SpIrrep vector reprise" << std::endl;
  for (int i=0; i<sp_irrep_vector.size(); ++i)
    std::cout << i << " " << sp_irrep_vector[i].DebugString();

  // examine irreps
  std::cout << "irreps (by sigma)" << std::endl;
  for (auto it = sigma_irrep_map.begin(); it != sigma_irrep_map.end(); ++it)
    {
      u3::U3 sigma = it->first;
      const sp3r::Sp3RSpace& irrep = it->second;

      std::cout << irrep.DebugString();
      std::cout << std::endl;
    }
  std::cout << "********************************" << std::endl;

  ////////////////////////////////////////////////////////////////
  // indicate dimensions
  ////////////////////////////////////////////////////////////////

  std::cout << "TotalU3Subspaces " << spncci::TotalU3Subspaces(sp_irrep_vector) << std::endl;
  std::cout << "TotalDimensionU3 " << spncci::TotalDimensionU3(sp_irrep_vector) << std::endl;
  std::cout << "TotalDimensionU3LS " << spncci::TotalDimensionU3LS(sp_irrep_vector) << std::endl;
  std::cout << "TotalDimensionU3LSJConstrained ";
  for (HalfInt J=0; J<10; ++J)
    std::cout << J << " " << spncci::TotalDimensionU3LSJConstrained(sp_irrep_vector,J) << "    ";
  std::cout << std::endl;
  std::cout << "TotalDimensionU3LSJAll " << spncci::TotalDimensionU3LSJAll(sp_irrep_vector) << std::endl;

} //main
