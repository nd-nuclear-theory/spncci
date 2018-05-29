  
#include <iostream>
#include <fstream>

#include "SymEigsSolver.h"  // from spectra
#include "cppformat/format.h"

#include "lgi/lgi_solver.h"
#include "mcutils/eigen.h"
#include "sp3rlib/u3coef.h"
#include "spncci/spncci_common.h"


int main(int argc, char **argv)
{
  std::cout<<"generating lgi expansions"<<std::endl;
  ////////////////////////////////////////////////////////////////
  // initialization
  ////////////////////////////////////////////////////////////////
  
  // SU(3) caching
  u3::U3CoefInit();
  u3::UCoefCache u_coef_cache;
  u3::PhiCoefCache phi_coef_cache;
  u3::g_u_cache_enabled = true;

  // parameters for certain calculations
  // lgi::g_zero_tolerance = 1e-6;
  spncci::g_suppress_zero_sectors = true;

  // Eigen OpenMP multithreading mode
  Eigen::initParallel();

  ////////////////////////////////////////////////////////////////
  // read lsu3shell basis (regroup into U3SPN subspaces)
  ////////////////////////////////////////////////////////////////
  std::string lsu3shell_basis_filename; //TODO
  HalfInt Nsigma0;
  std::string Brel_filename;
  std::string Nrel_filename;
  int A; 


  std::cout << "Read lsu3shell basis..." << std::endl;
  lsu3shell::LSU3ShellBasisTable lsu3shell_basis_table;
  lsu3shell::U3SPNBasisLSU3Labels lsu3shell_basis_provenance;
  u3shell::SpaceU3SPN lsu3shell_space;
  lsu3shell::ReadLSU3ShellBasis(
      Nsigma0,lsu3shell_basis_filename,lsu3shell_basis_table,
      lsu3shell_basis_provenance,lsu3shell_space
    );

  ////////////////////////////////////////////////////////////////
  // solve for LGIs
  ////////////////////////////////////////////////////////////////
  std::cout << "Solve for LGIs..." << std::endl;

  lgi::MultiplicityTaggedLGIVector lgi_families;
  basis::OperatorBlocks<double> lgi_expansions;
  std::vector<int> lsu3shell_index_lookup_table;

  lgi::GetLGIExpansion(
      lsu3shell_space,lsu3shell_basis_table,
      Brel_filename,Nrel_filename,A,Nsigma0,
      lgi_families, lgi_expansions,lsu3shell_index_lookup_table
    );

  // diagnostics
  std::cout << fmt::format("  LGI families {}",lgi_families.size()) << std::endl;
  
  for(int i=0; i<lgi_families.size(); ++i)
    {
      int Nex;
      u3::U3 sigma;
      HalfInt Sp,Sn,S;
      auto& lgi_family=lgi_families[i];
      std::tie(Nex,sigma,Sp,Sn,S)=lgi_family.irrep.Key();
      int gamma_max=lgi_family.tag;
        std::cout
          <<Nex
          <<"  "<<TwiceValue(sigma.N())<<"  "<<sigma.SU3().lambda()<<"  "<<sigma.SU3().mu()
          <<"  "<<TwiceValue(Sp)<<"  "<<TwiceValue(Sn)<<"  "<<TwiceValue(S)
          <<"  "<<gamma_max
          <<std::endl;     

      const spncci::OperatorBlock& lgi_expansion=lgi_expansions[i];
      std::cout<<mcutils::FormatMatrix(lgi_expansion,"13.6e")<<std::endl;

    }

  // if (true)
  //   lgi::WriteLGILabels(lgi_families_truncated,std::cout);


}
