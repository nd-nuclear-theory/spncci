/****************************************************************
get_spncci_seed_blocks.cpp

Code which generates the lgi families and their corresponding 
seed blocks for the spncci reccurence starting from 
unit tensor rmes computed using LSU3Shell/SU3RME


Output: All files are written to directory "seeds"

"seeds/lgi_families.dat": File containing a list of lgi family labels with 
  their corresponding index 
  
    i  twice_Nsigma lambda_sigma mu_sigma twice_Sp twice_Sn twice_S

For each pair of lgi families (i,j) with i>=j, there are 6
files in directory seeds containing:
  1. "seeds/seeds_00000i_00000j.rmes": All non-zero unit tensor seed blocks
      
      binary_format_code (1) and binary precision (4 or 8)
      num_rows, num_cols, rmes...

  2. "seeds/operators_00000i_00000j.dat": unit tensor labels corresponding
      to each seed block
      
      lambda0, mu0, 2S0, 2T0, etap, 2Sp, 2Tp, eta, 2S, 2T, rho0

  3. "seeds/obseeds_00000i_00000j.rmes": All non-zero one-body unit tensor seed blocks
      
      binary_format_code (1) and binary precision (4 or 8)
      num_rows, num_cols, rmes...

  4. "seeds/oboperators_00000i_00000j.dat": one-body unit tensor labels corresponding
      to each seed block
      
      lambda0, mu0, 2S0, etap, eta, Tz, rho

  5. "seeds/tbseeds_00000i_00000j.rmes": All non-zero two-body density seed blocks

      binary_format_code (1) and binary precision (4 or 8)
      num_rows, num_cols, rmes...

  4. "seeds/tboperators_00000i_00000j.dat": two-body density labels corresponding
      to each seed block

      N1, N2, N3, N4, lmf, muf, SSf, lmi, mui, SSi, rho0, lm0, mu0, SS0, Tz, rho
    
"seeds/lgi_expansions.dat": File containing expansion of lgis in lsu3shell basis.
    For each lsu3shell subspace, 
        rows, cols, rmes


  Anna E. McCoy
  TRIUMF

  SPDX-License-Identifier: MIT

12/29/17 (aem): Created. 
****************************************************************/
#include <fstream>  

#include "fmt/format.h"
#include "lgi/lgi_solver.h"
#include "lgi/lgi_unit_tensors.h"
#include "mcutils/eigen.h"
#include "mcutils/parsing.h"
#include "sp3rlib/u3coef.h"
#include "spncci/spncci_basis.h"
#include "u3shell/relative_operator.h"


// Testing function 
namespace lgi{}//namespace

int main(int argc, char **argv)
{
  if(argc<3)
  {
    std::cout<<"Syntax: Z N Nmax "<<std::endl;
    std::cout<<" or "<<std::endl;
    std::cout<<"Syntax: Z N Nmax lgi_filename lgi_expansion_filename"<<std::endl;
    std::exit(EXIT_FAILURE);
  }

  // nuclide 
  int Z=std::stoi(argv[1]);
  int N=std::stoi(argv[2]);

  // Basis parameters
  int Nmax=std::stoi(argv[3]);
  

  // zero tolerance 
  lgi::zero_tolerance=1e-6;
  
  // output mode
  lgi::binary_format_code = 1;
  lgi::binary_float_precision=8;

  nuclide::NuclideType nuclide; // proton and neutron numbers
  nuclide[0]=Z;
  nuclide[1]=N; 

  bool intrinsic=true;

  // su3rme output files
  std::string su3rme_filename_base="lsu3shell_rme";
  std::string lsu3shell_basis_filename=su3rme_filename_base+"/lsu3shell_basis.dat"; // Will need to include path to file

  // Generate Nsigma0 from nuclei and type 
  HalfInt Nsigma0 = nuclide::Nsigma0ForNuclide(nuclide,intrinsic);


  // Unit tensor parameters
  int J0=-1;
  int T0=-1;

  int N1v=nuclide::ValenceShellForNuclide(nuclide);

  //TODO
  std::string relative_unit_tensor_filename_template = su3rme_filename_base + "/relative_unit_{:06d}.rme";

  ////////////////////////////////////////////////////////////////
  // read lsu3shell basis
  ////////////////////////////////////////////////////////////////
  // std::cout << "Read lsu3shell basis..." << std::endl;
  // read lsu3shell basis (regroup into U3SPN subspaces)
  lsu3shell::LSU3ShellBasisTable lsu3shell_basis_table;
  lsu3shell::U3SPNBasisLSU3Labels lsu3shell_basis_provenance;
  u3shell::SpaceU3SPN lsu3shell_space;
  lsu3shell::ReadLSU3ShellBasis(
      Nsigma0, lsu3shell_basis_filename,lsu3shell_basis_table,
      lsu3shell_basis_provenance,lsu3shell_space
    );

  ////////////////////////////////////////////////////////////////
  // TODO: If lgi and lgi expansion file names not provided, solve for LGI
  // std::cout << "Solve for LGIs..." << std::endl;
  ////////////////////////////////////////////////////////////////
  // Operator parameters
  std::string Brel_filename=su3rme_filename_base+"/Brel.rme";
  std::string Nrel_filename=su3rme_filename_base+"/Nrel.rme";


  lgi::MultiplicityTaggedLGIVector lgi_families;
  lsu3shell::OperatorBlocks lgi_expansions;
  std::vector<int> lsu3hsell_index_lookup_table;

  lgi::GetLGIExpansion(
      lsu3shell_space,lsu3shell_basis_table,
      Brel_filename,Nrel_filename,Z+N, Nsigma0,
      lgi_families, lgi_expansions,
      lsu3hsell_index_lookup_table
    );
  
  std::string lgi_filename="seeds/lgi_families.dat";
  lgi::WriteLGILabels(lgi_families, lgi_filename);

  // std::cout<<"write expansion to file "<<std::endl;
  std::string lgi_expansion_filename="seeds/lgi_expansions.dat";
  lgi::WriteLGIExpansions(lgi_expansion_filename,lgi_expansions);



  ////////////////////////////////////////////////////////////////
  // Generate Seed blocks 
  ////////////////////////////////////////////////////////////////
  // std::cout<<"Generate tensors "<<std::endl;
  // Get list of unit tensor labels between lgi's 
  // TODO: Can we restrict N0?
  bool restrict_positive_N0=false;
  std::vector<u3shell::RelativeUnitTensorLabelsU3ST> lgi_unit_tensor_labels;
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(
      Nmax, N1v,lgi_unit_tensor_labels,J0,T0,
      restrict_positive_N0
    );

  // std::cout<<"comput seeds "<<std::endl;
  lgi::LGIGroupedSeedLabels lgi_grouped_seed_labels;
  std::vector<basis::OperatorBlocks<double>> unit_tensor_spncci_matrices_array;
  bool restrict_seeds=false;
  lgi::ComputeUnitTensorSeedBlocks(
      lgi_unit_tensor_labels,relative_unit_tensor_filename_template,
      lsu3shell_space, lsu3shell_basis_table,lgi_expansions,
      lgi_grouped_seed_labels,unit_tensor_spncci_matrices_array, 
      lsu3hsell_index_lookup_table,restrict_seeds
    );

  // std::cout<<"write seeds to file"<<std::endl;
  lgi::WriteSeedsToFile(lgi_grouped_seed_labels,unit_tensor_spncci_matrices_array);

//************************************************ Added by J.H. ********************************************
  int N1vp, N1vn;
  if(Z<=2)
    {
      N1vp=0;
    }
  else if(Z<=8)
    {
      N1vp=1;
    }
  else if(Z<=20)
    {
      N1vp=2;
    }
  if(N<=2)
    {
      N1vn=0;
    }
  else if(N<=8)
    {
      N1vn=1;
    }
  else if(N<=20)
    {
      N1vn=2;
    }

  std::string one_body_unit_tensor_filename_template = su3rme_filename_base + "/a+{}ta{}_{}_{}_{}_{}.rme";

  std::vector<u3shell::OneBodyUnitTensorLabelsU3S> lgi_one_body_unit_tensor_labels;
  u3shell::GenerateOneBodyUnitTensorLabelsU3S(
      Nmax,N1vp,N1vn,lgi_one_body_unit_tensor_labels
    );

  lgi::LGIGroupedOneBodySeedLabels lgi_grouped_one_body_seed_labels;
  std::vector<basis::OperatorBlocks<double>> one_body_unit_tensor_spncci_matrices_array;
  lgi::ComputeOneBodyUnitTensorSeedBlocks(
      lgi_one_body_unit_tensor_labels,one_body_unit_tensor_filename_template,
      lsu3shell_space, lsu3shell_basis_table,lgi_expansions,
      lgi_grouped_one_body_seed_labels,one_body_unit_tensor_spncci_matrices_array,
      lsu3hsell_index_lookup_table,false
    );

  lgi::WriteOneBodySeedsToFile(lgi_grouped_one_body_seed_labels,one_body_unit_tensor_spncci_matrices_array);

  std::string two_body_density_filename_template = su3rme_filename_base + "/{}_{}_{}_{}--{}_{}_{}--{}_{}_{}--{}_{}_{}_{}_{}.rme";

  std::vector<u3shell::TwoBodyDensityLabels> lgi_two_body_density_labels;
  u3shell::GenerateTwoBodyDensityLabels(
      Nmax,N1vp,N1vn,lgi_two_body_density_labels
    );

  lgi::LGIGroupedTwoBodySeedLabels lgi_grouped_two_body_seed_labels;

  std::vector<basis::OperatorBlocks<double>> two_body_density_spncci_matrices_array;
  lgi::ComputeTwoBodyDensitySeedBlocks(
      lgi_two_body_density_labels,two_body_density_filename_template,
      lsu3shell_space, lsu3shell_basis_table,lgi_expansions,
      lgi_grouped_two_body_seed_labels,two_body_density_spncci_matrices_array,
      lsu3hsell_index_lookup_table,false
    );

  lgi::WriteTwoBodySeedsToFile(lgi_grouped_two_body_seed_labels,two_body_density_spncci_matrices_array);
//***********************************************************************************************************

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Test code
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if(false)
    {
      for(auto it=lgi_grouped_seed_labels.begin(); it!=lgi_grouped_seed_labels.end(); ++it)
        { 
          //get operator file name  
          int bra_index, ket_index;
          std::tie(bra_index,ket_index)=it->first;
          std::string filename=fmt::format("seeds/operators_{:06d}_{:06d}.dat",bra_index,ket_index);
          
          // read in unit tensor sector labels 
          std::vector<u3shell::RelativeUnitTensorLabelsU3ST> unit_tensor_labels;
          std::vector<int> rho0_labels;
          lgi::ReadUnitTensorLabels(filename,unit_tensor_labels,rho0_labels);

          // Check that unit tensor labels read back in correctly
          for(int i=0; i<unit_tensor_labels.size();  ++i)
            {
              auto& labels2=it->second[i];
              assert(unit_tensor_labels[i]==std::get<0>(labels2));
              assert(rho0_labels[i]==std::get<1>(labels2));
            }

          // Reading in seeds 
          basis::OperatorBlocks<double> unit_tensor_seed_blocks(unit_tensor_labels.size());
          std::string seed_filename=fmt::format("seeds/seeds_{:06d}_{:06d}.rmes",bra_index,ket_index);
          
          // Read seeds from file
          int num_blocks=unit_tensor_labels.size();
          lgi::ReadBlocks(seed_filename, num_blocks, unit_tensor_seed_blocks);
          
          // Setting up i/o test
          auto& test=it->second;
          for(int i=0; i<unit_tensor_labels.size(); ++i)
            {
              int index1=std::get<2>(test[i]);
              int index2=std::get<3>(test[i]);

              auto& block1=unit_tensor_seed_blocks[i];
              auto& block2=unit_tensor_spncci_matrices_array[index1][index2];

              if(not mcutils::IsZero(block1-block2,1e-6))
                std::cout<<block1<<std::endl<<block2<<std::endl;
              // FINISH COMPARISON AND THEN CREATE FUNCTION FOR SPNCCI
            }
        }
        std::cout<<"all done checking "<<std::endl;
      }
    
}
