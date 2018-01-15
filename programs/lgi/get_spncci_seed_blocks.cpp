/****************************************************************
get_spncci_seed_blocks.cpp

Code which generates the lgi families and their corresponding 
seed blocks for the spncci reccurence starting from 
unit tensor rmes computed using LSU3Shell/SU3RME


Output: 

irrep_labels.dat: File containing a list of lgi family labels with 
  their corresponding index 
  
    i  twice_Nsigma lambda_sigma mu_sigma twice_Sp twice_Sn twice_S

For each pair of lgi families (i,j) with i>=j, there are two
files containing:
  1. seeds_00000i_00000j.rmes: All non-zero unit tensor seed blocks
      
      binary_format_code (1) and binary precision (4 or 8)
      num_rows, num_cols, rmes...

  2. operators_00000i_00000j.dat: unit tensor labels corresponding
      to each seed block
      
      lambda0, mu0, 2S0, 2T0, etap, 2Sp, 2Tp, eta, 2S, 2T, rho0
 

12/29/17 (aem): Created. 
****************************************************************/
#include <fstream>  
#include "cppformat/format.h"
#include "sp3rlib/u3coef.h"
#include "lgi/lgi_solver.h"
#include "mcutils/eigen.h"
#include "lgi/lgi_unit_tensors.h"
#include "spncci/spncci_basis.h"
#include "u3shell/relative_operator.h"

// Testing function 
namespace lgi{}

int main(int argc, char **argv)
{
  // zero tolerance 
  lgi::zero_tolerance=1e-6;
  
  // output mode
  lgi::binary_format_code = 1;
  lgi::binary_float_precision=8;

  // Basis parameters
  int Nmax=2;
  
  // nuclide 
  int Z=3;
  int N=3;
  
  std::array<int,2> nuclide; // proton and neutron numbers
  nuclide[0]=Z;
  nuclide[1]=N; 

  bool intrinsic=true;

  // su3rme output files
  std::string su3rme_filename_base="su3rme";
  std::string lsu3shell_basis_filename=su3rme_filename_base+"/lsu3shell_basis.dat"; // Will need to include path to file

  // Generate Nsigma0 from nuclei and type 
  HalfInt Nsigma0 = lgi::Nsigma0ForNuclide(nuclide,intrinsic);

  // Operator parameters
  std::string Brel_filename=su3rme_filename_base+"/Brel.rme";
  std::string Nrel_filename=su3rme_filename_base+"/Nrel.rme";

  // Unit tensor parameters
  int J0=-1;
  int T0=-1;

  int N1v=spncci::ValenceShellForNuclide(nuclide);

  //TODO
  std::string relative_unit_tensor_filename_template = su3rme_filename_base + "/relative_unit_{:06d}.rme";

  ////////////////////////////////////////////////////////////////
  // read lsu3shell basis
  ////////////////////////////////////////////////////////////////
  std::cout << "Read lsu3shell basis..." << std::endl;
  // read lsu3shell basis (regroup into U3SPN subspaces)
  lsu3shell::LSU3ShellBasisTable lsu3shell_basis_table;
  lsu3shell::U3SPNBasisLSU3Labels lsu3shell_basis_provenance;
  u3shell::SpaceU3SPN lsu3shell_space;
  lsu3shell::ReadLSU3ShellBasis(
      Nsigma0, lsu3shell_basis_filename,lsu3shell_basis_table,
      lsu3shell_basis_provenance,lsu3shell_space
    );


  std::cout<<lsu3shell_basis_filename<<std::endl;
  ////////////////////////////////////////////////////////////////
  // solve for LGIs
  ////////////////////////////////////////////////////////////////
  std::cout << "Solve for LGIs..." << std::endl;
  lgi::MultiplicityTaggedLGIVector lgi_families;
  lsu3shell::OperatorBlocks lgi_expansions;
  
  lgi::GetLGIExpansion(
      lsu3shell_space,lsu3shell_basis_table,
      Brel_filename,Nrel_filename,Z+N, Nsigma0,
      lgi_families, lgi_expansions
    );
  
  std::string lgi_filename="lgi_families.dat";
  lgi::WriteLGILabels(lgi_families, lgi_filename);
  ////////////////////////////////////////////////////////////////
  // Generate Seed blocks 
  ////////////////////////////////////////////////////////////////
  std::cout<<"Generate tensors "<<std::endl;
  // Get list of unit tensor labels between lgi's 
  // TODOL change to restrict N0 after testing
  bool restrict_positive_N0=false;
  std::vector<u3shell::RelativeUnitTensorLabelsU3ST> lgi_unit_tensor_labels;
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(
      Nmax, N1v,lgi_unit_tensor_labels,J0,T0,
      restrict_positive_N0
    );

  std::cout<<"comput seeds "<<std::endl;
  lgi::LGIGroupedSeedLabels lgi_grouped_seed_labels;
  std::vector<basis::MatrixVector> unit_tensor_spncci_matrices_array;
  bool restrict_seeds=false;
  lgi::ComputeUnitTensorSeedBlocks(
      lgi_unit_tensor_labels,relative_unit_tensor_filename_template,
      lsu3shell_space, lsu3shell_basis_table,lgi_expansions,
      lgi_grouped_seed_labels,unit_tensor_spncci_matrices_array, 
      restrict_seeds
    );

  std::cout<<"write seeds to file"<<std::endl;
  lgi::WriteSeedsToFile(lgi_grouped_seed_labels,unit_tensor_spncci_matrices_array);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Test code
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if(true)
    {
      for(auto it=lgi_grouped_seed_labels.begin(); it!=lgi_grouped_seed_labels.end(); ++it)
        { 
          //get operator file name  
          int bra_index, ket_index;
          std::tie(bra_index,ket_index)=it->first;
          std::string filename=fmt::format("operators_{:06d}_{:06d}.dat",bra_index,ket_index);
          
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
          basis::MatrixVector unit_tensor_seed_blocks(unit_tensor_labels.size());
          std::string seed_filename=fmt::format("seeds_{:06d}_{:06d}.rmes",bra_index,ket_index);
          
          // Read seeds from file
          int num_blocks=unit_tensor_labels.size();
          lgi::ReadSeedBlocks(seed_filename, num_blocks, unit_tensor_seed_blocks);
          
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
