/****************************************************************
  compute_unit_tensor_rmes.cpp

  compute relative unit tensor rmes in spncci basis using recurrance
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  12/3/16 (aem): Created.  Based on unit_tensor_test.cpp
****************************************************************/
#include <cstdio>
#include <fstream>
#include <sys/resource.h>
#include "cppformat/format.h"

#include "lgi/lgi.h"
#include "lgi/lgi_solver.h"
#include "lsu3shell/lsu3shell_basis.h"
#include "lsu3shell/lsu3shell_rme.h"
#include "sp3rlib/u3coef.h"
#include "sp3rlib/vcs.h" 
#include "spncci/unit_tensor.h"
#include "u3shell/relative_operator.h"
#include "u3shell/upcoupling.h"

int main(int argc, char **argv)
{
  u3::U3CoefInit();
  //unit tensor cache 
  u3::UCoefCache u_coef_cache;
	u3::g_u_cache_enabled = true;

  // For GenerateRelativeUnitTensors
  // should be consistant with lsu3shell tensors
  int T0=0;
  int J0=0;
  // parse arguments
  if (argc<8)
    {
      std::cout << "Syntax: A twice_Nsigma0 Nsigma0_ex_max N1B Nmax <basis filename> <nrel filename> <brel filename>" << std::endl;
      std::exit(1);
    }
  int A = std::stoi(argv[1]); 
  int twice_Nsigma0= std::stoi(argv[2]);
  int Nsigma0_ex_max=std::stoi(argv[3]);
  int N1b=std::stoi(argv[4]);
  int Nmax = std::stoi(argv[5]);
  std::string lsu3_filename = argv[6];
  std::string nrel_filename = argv[7];
  std::string brel_filename = argv[8];

  HalfInt Nsigma_0=HalfInt(twice_Nsigma0,2);
  //initializing map that will store map containing unit tensors {SpIrrep pair : {}  
  // inner map is keyed by unit tensor matrix element labels of type UnitTensorRME
  // SpIrrep pair -> NnpNn -> UnitTensorRME -> v'v subsector
  std:: map< 
    std::pair<int,int>, 
    std::map<std::pair<int,int>,spncci::UnitTensorSectorsCache >
    > unit_tensor_sector_cache;

  //////////////////////////////////////////////////////////////////////////////////////////
  // Get LGI's and populate sector map with LGI rme's
  //////////////////////////////////////////////////////////////////////////////////////////
  // Generating LGI matrix elements 
  lsu3shell::LSU3BasisTable basis_table;
  lsu3shell::U3SPNBasisLSU3Labels basis_provenance;
  u3shell::SpaceU3SPN space;
  // Read in lsu3shell basis 
  lsu3shell::ReadLSU3Basis(Nsigma_0,lsu3_filename, basis_table, basis_provenance, space);

  //Get LGI expansion by solving for null space of Brel and Ncm
  std::ifstream is_brel(brel_filename.c_str());
  std::ifstream is_nrel(nrel_filename.c_str());
  lgi::LGIVector lgi_vector;
  basis::MatrixVector lgi_expansion_matrix_vector(space.size());
  lgi::GenerateLGIExpansion(A,basis_table,space, is_brel,is_nrel,lgi_expansion_matrix_vector);
  // Extract LGI labels and store in lgi_vector
  lgi::GetLGILabels(Nsigma_0,space,lgi_expansion_matrix_vector, lgi_vector);
  for(auto lgi_tagged : lgi_vector)
    std::cout<<lgi_tagged.irrep.Str()<<"  "<<lgi_tagged.tag<<std::endl;

  // Setting up the symplectic basis containers
  spncci::SigmaIrrepMap sigma_irrep_map;
  spncci::SpIrrepVector sp_irrep_vector;
  spncci::NmaxTruncator truncator(Nsigma_0,Nmax);
  // Generating sp3r irreps
  spncci::GenerateSp3RIrreps(lgi_vector,truncator,sp_irrep_vector,sigma_irrep_map);

  // Labels of relative unit tensors computed between LGI's
  std::vector<u3shell::RelativeUnitTensorLabelsU3ST> LGI_unit_tensor_labels;  
  std::cout<<"Nsigma0_ex_max "<<Nsigma0_ex_max<<std::endl;
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nsigma0_ex_max, LGI_unit_tensor_labels,-1,0,true);
  // For each operator, transform from lsu3shell basis to spncci basis
  std::cout<<"Number of unit tensors "<<LGI_unit_tensor_labels.size()<<std::endl;
  
  for(int i=0; i<LGI_unit_tensor_labels.size(); ++i)
    {
      //Get unit tensor labels
      u3shell::RelativeUnitTensorLabelsU3ST unit_tensor(LGI_unit_tensor_labels[i]);
      u3shell::OperatorLabelsU3ST operator_labels(unit_tensor);
      // Generate sectors from labels
      u3shell::SectorsU3SPN sectors(space,operator_labels,false);
      // Read in lsu3shell rme's of unit tensor
      std::ifstream is_operator(fmt::format("relative_unit_{:06d}",i));
      // If operator is not found, print name and continue;
      if(not is_operator)
        {
          std::cout<<fmt::format("relative_unit_{:06d} not found",i)<<std::endl;
          continue;
        }
      // Read in operator from file
      std::cout<<fmt::format("Transforming relative_unit_{:06d}",i)<<std::endl;
      std::cout<< LGI_unit_tensor_labels[i].Str()<<std::endl;

      std::cout<<"  Reading in rmes"<<std::endl;
      basis::MatrixVector lsu3shell_operator_matrices(sectors.size());
      lsu3shell::ReadLSU3ShellRMEs(
        is_operator,operator_labels, basis_table,space, 
        sectors,lsu3shell_operator_matrices
        );
      // Basis transformation for each sector
      std::cout<<"  Transforming"<<std::endl;
      basis::MatrixVector spncci_operator_matrices;
      lgi::TransformOperatorToSpBasis(
        sectors,lgi_expansion_matrix_vector,
        lsu3shell_operator_matrices,spncci_operator_matrices
        );
      // Populate unit tensor map with sectors
      std::cout<<"Populating lgi sectors"<<std::endl;
      std::pair<int,int> N0_pair(0,0);
      for(int s=0; s<sectors.size(); ++s)
        {
          int i=sectors.GetSector(s).bra_subspace_index();
          int j=sectors.GetSector(s).ket_subspace_index();
          int rho0=sectors.GetSector(s).multiplicity_index();
          std::pair<int,int>sp_irrep_pair(i,j);
          u3::U3 sigmap=sp_irrep_vector[i].irrep.sigma();
          u3::U3 sigma=sp_irrep_vector[j].irrep.sigma();
          spncci::UnitTensorU3Sector u3_sector(sigmap,sigma,unit_tensor,rho0);
          unit_tensor_sector_cache[sp_irrep_pair][N0_pair][u3_sector]
              =lsu3shell_operator_matrices[s];
        }
    }

  for(auto it=unit_tensor_sector_cache.begin(); it!=unit_tensor_sector_cache.end(); ++it)
    {
      std::cout<<"Irreps "<<it->first.first<<" and "<<it->first.second<<std::endl;
    }

  //////////////////////////////////////////////////////////////////////////////////////////
  //Precomputing Kmatrices 
  //////////////////////////////////////////////////////////////////////////////////////////
  std::unordered_set<u3::U3,boost::hash<u3::U3> >sigma_set;
  for(int l=0; l<sp_irrep_vector.size(); l++)
    {
      sigma_set.insert(sp_irrep_vector[l].irrep.sigma());
    }
  std::unordered_map<u3::U3,vcs::MatrixCache, boost::hash<u3::U3>> k_matrix_map;
  for( const auto& s : sigma_set)
    {
      vcs::MatrixCache& k_map=k_matrix_map[s];
      // int Nex=int(s.N()-Nsigma_0);
      // vcs::GenerateKMatricesOpenMP(sigma_irrep_map[s], Nmax-Nex, k_map);
      vcs::GenerateKMatrices(sigma_irrep_map[s], k_map);
    }

  //////////////////////////////////////////////////////////////////////////////////////////
  // Computing unit tensor sectors
  //////////////////////////////////////////////////////////////////////////////////////////
  std::map<int,std::vector<u3shell::RelativeUnitTensorLabelsU3ST>> unit_tensor_labels;
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax, unit_tensor_labels, J0, T0, true);

  // std::map<std::pair<int,int>,std::vector<spncci::UnitTensorU3Sector>> unit_tensor_NpN_sector_map;
  for(auto it=unit_tensor_sector_cache.begin(); it!=unit_tensor_sector_cache.end(); ++it)
  {
    
    std::map<std::pair<int,int>,std::vector<spncci::UnitTensorU3Sector>> 
      unit_tensor_NpN_sector_map;

    // std::cout<<"Getting sector labes"<<std::endl;
    spncci::GenerateUnitTensorU3SectorLabels(
      N1b,Nmax,it->first,sp_irrep_vector,
      unit_tensor_labels,unit_tensor_NpN_sector_map);

    spncci::GenerateUnitTensorMatrix(
      N1b,Nmax,it->first,sp_irrep_vector,u_coef_cache,k_matrix_map,
      // unit_tensor_NpN_sector_map,
      unit_tensor_labels,
      unit_tensor_sector_cache);
  }

  //////////////////////////////////////////////////////////////////////////////////////////
  // Getting interaction
  //////////////////////////////////////////////////////////////////////////////////////////
  //TODO make input 
  std::string interaction_file="/Users/annamccoy/projects/spncci/data/trel_SU3_Nmax06.dat";
  std::ifstream interaction_stream(interaction_file.c_str());
  assert(interaction_stream);
  
  u3shell::RelativeRMEsU3ST interaction_rme_cache;
  u3shell::ReadRelativeOperatorU3ST(interaction_stream, interaction_rme_cache);


}
