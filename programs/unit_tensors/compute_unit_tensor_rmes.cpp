/****************************************************************
  compute_unit_tensor_rmes.cpp

  compute relative unit tensor rmes in spncci basis using recurrance
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  12/3/16 (aem): Created.  Based on unit_tensor_test.cpp
  1/6/17 (aem) : Added contraction of unit tensors and interaction
****************************************************************/
#include <cstdio>
#include <ctime>
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


namespace spncci
{
  typedef std:: map< std::pair<int,int>, std::map<std::pair<int,int>,spncci::UnitTensorSectorsCache >
                      > LGIUnitTensorSectorCache;

  void 
  Contracting(int Nmax, int N1b,
              const std::vector<spncci::SectorLabelsU3S>& sector_labels_vector,
              const u3shell::RelativeRMEsU3ST& interaction_rme_cache,
              const spncci::SpaceU3S& space,
              LGIUnitTensorSectorCache& unit_tensor_sector_cache,
              basis::MatrixVector& matrix_vector
              )
  // Args:
  //  Nmax (input) : Basis truncation parameter
  //  N1b (input) : Basis single particle cutoff for Nmax=0
  //  sector_labels_vector (input) : vector of sector labels 
  //  interaction_rme_cache (input) : Container holding interaction rme's keyed
  //     by RelativeUnitTensorU3ST labels
  //  space (input) : space of omegaS subspaces
  //  unit_tensor_sector_cache (input) : nested container holding unit tensor rmes
  //  matrix_vector (output) : vector of U3S sectors indexed by U3S labels and kappa0,L0
  {
    // Contracting interaction matrix elements with unit tensor matrix elements
    std::cout<<"contracting"<<std::endl;
    // Initial sectors to zero matrices
     matrix_vector.resize(sector_labels_vector.size());
    for(int s=0; s<matrix_vector.size(); ++s)
    {
      const spncci::SectorLabelsU3S& sector=sector_labels_vector[s];
      const spncci::SubspaceU3S& ket_subspace=space.GetSubspace(sector.ket_index());
      const spncci::SubspaceU3S& bra_subspace=space.GetSubspace(sector.bra_index());
      int sector_dim_bra=bra_subspace.sector_dim();
      int sector_dim_ket=ket_subspace.sector_dim();
      matrix_vector[s]=Eigen::MatrixXd::Zero(sector_dim_bra,sector_dim_ket);
    }

    // iterate over interaction get unit tensor,kappa0,L0
    // iterate over sectors, get i,j, omega'S', rho0, omegaS->Nn etc. 
    // get corresponding unit tensor sector   
    for(auto it=interaction_rme_cache.begin(); it!=interaction_rme_cache.end(); ++it)
      {
        // Extract labels 
        int kappa0,L0;
        u3shell::RelativeUnitTensorLabelsU3ST tensor_u3st;
        std::tie(tensor_u3st,kappa0,L0)=it->first;
        double interaction_rme=it->second;
        //Check that unit tensor has rme between states in Nmax truncated basis
        int rp=tensor_u3st.bra().eta();
        int r=tensor_u3st.ket().eta();
        if((r>Nmax+2*N1b)||(rp>Nmax+2*N1b))
          continue;
        // Iterate over U3 sectors and find target sectors
        #pragma omp parallel for schedule(runtime)
        for(int s=0; s<sector_labels_vector.size(); ++s)
        {
          const spncci::SectorLabelsU3S& sector=sector_labels_vector[s];
          
          // Checking if sector is target sector
          bool allowed=sector.operator_labels()==u3shell::OperatorLabelsU3S(tensor_u3st.operator_labels());
          allowed&=sector.kappa0()==kappa0;
          allowed&=(sector.L0()==L0);
          if(not allowed)
              continue;
          
          //get subspace labels
          const spncci::SubspaceU3S& ket_subspace=space.GetSubspace(sector.ket_index());
          const spncci::SubspaceU3S& bra_subspace=space.GetSubspace(sector.bra_index());
          const u3::U3& omegap=bra_subspace.GetSubspaceLabels().U3();
          const u3::U3& omega=ket_subspace.GetSubspaceLabels().U3();
          int rho0=sector.rho0();
          // iterate through U3S subspaces, states correspond to lgi sets
          for(int i=0; i<bra_subspace.size(); ++i)
            for(int j=0; j<ket_subspace.size(); ++j)
              {
                // extracting subsector information
                // (dimp,dim)->size of subsector
                // (indexp,index)-> position of upper left corner of subsector in sector
                int dimp, dim, indexp,index,gammap,gamma;
                u3::U3 sigma,sigmap;
                std::tie(gammap,sigmap,dimp,indexp)=bra_subspace.GetStateLabels(i);
                std::tie(gamma,sigma,dim,index)=ket_subspace.GetStateLabels(j);

                //Keys for looking up subsector in unit tensor cache
                std::pair<int,int> lgi_pair(gammap,gamma);
                std::pair<int,int> NnpNn(int(omegap.N()-sigmap.N()),int(omega.N()-sigma.N()));              
                spncci::UnitTensorU3Sector unit_sector(omegap,omega,tensor_u3st,rho0);
                // Get cache containing unit tensor sector 
                if(not unit_tensor_sector_cache.count(lgi_pair))
                  continue;
                if(not unit_tensor_sector_cache[lgi_pair].count(NnpNn)) 
                  continue;
                spncci::UnitTensorSectorsCache& cache=unit_tensor_sector_cache[lgi_pair][NnpNn];
                
                //REMOVE : turn back on when finished debugging
                // Diagonistic: check if expected unit tensor sectors are found
                // if(not cache.count(unit_sector))
                  // for(auto t=cache.begin(); t!=cache.end(); ++t)
                    // std::cout<<t->first.Str()<<std::endl;
                #pragma omp critical
                {
                  if(cache.count(unit_sector))
                    matrix_vector[s].block(indexp,index,dimp,dim)+=interaction_rme*cache[unit_sector];
                }
              }
        }
      }
  }

}// end namespace

int main(int argc, char **argv)
{
      //REMOVE
      int testb=3;
      int testk=3;

  u3::U3CoefInit();
  //unit tensor cache 
  u3::UCoefCache u_coef_cache;
  u3::PhiCoefCache phi_coef_cache;

	u3::g_u_cache_enabled = true;
  double zero_threshold=1e-6;
  // For GenerateRelativeUnitTensors
  // should be consistant with lsu3shell tensors
  int T0=0;
  int J0=0;
  // parse arguments
  if (argc<8)
    {
      std::cout << "Syntax: A twice_Nsigma0 Nsigma0_ex_max N1B Nmax <basis filename> <nrel filename> <brel filename>" 
                << std::endl;
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
  spncci::LGIUnitTensorSectorCache unit_tensor_sector_cache;

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
  lgi::GenerateLGIExpansion(A,Nsigma_0,basis_table,space, is_brel,is_nrel,lgi_vector,lgi_expansion_matrix_vector);
  // Extract LGI labels and store in lgi_vector
  // lgi::GetLGILabels(Nsigma_0,space,lgi_expansion_matrix_vector, lgi_vector);
  int i=0;
  for(auto lgi_tagged : lgi_vector)
  {
    std::cout<<i<<"  "<<lgi_tagged.irrep.Str()<<"  "<<lgi_tagged.tag<<std::endl;
    ++i;
  }

  // Setting up the symplectic basis containers
  spncci::SigmaIrrepMap sigma_irrep_map;
  spncci::SpIrrepVector sp_irrep_vector;
  spncci::NmaxTruncator truncator(Nsigma_0,Nmax);
  // Generating sp3r irreps
  spncci::GenerateSp3RIrreps(lgi_vector,truncator,sp_irrep_vector,sigma_irrep_map);

  // Labels of relative unit tensors computed between LGI's
  std::vector<u3shell::RelativeUnitTensorLabelsU3ST> LGI_unit_tensor_labels;  
  std::cout<<"Nsigma0_ex_max "<<Nsigma0_ex_max<<std::endl;
  // -1 is all J0, T0=0 and false->don't restrict N0 to positive (temp)
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nsigma0_ex_max+2*N1b, LGI_unit_tensor_labels,-1,0,false);
  // For each operator, transform from lsu3shell basis to spncci basis
  std::cout<<"Number of unit tensors "<<LGI_unit_tensor_labels.size()<<std::endl;
  
  for(int i=0; i<LGI_unit_tensor_labels.size(); ++i)
    {
      //Get unit tensor labels
      u3shell::RelativeUnitTensorLabelsU3ST unit_tensor(LGI_unit_tensor_labels[i]);
      // std::cout<<"unit tensor "<<i<<"  "<<unit_tensor.Str()<<std::endl;
      u3shell::OperatorLabelsU3ST operator_labels(unit_tensor);
      // Generate sectors from labels
      u3shell::SectorsU3SPN sectors(space,operator_labels,false);
      // Read in lsu3shell rme's of unit tensor
      std::ifstream is_operator(fmt::format("relative_unit_{:06d}.rme",i));
      // If operator is not found, print name and continue;
      if(not is_operator)
        {
          std::cout<<fmt::format("relative_unit_{:06d}.rme not found",i)<<std::endl;
          continue;
        }

      // Read in operator from file
      // std::cout<<fmt::format("Transforming relative_unit_{:06d}.rme",i)<<std::endl;
      // std::cout<< LGI_unit_tensor_labels[i].Str()<<std::endl;
      // std::cout<<"  Reading in rmes"<<std::endl;
      basis::MatrixVector lsu3shell_operator_matrices(sectors.size());
      lsu3shell::ReadLSU3ShellRMEs(
        is_operator,operator_labels, basis_table,space, 
        sectors,lsu3shell_operator_matrices
        );

      // Basis transformation for each sector
      // std::cout<<"  Transforming"<<std::endl;

      // for(auto matrix : lsu3shell_operator_matrices)
      //   std::cout<<matrix<<std::endl<<std::endl;

      basis::MatrixVector spncci_operator_matrices;
      lgi::TransformOperatorToSpBasis(
        sectors,lgi_expansion_matrix_vector,
        lsu3shell_operator_matrices,spncci_operator_matrices
        );

      // std::cout<<"lgi expansions"<<std::endl;
      // for(auto matrix :lgi_expansion_matrix_vector)
      //   std::cout<<matrix<<std::endl<<std::endl;

      // Populate unit tensor map with sectors
      // std::cout<<"Populating lgi sectors"<<std::endl;
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
          if(not CheckIfZeroMatrix(spncci_operator_matrices[s],zero_threshold))
            unit_tensor_sector_cache[sp_irrep_pair][N0_pair][u3_sector]
                =spncci_operator_matrices[s];
          // std::cout<<i<<"  "<<j<<std::endl;
          // std::cout<<spncci_operator_matrices[s]<<std::endl;
        }
    }

  // for(auto it=unit_tensor_sector_cache.begin(); it!=unit_tensor_sector_cache.end(); ++it)
  //   {
  //     std::cout<<"Irreps "<<it->first.first<<" and "<<it->first.second<<std::endl;
  //   }

  // std::cout<<"traversing the map"<<std::endl;
  // auto it_ending=unit_tensor_sector_cache.begin();
  // for(int i=0; i<unit_tensor_sector_cache.size(); ++i) 
  //   it_ending++;

  // for(auto it=unit_tensor_sector_cache.begin(); it!=it_ending; ++it)
  //   {
  //     // if((it->first.first!=2)||(it->first.second!=5))
  //     //   continue;
  //     std::cout<<it->first.first<<"  "<<it->first.second<<std::endl;

  //     for(auto it2=it->second.begin(); it2!=it->second.end(); ++it2)
  //       {
  //         std::cout<<"  "<<it2->first.first<<"  "<<it2->first.second<<std::endl;
  //         for(auto it3=it2->second.begin(); it3!=it2->second.end();++ it3)
  //           {
  //             std::cout<<it3->first.Str()<<std::endl;
  //             // std::cout<<it3->second<<std::endl;
  //           }
  //       }
  //   }

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
      // vcs::GenerateKMatricesOpenMP(sigma_irrep_map[s], k_map);
      vcs::GenerateKMatrices(sigma_irrep_map[s], k_map);

      // std::cout<<"kmap "<<s.Str()<<std::endl;
      // for(auto it=k_map.begin(); it!=k_map.end(); ++it)
      //   std::cout<<it->first.Str()<<"  "<<it->second<<std::endl;
    }

  //////////////////////////////////////////////////////////////////////////////////////////
  // Computing unit tensor sectors
  //////////////////////////////////////////////////////////////////////////////////////////
  if(Nmax!=0)
  {
    std::map<int,std::vector<u3shell::RelativeUnitTensorLabelsU3ST>> unit_tensor_labels;
    // if true, unit tensors restricted to N0>=0
    // Nrel_max=Nmax+2N1b
    u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax+2*N1b, unit_tensor_labels, J0, T0, false);

    //REMOVE
    // std::map<std::pair<int,int>,std::vector<spncci::UnitTensorU3Sector>> unit_tensor_NpN_sector_map;
    // auto it_end=unit_tensor_sector_cache.begin();
    // for(int a=0; a<3; a++)
    //   ++it_end;

    // for(auto it=unit_tensor_sector_cache.begin(); it!=it_end; ++it)

    // std::tuple<accumulated time, std::vector<lgi_pairs>,std::vector<num non-zero unit tensors>>
    std::map<int,double> timing_map;
    std::map<int,std::vector<std::pair<int,int>>> lgi_distribution;
    std::map<int,std::vector<int>>  sector_count_map;
    std::map<int,std::vector<double>> individual_times;
    std::map<int,u3::UCoefCache> u_cache_map;
    std::pair<int,int> N0_pair(0,0);
    int num_nodes=20;
    int counter=0; 
    for(int n=0; n<num_nodes; ++n)
      timing_map[n]=0;

    for(auto it=unit_tensor_sector_cache.begin(); it!=unit_tensor_sector_cache.end(); ++it)
    {
      clock_t start_time=std::clock();

      // spncci::GenerateUnitTensorMatrix(
      //   N1b,Nmax,it->first,sp_irrep_vector,u_coef_cache,k_matrix_map,
      //   unit_tensor_labels,unit_tensor_sector_cache);
      int node=counter%num_nodes;
      assert(node<num_nodes);
      //Timing data
      spncci::GenerateUnitTensorMatrix(
        N1b,Nmax,it->first,sp_irrep_vector,u_cache_map[node], phi_coef_cache,k_matrix_map,
        unit_tensor_labels,unit_tensor_sector_cache);

      double duration=(std::clock()-start_time)/(double) CLOCKS_PER_SEC;
      int num_sector=unit_tensor_sector_cache[it->first][N0_pair].size();

      timing_map[node]+=duration;
      individual_times[node].push_back(duration);
      lgi_distribution[node].push_back(it->first);
      sector_count_map[node].push_back(num_sector);
      counter++;

    }
    std::cout<<"individual pairs"<<std::endl;
    for(int n=0; n<num_nodes; ++n)
      {
        std::cout<<"node"<<n<<std::endl;
        int i_stop=lgi_distribution[n].size();
        for(int i=0; i<i_stop; ++i)
          {
            std::cout<<lgi_distribution[n][i].first<<" "<<lgi_distribution[n][i].second<<"  "
            <<sector_count_map[n][i]<<"  "<<individual_times[n][i]<<std::endl;
          }
        std::cout<<"  "<<std::endl;

      }

    std::cout<<"summarizing"<<std::endl;
    for(int n=0; n<num_nodes; ++n)
      {
        std::cout<<"node "<<n<<std::endl;
        std::cout<<" time "<<timing_map[n]<<std::endl;
        std::cout<<" U cache size "<<u_cache_map[n].size()<<std::endl;
      } 
  }
  // REMOVE
  // auto it_stop=unit_tensor_sector_cache.begin();
  // int all=unit_tensor_sector_cache.size();
  // for(int a=0; a<all; ++a)
  //   ++it_stop;
  
  // std::cout<<"traversing the map"<<std::endl;
  // for(auto it=unit_tensor_sector_cache.begin(); it!=unit_tensor_sector_cache.end(); ++it)
  //   {
  //     // if((it->first.first!=2)||(it->first.second!=5))
  //     //   continue;
  //     std::cout<<it->first.first<<"  "<<it->first.second<<std::endl;
  //     for(auto it2=it->second.begin(); it2!=it->second.end(); ++it2)
  //       {
  //         std::cout<<"  "<<it2->first.first<<"  "<<it2->first.second<<std::endl;
  //         for(auto it3=it2->second.begin(); it3!=it2->second.end();++ it3)
  //           {
  //             std::cout<<it3->first.Str()<<std::endl;
  //             std::cout<<it3->second<<std::endl;
  //           }
  //       }
  // }

  //////////////////////////////////////////////////////////////////////////////////////////
  // // Getting interaction
  // //////////////////////////////////////////////////////////////////////////////////////////
  // //TODO make input 
  // //std::string interaction_file="/Users/annamccoy/projects/spncci/data/trel_SU3_Nmax06.dat";
  std::string interaction_file="id_SU3_Nmax16.dat";
  std::ifstream interaction_stream(interaction_file.c_str());
  assert(interaction_stream);
  
  u3shell::RelativeRMEsU3ST interaction_rme_cache;
  u3shell::ReadRelativeOperatorU3ST(interaction_stream, interaction_rme_cache);

  // for(auto it=interaction_rme_cache.begin(); it!=interaction_rme_cache.end(); ++it)
  //   std::cout<<it->second<<std::endl;

  std::vector<u3shell::IndexedOperatorLabelsU3S> operator_u3s_list;
  u3shell::GetInteractionTensorsU3S(interaction_rme_cache,operator_u3s_list);

  // Get U3S space 
  spncci::SpaceU3S u3s_space(sp_irrep_vector);
  // Storage for sectors, value gives sector index
  
  // spncci::SectorLabelsU3SCache u3s_sectors;
  std::vector<spncci::SectorLabelsU3S> u3s_sector_vector;
  // for(auto tensor : operator_u3s_list)
  //   {
  //     int kappa0, L0;
  //     u3shell::OperatorLabelsU3S tensor_u3st;
  //     std::tie(tensor_u3st,kappa0,L0)=tensor;
  //     std::cout<<tensor_u3st.Str()<<"  "<<kappa0<<"  "<<L0<<std::endl;
  //   }

  spncci::GetSectorsU3S(u3s_space,operator_u3s_list,u3s_sector_vector);

  //////////////////////////////////////////////////////////////////////////////////////////////
  // Contracting
  //////////////////////////////////////////////////////////////////////////////////////////////
  basis::MatrixVector matrix_vector;

  spncci::Contracting(Nmax, N1b,u3s_sector_vector,interaction_rme_cache,
                      u3s_space,unit_tensor_sector_cache, matrix_vector);

  // ZeroOutMatrix(matrix_vector,1e-6);
  // std::cout<<"printing"<<std::endl;
  // for(int i=0; i<u3s_space.size(); ++i)
  //   std::cout<<i<<"  "<<space.GetSubspace(i).GetSubspaceLabels().Str()<<std::endl;
  // for(int s=0; s<matrix_vector.size();  ++s)
  //   {
  //     if (not CheckIfZeroMatrix(matrix_vector[s], zero_threshold))
  //     {
  //       std::cout<<u3s_sector_vector[s].Str()<<std::endl;
  //       std::cout<<matrix_vector[s]<<std::endl;
  //     }
  //   }
}
