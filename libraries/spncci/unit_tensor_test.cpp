/****************************************************************
  unit_tensor_test.cpp

  Unit tensor algorithms
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/18/16 (aem,mac): Created.
****************************************************************/
#include <cstdio>
#include <fstream>
#include <sys/resource.h>

#include "cppformat/format.h"


#include "lgi/lgi.h"
#include "sp3rlib/u3coef.h"
#include "sp3rlib/vcs.h" 
#include "spncci/unit_tensor.h"
#include "spncci/spncci_branching_u3s.h"

int main(int argc, char **argv)
{
  u3::U3CoefInit();

  //unit tensor cache 
  u3::UCoefCache u_coef_cache;
  u3::PhiCoefCache phi_coef_cache;
  std::unordered_map<u3::U3,vcs::MatrixCache, boost::hash<u3::U3> >  k_matrix_map;
  spncci::SpNCCISpace sp_irrep_vector;
  // general diagnostics
  // stack size
  {
    struct rlimit rl;
    if (!getrlimit(RLIMIT_STACK, &rl))
      std::cout << "stack size (Mb) " << rl.rlim_cur/double(1024*1024) << std::endl;
  }

  // control caching status
  #ifdef USE_U_COEF_CACHE
  u3::g_u_cache_enabled = true;
  #else
  u3::g_u_cache_enabled = false;
  #endif
  std::cout << "u3::g_u_cache_enabled " << u3::g_u_cache_enabled << std::endl;

  // parse arguments
  if (argc<2)
    {
      std::cout << "Syntax: lgi_file_name Nmax" << std::endl;
      std::exit(1);
    }
  std::string filename(argv[1]);
  int Nmax = std::stoi(argv[2]);

  // For generating the lgi_vector, using Li-6 as example;
  HalfInt Nsigma_0 = HalfInt(11,1);
  int N1b=1;

  // Generate vector of LGI's from input file 
  lgi::MultiplicityTaggedLGIVector lgi_vector;
  lgi::ReadLGISet(lgi_vector,filename);
  // spncci::SpNCCISpace sp_irrep_vector;
  spncci::SigmaIrrepMap sigma_irrep_map;
  spncci::NmaxTruncator truncator(Nsigma_0,Nmax);
  std::cout<<"Generating irreps"<<std::endl;
  spncci::GenerateSpNCCISpace(lgi_vector,truncator,sp_irrep_vector,sigma_irrep_map);
  // Generate list of LGI's for which two-body operators will have non-zero matrix elements 
  std::cout<<"Sp irrep vector size "<<sp_irrep_vector.size()<<std::endl;
  // spncci::SpNCCIIrrepFamily sp_irrep=sp_irrep_vector[0].irrep;
  // sp3r::Sp3RSpace irrep=sp_irrep.Sp3RSpace();

  // enumerate pairs of families connected under two-body allowed spin
  // selection rules
  //
  // bra<=ket
  std::vector< std::pair<int,int> > spncci_irrep_family_pair_vector
    = GenerateSpNCCIIrrepFamilyPairs(sp_irrep_vector);

  //Generate list of sigma's and count of (sigma,S)
  // spncci::U3SCount sigma_S_count;
  std::unordered_set<u3::U3,boost::hash<u3::U3> >sigma_set;
  for(int l=0; l<sp_irrep_vector.size(); l++)
      sigma_set.insert(sp_irrep_vector[l].sigma());

  /////////////////////////////////////////////////////////////////////////////////////  
  // Generate Kmatrices 
  std::cout<<"Generating K matrices"<<std::endl;
  for( const auto& s : sigma_set)
    {
      vcs::MatrixCache K_map;
      // int Nex=int(s.N()-Nsigma_0);
      // vcs::GenerateKMatricesOpenMP(sigma_irrep_map[s], Nmax-Nex, K_map);
      vcs::GenerateKMatrices(sigma_irrep_map[s], K_map);

      k_matrix_map[s]=K_map;
    }

  // generate map that stores unit tensor labels keyed by N0
  std::map< int, std::vector<u3shell::RelativeUnitTensorLabelsU3ST> > unit_sym_map;
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax+2*N1b, unit_sym_map, 0, 0, false);

  // for(auto it=unit_sym_map.begin(); it!=unit_sym_map.end(); ++it)
  //   std::cout<<"ja "<<unit_sym_map[it->first].size()<<std::endl;
  //   // for(auto tensor : unit_sym_map[it->first])
  //   //   {
  //   //     std::cout<< it->first <<"  "<<tensor.Str()<<std::endl;
  //   //   }

  // spncci::GenerateUnitTensors(Nmax,unit_sym_map);

  //initializing map that will store map containing unit tensors.  Outer map is keyed SpNCCIIrrepFamily pair. 
  // inner map is keyed by unit tensor matrix element labels of type UnitTensorRME
  // SpNCCIIrrepFamily pair -> UnitTensorRME -> Matrix of reduced unit tensor matrix elements for v'v subsector
  std:: map< 
    std::pair<int,int>, 
    std::map<std::pair<int,int>,spncci::UnitTensorSectorsCache >
    > sp_irrep_unit_tensor_rme_map;
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Filling out sp_irrep_unit_tensor_rme_map 
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  
  u3::U3 sigma_old(0,0,0), sigmap_old(0,0,0);
  HalfInt Sp_old, S_old;
  for (int i=0; i<spncci_irrep_family_pair_vector.size(); i++)
    {
      std::pair<int,int> sp_irrep_pair=spncci_irrep_family_pair_vector[i];
      spncci::SpNCCIIrrepFamily sp_irrepp=sp_irrep_vector[sp_irrep_pair.first];
      spncci::SpNCCIIrrepFamily sp_irrep=sp_irrep_vector[sp_irrep_pair.second];

      int dimp=sp_irrep_vector[sp_irrep_pair.first].gamma_max();
      int dim=sp_irrep_vector[sp_irrep_pair.second].gamma_max();
      std::cout <<"SpNCCIIrrepFamily pair"<< sp_irrepp.Str()<<"  "<<sp_irrep.Str()<<std::endl;

      u3::U3 sigmap=sp_irrepp.sigma();  		
      u3::U3 sigma=sp_irrep.sigma();

      spncci::UnitTensorSectorsCache temp_unit_map;

      // operator boson number between two lgi's
      int N0=int(sigmap.N()-sigma.N());
      //////////////////////////////////////////////////////////////////////////////////////////////
      // Initializing the unit_tensor_rme_map with LGI rme's 
      //////////////////////////////////////////////////////////////////////////////////////////////
      // std::pair<int,int> N0_pair(sp_irrepp.Nex(),sp_irrep.Nex());
      // std::cout<<"Populating LGI rmes"<<std::endl;
      Eigen::MatrixXd temp_matrix=Eigen::MatrixXd::Constant(dimp,dim,1); 
      std::pair<int,int> N0_pair(0,0);
      for (int j=0; j<unit_sym_map[N0].size(); j++)
        {
          u3shell::RelativeUnitTensorLabelsU3ST unit_tensor=unit_sym_map[N0][j];
          u3::SU3 x0=unit_tensor.x0();
          HalfInt S0=unit_tensor.S0();
          int rp=unit_tensor.bra().eta();
          int r=unit_tensor.ket().eta();
          
          int rho0_max=u3::OuterMultiplicity(sigma.SU3(),x0, sigmap.SU3());
          for (int rho0=1; rho0<=rho0_max; rho0++)
            {
              if (rp<=(2*N1b+sp_irrepp.Nex()) && r<=(2*N1b+sp_irrep.Nex()))
                {
                  //std::cout<<unit_tensor.Str()<<std::endl;
                  // std::cout<<spncci::UnitTensorU3Sector(sigmap,sigma,unit_tensor,rho0).Str()<<std::endl;
                  // std::cout<<temp_matrix<<std::endl;    
                  temp_unit_map[spncci::UnitTensorU3Sector(sigmap,sigma,unit_tensor,rho0)]=temp_matrix;	
                }
            }
        }
      if(temp_unit_map.size())
        sp_irrep_unit_tensor_rme_map[sp_irrep_pair][N0_pair]=temp_unit_map;
      // std::cout<<"number of sp_irrep sectors "<<temp_unit_map.size()<<std::endl;;
      else
        continue;
      //////////////////////////////////////////////////////////////////////////////////////////////
      // Generating the rme's of the unit tensor for each SpNCCIIrrepFamily
      std::cout<<"Generating unit tensor sectors"<<std::endl;
      spncci::GenerateUnitTensorMatrix(
        N1b,Nmax,sp_irrep_pair,sp_irrep_vector,
        u_coef_cache,phi_coef_cache,k_matrix_map,
        unit_sym_map,
        sp_irrep_unit_tensor_rme_map);


      // for(auto it=sp_irrep_unit_tensor_rme_map[sp_irrep_pair].begin(); it!=sp_irrep_unit_tensor_rme_map[sp_irrep_pair].end(); ++it)
      //   {
      //     std::cout<<"Sp irrep pair "<<sp_irrep_pair.first<<"  "<<sp_irrep_pair.second<<std::endl;
      //     const std::pair<int,int>& NNpair=it->first;
      //     std::cout<<"Nnp and N "<<NNpair.first<<"  "<<NNpair.second<<std::endl;
      //     const spncci::UnitTensorSectorsCache& cache=it->second;
      //     for(auto it2=cache.begin(); it2!=cache.end(); ++it2)
      //       {
      //         std::cout<<it2->first.Str()<<std::endl;
      //         std::cout<<it2->second<<std::endl;
      //       }
      //   }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////
  // Getting interaction and setting up sectors 
  //////////////////////////////////////////////////////////////////////////////////////////////
  std::cout<<"Getting interaction"<<std::endl;
  // std::tuple<u3shell::OperatorLabelsU3S,int,int> IndexedOperatorLabelsU3S;
  std::string interaction_file="/Users/annamccoy/projects/spncci/data/trel_SU3_Nmax06.dat";
  std::ifstream interaction_stream(interaction_file.c_str());
  assert(interaction_stream);
  
  u3shell::RelativeRMEsU3ST interaction_rme_cache;
  u3shell::ReadRelativeOperatorU3ST(interaction_stream, interaction_rme_cache);
  
  // Get list of operators for U3S sector construction
  // From RelativeOperators list, get IndexedOperatorsLabelsU3S for sector construction
  std::unordered_set<u3shell::IndexedOperatorLabelsU3S,boost::hash<u3shell::IndexedOperatorLabelsU3S>>
      temp_operator_u3s_list;
  for(auto it=interaction_rme_cache.begin(); it!=interaction_rme_cache.end(); ++it)
    {
      int kappa0,L0;
      u3shell::RelativeUnitTensorLabelsU3ST tensor_u3st;
      std::tie(tensor_u3st,kappa0,L0)=it->first;
      u3shell::OperatorLabelsU3S operator_labels_u3s(tensor_u3st.operator_labels());
      temp_operator_u3s_list.emplace(operator_labels_u3s,kappa0,L0);
      // std::cout<<"operator labels "<<operator_labels_u3s.Str()<<"  "<<kappa0<<"  "<<L0<<std::endl;
    }
  std::vector<u3shell::IndexedOperatorLabelsU3S> operator_u3s_list(temp_operator_u3s_list.size());
  for(auto tensor : temp_operator_u3s_list)
      operator_u3s_list.push_back(tensor);
  // Get U3S space 
  spncci::SpaceU3S space(sp_irrep_vector);
  // Storage for sectors, value gives sector index
  
  // spncci::SectorLabelsU3SCache u3s_sectors;
  std::vector<spncci::SectorLabelsU3S> u3s_sector_vector;

  spncci::GetSectorsU3S(space,operator_u3s_list,u3s_sector_vector);

  //////////////////////////////////////////////////////////////////////////////////////////////
  // Contracting
  //////////////////////////////////////////////////////////////////////////////////////////////
  std::cout<<"contracting"<<std::endl;
  std::unordered_map<
    u3shell::IndexedOperatorLabelsU3S,
    Eigen::MatrixXd,boost::hash<u3shell::IndexedOperatorLabelsU3S>
    > contract_rme_cache;


  // iterate over interaction get unit tensor,kappa0,L0
  // iterate over sectors, get i,j, omega'S', rho0, omegaS->Nn etc. 
  // get corresponding unit tensor sector 

  basis::MatrixVector matrix_vector(u3s_sector_vector.size());
  for(auto it=interaction_rme_cache.begin(); it!=interaction_rme_cache.end(); ++it)
    {
      // Extract labels 
      int kappa0,L0;
      u3shell::RelativeUnitTensorLabelsU3ST tensor_u3st;
      std::tie(tensor_u3st,kappa0,L0)=it->first;
      double interaction_rme=it->second;
      int rp=tensor_u3st.bra().eta();
      int r=tensor_u3st.ket().eta();
      if((r>Nmax)||(rp>Nmax))
        continue;
      // Iterate over U3 sectors
      for(int s=0; s<u3s_sector_vector.size(); ++s)
      {
        const spncci::SectorLabelsU3S& sector=u3s_sector_vector[s];
        bool allowed=sector.operator_labels()==u3shell::OperatorLabelsU3S(tensor_u3st.operator_labels());
        allowed&=sector.kappa0()==kappa0;
        allowed&=(sector.L0()==L0);
        if(not allowed)
            continue;
        //otherwise
        //get subspace labels
        const spncci::SubspaceU3S& ket_subspace=space.GetSubspace(sector.ket_index());
        const spncci::SubspaceU3S& bra_subspace=space.GetSubspace(sector.bra_index());
        u3::U3 omegap=bra_subspace.GetSubspaceLabels().U3();
        u3::U3 omega=ket_subspace.GetSubspaceLabels().U3();
        int rho0=sector.rho0();
        int sector_dim_bra=bra_subspace.sector_dim();
        int sector_dim_ket=ket_subspace.sector_dim();

        matrix_vector[s]=Eigen::MatrixXd::Zero(sector_dim_bra,sector_dim_ket);
        // std::cout<<"sector size "<<sector_dim_bra<<" "<<sector_dim_ket<<std::endl;
        // std::cout<<matrix_vector[s]<<std::endl;
        // iterate over lgi multiplicities
        for(int i=0; i<bra_subspace.size(); ++i)
          for(int j=0; j<ket_subspace.size(); ++j)
            {
              int dimp, dim, indexp,index,gammap,gamma;
              u3::U3 sigma,sigmap;
              std::tie(gammap,sigmap,dimp,indexp)=bra_subspace.GetStateLabels(i);
              std::tie(gamma,sigma,dim,index)=ket_subspace.GetStateLabels(j);
              // std::cout<<fmt::format("   {} {} {} {}", indexp,index,dimp,dim)<<std::endl;
              std::pair<int,int> lgi_pair(gammap,gamma);
              std::pair<int,int> NnpNn(int(omegap.N()-sigmap.N()),int(omega.N()-sigma.N()));
              spncci::UnitTensorU3Sector unit_sector(omegap,omega,tensor_u3st,rho0);
              // std::cout<<"block"<<std::endl;
              // std::cout<<matrix_vector[s].block(indexp,index,dimp,dim)<<std::endl;
              // std::cout<<"subsector"<<std::endl;
              // std::cout<<interaction_rme*sp_irrep_unit_tensor_rme_map[lgi_pair][NnpNn][unit_sector]<<std::endl;
              
              // std::cout<<unit_sector.Str()<<"  "<<gammap<<"  "<<gamma<<std::endl;
              spncci::UnitTensorSectorsCache& cache=sp_irrep_unit_tensor_rme_map[lgi_pair][NnpNn];
              // std::cout<<"count  "<<cache.count(unit_sector)<<std::endl;
              if(not cache.count(unit_sector))
                for(auto t=cache.begin(); t!=cache.end(); ++t)
                  std::cout<<"    "<<t->first.Str()<<std::endl;
              if(cache.count(unit_sector))
                matrix_vector[s].block(indexp,index,dimp,dim)
                  +=cache[unit_sector];
            }
      }
    }
    std::cout<<"printing"<<std::endl;
  for(int s=0; s<matrix_vector.size();  ++s)
    {
      std::cout<<u3s_sector_vector[s].Str()<<std::endl;
      std::cout<<matrix_vector[s]<<std::endl;
    }







  // for(auto matrix : matrix_vector)
  //   std::cout<<matrix<<std::endl;

  // iterate over i,j
  //  iterate over NnpN
  //    iterate over omega'omega unit rho sectors
  //      get
  // if corresponding unit tensor in interaction
  //    Find correct U3S sector and accumulate 
  
  // for each component of interaction, 
  //  extract corresponding unit tensor rme
  //  multiply matrix by interaction rme
  //  sum over T, T0, Tp



  std::cout<<"all done"<<std::endl;

}
// end main 
