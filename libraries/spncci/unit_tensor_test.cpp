/****************************************************************
  unit_tensor_test.cpp

  Unit tensor algorithms
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/18/16 (aem,mac): Created.
****************************************************************/
#include <cstdio>
#include <sys/resource.h>

#include "lgi/lgi.h"
#include "sp3rlib/u3coef.h"
#include "sp3rlib/vcs.h" 
#include "spncci/unit_tensor.h"


int main(int argc, char **argv)
{
  u3::U3CoefInit();

  //unit tensor cache 
  u3::UCoefCache u_coef_cache;
  std::unordered_map<u3::U3,vcs::MatrixCache, boost::hash<u3::U3> >  k_matrix_map;
  spncci::SpIrrepVector sp_irrep_vector;
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
  lgi::LGIVector lgi_vector;
  lgi::ReadLGISet(lgi_vector,filename);
  // spncci::SpIrrepVector sp_irrep_vector;
  spncci::SigmaIrrepMap sigma_irrep_map;
  spncci::NmaxTruncator truncator(Nsigma_0,Nmax);
  spncci::GenerateSp3RIrreps(lgi_vector,truncator,sp_irrep_vector,sigma_irrep_map);
  // Generate list of LGI's for which two-body operators will have non-zero matrix elements 
  std::vector< std::pair<int,int> > sp_irrep_pair_vector
    =spncci::GenerateSpIrrepPairs(sp_irrep_vector);

  //Generate list of sigma's and count of (sigma,S)
  // spncci::U3SCount sigma_S_count;
  std::unordered_set<u3::U3,boost::hash<u3::U3> >sigma_set;
  for(int l=0; l<sp_irrep_vector.size(); l++)
    {
      sigma_set.insert(sp_irrep_vector[l].irrep.sigma());
  //     std::pair<u3::U3,HalfInt> count_key(sp_irrep_vector[l].sigma(), sp_irrep_vector[l].S());
  //     sigma_S_count[count_key]+=1;
    }

  // for(auto it=sigma_S_count.begin(); it!=sigma_S_count.end(); ++it)
  //   {
  //     u3::U3 s;
  //     HalfInt SS;
  //     std::tie(s,SS)=it->first;
  //     int count=it->second;
  //     std::cout<<"sigma S pairs "<<s.Str()<<"  "<<SS<<"  "<<count<<std::endl;
  //   }
  /////////////////////////////////////////////////////////////////////////////////////  
  // Generate Kmatrices 
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
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax, unit_sym_map, 0, 0, true);

  // spncci::GenerateUnitTensors(Nmax,unit_sym_map);

  //initializing map that will store map containing unit tensors.  Outer map is keyed SpIrrep pair. 
  // inner map is keyed by unit tensor matrix element labels of type UnitTensorRME
  // SpIrrep pair -> UnitTensorRME -> Matrix of reduced unit tensor matrix elements for v'v subsector
  std:: map< 
    std::pair<int,int>, 
    std::map<std::pair<int,int>,spncci::UnitTensorSectorsCache >
    > sp_irrep_unit_tensor_rme_map;
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Filling out sp_irrep_unit_tensor_rme_map 
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  u3::U3 sigma_old(0,0,0), sigmap_old(0,0,0);
  HalfInt Sp_old, S_old;
  for (int i=0; i<sp_irrep_pair_vector.size(); i++)
    {
      std::pair<int,int> sp_irrep_pair=sp_irrep_pair_vector[i];
      spncci::SpIrrep sp_irrepp=sp_irrep_vector[sp_irrep_pair.first].irrep;
      spncci::SpIrrep sp_irrep=sp_irrep_vector[sp_irrep_pair.second].irrep;

      std::cout <<"SpIrrep pair"<< sp_irrepp.Str()<<"  "<<sp_irrep.Str()<<std::endl;

      u3::U3 sigmap=sp_irrepp.sigma();  		
      u3::U3 sigma=sp_irrep.sigma();

      spncci::UnitTensorSectorsCache temp_unit_map;

      // operator boson number between two lgi's
      int N0=int(sigmap.N()-sigma.N());
      //////////////////////////////////////////////////////////////////////////////////////////////
      // Initializing the unit_tensor_rme_map with LGI rme's 
      //////////////////////////////////////////////////////////////////////////////////////////////
      // std::pair<int,int> N0_pair(sp_irrepp.Nex(),sp_irrep.Nex());
      std::pair<int,int> N0_pair(0,0);
      for (int j=0; j<unit_sym_map[N0].size(); j++)
        {
          Eigen::MatrixXd temp_matrix(1,1);
          temp_matrix(0,0)=1;
					
          u3shell::RelativeUnitTensorLabelsU3ST unit_tensor=unit_sym_map[N0][j];
          u3::SU3 x0=unit_tensor.x0();
          HalfInt S0=unit_tensor.S0();
          int rp=unit_tensor.bra().eta();
          int r=unit_tensor.ket().eta();
          
          int rho0_max=u3::OuterMultiplicity(sigma.SU3(),x0, sigmap.SU3());
          for (int rho0=1; rho0<=rho0_max; rho0++)
            {
              if (rp<=(N1b+sp_irrepp.Nex()) && r<=(N1b+sp_irrep.Nex()))
                {
                  //std::cout<<unit_tensor.Str()<<std::endl;
                  temp_unit_map[spncci::UnitTensorU3Sector(sigmap,sigma,unit_tensor,rho0)]=temp_matrix;	
                }
            }
        }
      sp_irrep_unit_tensor_rme_map[sp_irrep_pair][N0_pair]=temp_unit_map;
      // std::cout<<"number of sp_irrep sectors "<<temp_unit_map.size()<<std::endl;;

      //////////////////////////////////////////////////////////////////////////////////////////////
      // Generating the rme's of the unit tensor for each SpIrrep
      std::map<std::pair<int,int>,std::vector<spncci::UnitTensorU3Sector>> unit_tensor_NpN_sector_map;
      GenerateUnitTensorU3SectorLabels(
        N1b,Nmax,sp_irrep_pair,sp_irrep_vector,
        unit_sym_map,unit_tensor_NpN_sector_map);

      spncci::GenerateUnitTensorMatrix(
        N1b,Nmax,sp_irrep_pair,sp_irrep_vector,u_coef_cache,k_matrix_map,
        unit_tensor_NpN_sector_map,sp_irrep_unit_tensor_rme_map[sp_irrep_pair]);
  }
  std::cout<<"all done"<<std::endl;
}
// end main 
