/****************************************************************
  recurrence_indexing.cpp

  Anna E. McCoy and Patrick J. Fasano
  Institute for Nuclear Theory  
  and University of Notre Dame

****************************************************************/

#include "spncci/recurrence_indexing.h"
#include "basis/basis.h"
#include <fstream>
#include <iostream>
#include <algorithm>
// #include "mcutils/parsing.h"
#include "fmt/format.h"
// #include "am/halfint_fmt.h"

// #include "sp3rlib/vcs.h"

namespace spncci
{
  namespace spin
  {
    SpNCCISpinSubspace::SpNCCISpinSubspace(
      const HalfInt& S, 
      const MultiplicityTagged<SpSn>::vector& spin_vector
    )
    {
      labels_ = {S};
      for(const auto& [SpSn,gamma_max] : spin_vector)
        PushStateLabels(SpSn,gamma_max); 
    }

    LGISpace::LGISpace(
      const u3::U3& sigma, 
      const std::map<HalfInt,MultiplicityTagged<SpSn>::vector>& spin_map
      )
    {
      labels_ = {sigma};
      for(const auto& [S,spin_vector] : spin_map)
        PushSubspace(SpNCCISpinSubspace(S,spin_vector));
    }
  
    Space::Space(const lgi::MultiplicityTaggedLGIVector& lgi_vector)
    {
      std::map<u3::U3,std::map<HalfInt,MultiplicityTagged<SpSn>::vector>> sigma_spins_map;
      for(const auto& [lgi,gamma_max] : lgi_vector)
        {
          const auto& [Nex,sigma,Sp,Sn,S] = lgi.Key();
          
          if (Nex==0)
            Nsigma0_=sigma.N();

          sigma_spins_map[sigma][S].emplace_back(SpSn{Sp,Sn},gamma_max);
        }
  
      for(const auto& [sigma,spin_map] : sigma_spins_map)
          PushSubspace(LGISpace(sigma,spin_map));
    }
  }
  //////////////////////////////////////////////////////////////////////////////////////////
  namespace spatial
    {      
      U3Subspace::U3Subspace(const sp3r::U3Subspace& u3subspace)
        {
          const u3::U3& omega=u3subspace.labels();
          for(int i=0; i<u3subspace.size(); ++i)
            {
              const auto&n_rho = u3subspace.GetStateLabels(i);
              PushStateLabels(n_rho);
            }

        }

      LGISpace::LGISpace(const u3::U3& sigma, const int Nn_max)
        {
          sp3r::Sp3RSpace sp3r_space(sigma, Nn_max);
          for(int i=0; i<sp3r_space.size(); i++ )
            {
              const auto& u3subspace=sp3r_space.GetSubspace(i);
              PushSubspace(U3Subspace(u3subspace));
            }
        }

      Space::Space(
            const std::vector<u3::U3>& sigma_vector, 
            const int N0,
            const int Nmax
        )
        {
         for(const auto& sigma : sigma_vector)
          {
            int Nn_max = Nmax-int(sigma.N()-N0);
            PushSubspace(LGISpace(sigma,Nn_max));
          } 
        }

      Space::Space(const spin::Space& spin_space, const int Nmax)
        {
          HalfInt Nsigma0=spin_space.Nsigma0(); 
          for(int i=0; i<spin_space.size(); ++i)
            {
              const u3::U3& sigma = spin_space.GetSubspace(i).sigma();
              int Nn_max = Nmax-int(sigma.N()-Nsigma0);
              PushSubspace(LGISpace(sigma,Nn_max));
            }
        }
    }


}

