/****************************************************************
  recurrence_indexing.cpp

  Anna E. McCoy[1] and Patrick J. Fasano[2,3]
  [1] Institute for Nuclear Theory  
  [2] University of Notre Dame
  [3] Lawrence Berkeley National Laboratory

  SPDX-License-Identifier: MIT
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
    SpinSubspace::SpinSubspace(
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
        PushSubspace(SpinSubspace(S,spin_vector));
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
            const int& upsilon_max=u3subspace.size();
            // Nn_max_=Nn_max;
            PushSubspace(U3Subspace(u3subspace),upsilon_max);
          }
      }

    Space::Space(
          const std::vector<u3::U3>& sigma_vector, 
          const int Nsigma0,
          const int Nmax
      )
      {
       for(const auto& sigma : sigma_vector)
        {
          // Nsigma0_=Nsigma0;
          int Nn_max = Nmax-int(sigma.N()-Nsigma0);
          PushSubspace(LGISpace(sigma,Nn_max));
        } 
      }

    Space::Space(const spin::Space& spin_space, const int Nmax,const int Nsigma0)
      {
        // Nsigma0_=spin_space.Nsigma0(); 
        for(int i=0; i<spin_space.size(); ++i)
          {
            const u3::U3& sigma = spin_space.GetSubspace(i).sigma();
            int Nn_max = Nmax-int(sigma.N()-Nsigma0);
            PushSubspace(LGISpace(sigma,Nn_max));
          }

      }

    /////////////////////////////////////////////////////////////////////////////////////////////////////// 
    // Recurrence indexing (spatial)
    /////////////////////////////////////////////////////////////////////////////////////////////////////// 
    RecurrenceOperatorSubspace::RecurrenceOperatorSubspace(
      const u3::SU3& x0, 
      const std::vector<std::tuple<int,int>>& Nbar_pairs
    )
    {
      labels_=x0;
      for(const auto& Nbar_pair : Nbar_pairs)
        PushStateLabels(Nbar_pair);
    }

    RecurrenceU3Space::RecurrenceU3Space(
      const std::tuple<u3::U3,u3::U3> omega_pair, 
      const std::vector<std::tuple<u3::SU3,int,int>>& spatial_unit_tensors
    )
    {
      labels_=omega_pair;
      std::map<u3::SU3,std::vector<std::tuple<int,int>>> x0_Nbar_pairs;
      for(const auto& [x0,Nbar,Nbar_p] : spatial_unit_tensors)
          x0_Nbar_pairs[x0].push_back({Nbar,Nbar_p});
        
      
      for(const auto& [x0,Nbar_pairs] : x0_Nbar_pairs)
        {
         int rho0_max=u3::OuterMultiplicity(omega_ket().SU3(),x0,omega_bra().SU3());
         PushSubspace(RecurrenceOperatorSubspace(x0,Nbar_pairs),rho0_max); 
        }
    } 


    RecurrenceNnsumSpace::RecurrenceNnsumSpace(
      int Nnsum,
      const std::vector<std::tuple<int,int>> partition,
      const LGISpace& lgi_space_ket,
      const LGISpace& lgi_space_bra,
      const std::vector<std::tuple<u3::SU3,int,int>>& spatial_unit_tensors //pass through to Operator subspaces
    )
    {
      labels_=Nnsum;
      for(const auto&[i_ket,i_bra] : partition)
        {
          const u3::U3& omega_ket=lgi_space_ket.GetSubspace(i_ket).omega();
          const u3::U3& omega_bra=lgi_space_bra.GetSubspace(i_bra).omega();

          const int& upsilon_max_ket=lgi_space_ket.GetSubspaceDegeneracy(i_ket);
          const int& upsilon_max_bra=lgi_space_bra.GetSubspaceDegeneracy(i_bra);

          upsilon_pairs_.push_back({upsilon_max_ket,upsilon_max_bra});
          PushSubspace(
            RecurrenceU3Space({omega_ket,omega_bra},spatial_unit_tensors),
            upsilon_max_bra*upsilon_max_ket
            );
          }
    } 


    RecurrenceLGISpace::RecurrenceLGISpace(
      const LGISpace& lgi_space_ket,
      const LGISpace& lgi_space_bra,
      const std::vector<std::tuple<u3::SU3,int,int>>& spatial_unit_tensors //pass through to Operator subspaces
    )
    {
      const u3::U3& sigma_ket=lgi_space_ket.sigma();
      const u3::U3& sigma_bra=lgi_space_bra.sigma();
      labels_={sigma_ket,sigma_bra};

      // const HalfInt& Nsigma0_ket=lgi_space_ket.Nsigma0();
      // const HalfInt& Nsigma0_bra=lgi_space_bra.Nsigma0();

      std::map<int,std::vector<std::tuple<int,int>>> Nnsum_partition;
      for(int i_ket=0; i_ket<lgi_space_ket.size(); ++i_ket)
        for(int i_bra=0; i_bra<lgi_space_bra.size(); ++i_bra)
          {
            const u3::U3& omega_ket=lgi_space_ket.GetSubspace(i_ket).omega();

            const u3::U3& omega_bra=lgi_space_bra.GetSubspace(i_bra).omega();


            int Nnsum=int(omega_ket.N()-sigma_ket.N() + omega_bra.N()-sigma_bra.N());
            Nnsum_partition[Nnsum].push_back({i_ket,i_bra});
          }
      for(const auto& [Nnsum,partition] : Nnsum_partition)
        PushSubspace(
          RecurrenceNnsumSpace(
            Nnsum,partition,
            lgi_space_ket,lgi_space_bra,
            spatial_unit_tensors
            )
          );
    }

    RecurrenceSpace::RecurrenceSpace(
      const spatial::Space& space_ket,
      const spatial::Space& space_bra,
      const std::vector<std::tuple<u3::SU3,int,int>>& spatial_unit_tensors //pass through to Operator subspaces
    )
    {
      for(int i_ket=0; i_ket<space_ket.size(); ++i_ket)
        for(int i_bra=0; i_bra<space_bra.size(); ++i_bra)
          {
            const LGISpace& lgi_space_ket=space_ket.GetSubspace(i_ket);
            const LGISpace& lgi_space_bra=space_bra.GetSubspace(i_bra);
            PushSubspace(RecurrenceLGISpace(lgi_space_ket,lgi_space_bra,spatial_unit_tensors));
          }
    }




  } //end spatial namespace



}

