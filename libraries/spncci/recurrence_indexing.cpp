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
  
    Space::Space(const lgi::MultiplicityTaggedLGIVector& lgi_vector,int Nmax)
    {
      std::map<u3::U3,std::map<HalfInt,MultiplicityTagged<SpSn>::vector>> sigma_spins_map;
      for(const auto& [lgi,gamma_max] : lgi_vector)
        {
          const auto& [Nex,sigma,Sp,Sn,S] = lgi.Key();
          if(Nex>Nmax)
            continue;

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
            labels_=omega;
            const auto&n_rho = u3subspace.GetStateLabels(i);
            PushStateLabels(n_rho);
          }

      }

    LGISpace::LGISpace(const u3::U3& sigma, const int Nn_max)
      {
        sp3r::Sp3RSpace sp3r_space(sigma, Nn_max);
        for(int i=0; i<sp3r_space.size(); i++ )
          {
            labels_=sigma;
            const auto& u3subspace=sp3r_space.GetSubspace(i);
            PushSubspace(U3Subspace(u3subspace));
          }
      }

    Space::Space(
          const std::vector<u3::U3>& sigma_vector, 
          const HalfInt& Nsigma0,
          const int Nmax
      )
      {
       for(const auto& sigma : sigma_vector)
        {
          // Nsigma0_=Nsigma0;
          int Nn_max = Nmax-int(sigma.N()-Nsigma0);
          if (Nn_max>=0)
            PushSubspace(LGISpace(sigma,Nn_max));
        } 
      }

    Space::Space(const spin::Space& spin_space, const int& Nmax,const HalfInt& Nsigma0)
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
          const std::tuple<u3::U3,u3::U3>& omega_pair,
          const UnitTensorParameters& unit_tensor_parameters
        )
    {
      labels_=omega_pair;
      const auto&[omega,omega_p]=omega_pair;
      std::map<u3::SU3,std::vector<std::tuple<int,int>>> x0_Nbar_pairs;
      
      ////////////////////////////////////////////////////////////////////////////////
      // Create list of spatial unit tensors to pass through to operator subspace
      int Nbar_max(2*unit_tensor_parameters.N1v+omega.N()
                      -unit_tensor_parameters.Nsigma0);
      
      int Nbar_p_max(2*unit_tensor_parameters.N1v+omega_p.N()
                      -unit_tensor_parameters.Nsigma0);

      int N0(omega_p.N()-omega.N());
      int Nbar_min=unit_tensor_parameters.state_parity;
      for(int Nbar=Nbar_min; Nbar<=Nbar_max; Nbar+=2)
          {
            int Nbar_p=N0+Nbar;
            if (Nbar_p<=Nbar_p_max & Nbar_p>=0)
              {
                MultiplicityTagged<u3::SU3>::vector possible_x0
                      =u3::KroneckerProduct({Nbar_p,0},{0,Nbar});
              
                for(const auto& [x0,rho] : possible_x0)
                  x0_Nbar_pairs[x0].push_back({Nbar,Nbar_p});
              
              }
          }

      ////////////////////////////////////////////////////////////////////////////////
      // Construct operator subspaces 
      ////////////////////////////////////////////////////////////////////////////////
      for(const auto& [x0,Nbar_pairs] : x0_Nbar_pairs)
        {
         int rho0_max=u3::OuterMultiplicity(omega_ket().SU3(),x0,omega_bra().SU3());
         if(rho0_max>0)
          PushSubspace(RecurrenceOperatorSubspace(x0,Nbar_pairs),rho0_max); 
        }
    } 


    RecurrenceNnsumSpace::RecurrenceNnsumSpace(
      int Nnsum,
      const std::vector<std::tuple<int,int>> u3subspace_index_pairs,
      const LGISpace& lgi_space_ket,
      const LGISpace& lgi_space_bra,
      const UnitTensorParameters& unit_tensor_parameters
    )
    {
      labels_=Nnsum;
      unit_tensor_state_parity_=unit_tensor_parameters.state_parity;

      for(const auto&[i_ket,i_bra] : u3subspace_index_pairs)
        {
          const auto& u3subspace_ket=lgi_space_ket.GetSubspace(i_ket);
          const auto& u3subspace_bra=lgi_space_bra.GetSubspace(i_bra);
          const u3::U3& omega_ket=u3subspace_ket.omega();
          const u3::U3& omega_bra=u3subspace_bra.omega();

          const int& upsilon_max_ket=u3subspace_ket.size();
          const int& upsilon_max_bra=u3subspace_bra.size();

          upsilon_pairs_.push_back({upsilon_max_ket,upsilon_max_bra});

          PushSubspace(
            RecurrenceU3Space({omega_ket,omega_bra},unit_tensor_parameters),
            upsilon_max_bra*upsilon_max_ket
            );
          }
    } 


    RecurrenceLGISpace::RecurrenceLGISpace(
      const LGISpace& lgi_space_ket,
      const LGISpace& lgi_space_bra,
      const int& N1v, const HalfInt& Nsigma0
    )
    {
      const u3::U3& sigma_ket=lgi_space_ket.sigma();
      const u3::U3& sigma_bra=lgi_space_bra.sigma();
      // std::cout<<sigma_bra.Str()<<"  "<<sigma_ket.Str()<<std::endl;
      labels_={sigma_ket,sigma_bra};

      // Partition pairs of omega',omega by Nnsum
      std::map<int,std::vector<std::tuple<int,int>>> Nnsum_partition;
      for(int i_ket=0; i_ket<lgi_space_ket.size(); ++i_ket)
        for(int i_bra=0; i_bra<lgi_space_bra.size(); ++i_bra)
          {
            const u3::U3& omega_ket=lgi_space_ket.GetSubspace(i_ket).omega();
            const u3::U3& omega_bra=lgi_space_bra.GetSubspace(i_bra).omega();

            int Nnsum=int(omega_ket.N()-sigma_ket.N() + omega_bra.N()-sigma_bra.N());
            Nnsum_partition[Nnsum].push_back({i_ket,i_bra});
          }

      // Create RecurrenceNnsumSpaces.  On for each unit tensor state parity
      for(int unit_tensor_state_parity=0; unit_tensor_state_parity<=1; ++unit_tensor_state_parity)
        {
          UnitTensorParameters unit_tensor_parameters(N1v,Nsigma0,unit_tensor_state_parity);
          for(const auto& [Nnsum,partition] : Nnsum_partition)
            PushSubspace(
              RecurrenceNnsumSpace(
                Nnsum,partition,
                lgi_space_ket,lgi_space_bra,
                unit_tensor_parameters
                )
              );
        }
    }

    RecurrenceSpace::RecurrenceSpace(
      const spatial::Space& space_ket,
      const spatial::Space& space_bra,
      const int& N1v, const HalfInt& Nsigma0
    )
    {
      for(int i_ket=0; i_ket<space_ket.size(); ++i_ket)
        for(int i_bra=0; i_bra<space_bra.size(); ++i_bra)
          {
            const LGISpace& lgi_space_ket=space_ket.GetSubspace(i_ket);
            const LGISpace& lgi_space_bra=space_bra.GetSubspace(i_bra);
            PushSubspace(RecurrenceLGISpace(lgi_space_ket,lgi_space_bra,N1v,Nsigma0));
          }
    }




  } //end spatial namespace



}

