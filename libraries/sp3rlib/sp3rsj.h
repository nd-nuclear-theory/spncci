/****************************************************************
  sp3rsj.h

  Sp(3,R)xSU(2)>SU(2) labeling and branching restricted to definite S and J.

  Colin V. Coane

  SPDX-License-Identifier: MIT

  10/17/23 (cvc): Created.
  01/25/24 (cvc): Added class for operator space and subspace structure.
  06/03/24 (cvc): Added class for Sp(3,R)xSU(2)>SU(2) observable.
  06/10/24 (cvc): Fixed multiplicity indexing for hypersectors.
  07/05/24 (cvc): Added bookkeeping for observables as linear combinations of operators.
  07/08/24 (cvc): Added ability to compute Sp(3,R)xSU(2)>SU(2) RMEs to OperatorInfoSp3RSJ.

****************************************************************/

#ifndef SP3RSJ_H_
#define SP3RSJ_H_

#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "basis/basis.h"
#include "basis/degenerate.h"
#include "basis/hypersector.h"
#include "basis/operator.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3boson.h"
#include "sp3rlib/sp3r.h"
#include "am/halfint.h"
#include "am/am.h"

#include <functional>
#include <fstream>
#include <iterator>

namespace sp3r
{
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Class for carrier space of Sp(3,R)>U(3)>SO(3)xSU(2)>SU(2) irrep restricted to S and J
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  class Sp3RSJSpace;

  class Sp3RSJSpace
      : public basis::BaseSpace<Sp3RSJSpace, sp3r::U3Subspace, std::tuple<u3::U3>>
  {
  public:
    Sp3RSJSpace() = default;
    Sp3RSJSpace(
        const u3::U3& sigma,
        unsigned int Nn_max,
        const HalfInt& s_const,
        const HalfInt& j_const,
        const bool cache_Kmatrices = true,
        // const bool subspace_labels_only = false,
        const bool branch_to_so3 = true
      );

    // accessors
    const u3::U3& sigma() const { return std::get<0>(labels()); }

    unsigned int Nn_max() const { return Nn_max_; }

    HalfInt S() const { return s_const_; }
    HalfInt J() const { return j_const_; }

    // Construct and return unrestricted Sp3RSpace
    sp3r::Sp3RSpace Sp3RSpace() const { return sp3r::Sp3RSpace(sigma(),Nn_max()); }
    //sp3r::Sp3RSpace Sp3RSpace() const { return sp3r_space_; }
    
    // Get dimension in fully branched basis
    std::size_t GetBranchedDimension() const
    {
      // Initialize full space dimension
      std::size_t full_dimension = 0;
      for (int i_subspace = 0; i_subspace < size(); i_subspace++)
      {
        // Loop through subspaces
        const sp3r::U3Subspace& subspace = GetSubspace(i_subspace);
        full_dimension += subspace.total_states()*subspace.dimension();
      }
      // Return dimension of space
      return full_dimension;
    }

    // Get subspace offset in fully branched basis
    std::size_t GetBranchedSubspaceOffset(std::size_t subspace_idx) const 
    {
      // Initialize branched offset
      std::size_t subspace_offset = 0;
      for (std::size_t i_subspace = 0; i_subspace < subspace_idx; i_subspace++)
      {
        // Continue if reached end of space
        if (i_subspace >= size())
          continue;
        // Loop through subspaces
        const sp3r::U3Subspace& subspace = GetSubspace(i_subspace);
        subspace_offset += subspace.total_states()*subspace.dimension();
      }
      return subspace_offset;
    }

    // diagnostic output
    std::string DebugStr() const;

  private:
    HalfInt s_const_;
    HalfInt j_const_;
    unsigned int Nn_max_;
    //sp3r::Sp3RSpace sp3r_space_; // Unrestricted Sp3RSpace for given sigma, Nn_max
  };


  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Sectors: U3Subspaces connected by operator w0
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  class Sp3RSJSectors
      : public basis::BaseSectors<Sp3RSJSpace>
  {
  public:
    // Default constructor
    Sp3RSJSectors() = default;

    // Constructor
    Sp3RSJSectors(
        const Sp3RSJSpace& space, const u3::U3& omega0, const bool& su3_generator = false
      );

  private:
    u3::U3 omega0_;
  };

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Operator Space and Subspaces
  // sp3r::OperatorSpaceSp3RSJ [J0,S0]
  //    -> sp3r::OperatorU3Subspace [N0,x0]
  //      -> sp3r::OperatorStateSp3RSJ [L, kappa]
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  class OperatorU3Subspace;
  class OperatorSpaceSp3RSJ;
  class OperatorStateSp3RSJ;


  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Operator Subpace: U(1)xSU(3) labels restricted to good J0,S0 
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  class OperatorU3Subspace
      : public basis::
            BaseDegenerateSubspace<OperatorU3Subspace, std::tuple<u3::U3>, OperatorStateSp3RSJ, std::tuple<unsigned int, unsigned int>>
  {
  public:
    OperatorU3Subspace() = default;
    // constructor

    inline void GenerateOpStates(
      const u3::U3& omega,
      const std::pair<unsigned int,unsigned int>& L_min_max
    )
    {
      const auto& L_kappa_vector = u3::BranchingSO3(omega.SU3());
      const auto& [Lmin,Lmax] = L_min_max;
      for (const auto& [L, kappa_max] : L_kappa_vector)
      {
        if(L >= Lmin && L<= Lmax)
        {
          // Push states with all kappa <= kappa_max
          for (int kappa = 1; kappa <= kappa_max; kappa++)
          {
            // make L,kappa tuple
            std::tuple<unsigned int, unsigned int> L_k = {L,kappa};
            // push SU(3)>SO(3) labels and enforce multiplicity = 1
            // creates non-degenerate states with multiplicity as part of label
            // change this to be a non-degenerate subspace at some point??
            PushStateLabels(L_k, 1);
          }
        }
      }
    }

    inline OperatorU3Subspace(
      const HalfInt& s_op,
      const HalfInt& j_op,
      const u3::U3& u3_op,
      const bool branch_to_states,
      const std::pair<unsigned int,unsigned int>& L_min_max={0,so3::kNone}
    )
      : BaseDegenerateSubspace{u3_op}, s_op_(s_op), j_op_(j_op), l_pair_(L_min_max)
      // : BaseDegenerateSubspace{SubspaceLabelsType(s_op,j_op,u3_op)} // is this right???
    {
      if (branch_to_states)
      {

        GenerateOpStates(u3_op,L_min_max);
      }
    }

    // accessors
    const u3::U3& omega() const { return std::get<0>(labels()); }
    const u3::U3& U3() const { return omega(); }
    HalfInt S0() const { return s_op_; }
    HalfInt J0() const { return j_op_; }
    MultiplicityTagged<unsigned int>::vector GetOpStateLabels() const
    {
      // Branch without L x S -> J coupling
      const u3::U3 omega_op = omega();
      const auto& L_kappa_vector = u3::BranchingSO3(omega_op.SU3());
      // Get bounds for L value range
      const auto& [Lmin,Lmax] = l_pair_;

      // allocate container for product
      MultiplicityTagged<unsigned int>::vector branching;
      // int max_entries = static_cast<std::size_t>(L_max - L_min + 1);
      branching.reserve(static_cast<std::size_t>(Lmax - Lmin + 1));

      for (const auto& [L, kappa_max] : L_kappa_vector)
      {
        if(L >= Lmin && L<= Lmax)
        {
          // Push states
          branching.push_back(MultiplicityTagged<unsigned int>(L, kappa_max));
        }
      }

      // Return Operator State Labels
      return branching;
    };

    std::string DebugStr() const;

  private:
    std::pair<unsigned int,unsigned int> l_pair_;
    HalfInt s_op_;
    HalfInt j_op_;

  };

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Operator Space: Good J0, S0, and vector of operator U(1)xSU(3) labels
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  class OperatorSpaceSp3RSJ
      : public basis::BaseSpace<OperatorSpaceSp3RSJ, OperatorU3Subspace, std::tuple<HalfInt,HalfInt>>
  {
  public:
    OperatorSpaceSp3RSJ() = default;

    inline OperatorSpaceSp3RSJ(
      const HalfInt& s_op,
      const HalfInt& j_op,
      const std::vector<u3::U3>& u3_op_labels,
      bool branch_to_states = true
    ) 
      : BaseSpace{std::make_tuple(s_op,j_op)}, s_op_(s_op), j_op_(j_op)
    {
      //
      const HalfInt::pair l_range = am::ProductAngularMomentumRange(s_op,j_op);
      // Check allowed coupling of S and J results in integer values of L
      assert(IsInteger(l_range.second));
      // Allowed angular momentum range
      const std::pair<unsigned int,unsigned int>& l_pair={l_range.first.TwiceValue()/2,l_range.second.TwiceValue()/2};
      const auto& [Lmin,Lmax] = l_pair;
      // Loop through U(3) labels and push subspaces
      for (auto& u3_label : u3_op_labels)
      {
        // Check branching
        bool non_empty_subspace = false;
        const auto& subspace_so3_labels = u3::BranchingSO3(u3_label.SU3());
        // Check if subspace is empty
        for (const auto& [L, kappa_max] : subspace_so3_labels)
        {
          if (L >= Lmin && L<= Lmax)
            non_empty_subspace = true;
        };
        // Push subspace if non-empty
        if (non_empty_subspace)
        {
          PushSubspace(OperatorU3Subspace(s_op,j_op,u3_label,branch_to_states,l_pair));
        }
      }
    }

    // accessors
    HalfInt S0() const { return s_op_; }
    HalfInt J0() const { return j_op_; }
    //const std::vector<u3::U3>& U3Labels() const { return std::get<2>(labels()); }
    //HalfInt S0() const { return std::get<0>(labels()); }
    //HalfInt J0() const { return std::get<1>(labels()); }
    //const std::vector<u3::U3>& U3Labels() const { return std::get<2>(labels()); }

    // diagnostic output
    std::string DebugStr() const;

  private:
    HalfInt s_op_;
    HalfInt j_op_;

  };


  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Operator State: SO(3) and multiplicity labels for state
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  class OperatorStateSp3RSJ
      : public basis::BaseDegenerateState<OperatorU3Subspace>
  {
  public:
    OperatorStateSp3RSJ() = default;
    // pass-through constructors

    OperatorStateSp3RSJ(const SubspaceType& subspace, std::size_t index)
    // Construct state by index.
        : basis::BaseDegenerateState<OperatorU3Subspace>(subspace, index)
    {}

    OperatorStateSp3RSJ(
        const SubspaceType& subspace,
        const typename SubspaceType::StateLabelsType& state_labels
      )
    // Construct state by reverse lookup on labels.
        : basis::BaseDegenerateState<OperatorU3Subspace>(subspace, state_labels)
    {}

    // pass-through accessors for subspace labels
    unsigned int L() const { return std::get<0>(labels()); }
    unsigned int kappa() const { return std::get<1>(labels()); }

  private:

  };

  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Physical Operator and Observable Content
  //
  // sp3r::OperatorInfoSp3RSJ [OperatorSpaceSp3RSJ, OperatorRME]
  //    -> Stores operator space and RME function for a single operator
  //
  // sp3r::ObservableSp3RSJ 
  //    -> Stores Operators and coefficients as a vector
  //    -> Observable = sum_i {coef_i * op_i}
  //
  /////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Stores operator space and function for evaluating RMEs
  class OperatorInfoSp3RSJ
  {
  public:
    OperatorInfoSp3RSJ() = default;

    // Constructor for Sp(3,R) operator / U(3) tensor
    OperatorInfoSp3RSJ(
      const OperatorSpaceSp3RSJ& op_space,
      std::function<basis::OperatorBlock<double>(
      const u3::U3&,
      int,
      const sp3r::Sp3RSpace&,
      const sp3r::U3Subspace&,
      const sp3r::U3Subspace&,
      u3::UCoefCache&
      )> obs_func
    ) : op_space_(op_space), obs_func_(obs_func), type_obs_(1)
    {}

    // Constructor for Sp(3,R)xSU(2)>SU(2) operator (S,J dependent)
    OperatorInfoSp3RSJ(
      const OperatorSpaceSp3RSJ& op_space,
      std::function<basis::OperatorBlock<double>(
      const u3::U3&,
      int,
      const sp3r::Sp3RSpace&,
      const sp3r::U3Subspace&,
      const sp3r::U3Subspace&,
      u3::UCoefCache&,
      const HalfInt&,
      const HalfInt&,
      const HalfInt&,
      const HalfInt&
      )> obs_func_sj
    ) : op_space_(op_space), obs_func_sj_(obs_func_sj), type_obs_(2)
    {}

    // accessors
    // Return operator space labels
    const sp3r::OperatorSpaceSp3RSJ OperatorSpaceSp3RSJ() const { return op_space_; }

    // Compute RMEs for given operator function
    basis::OperatorBlock<double> ComputeRMESector(
      const u3::U3& omega,
      int rho,
      const sp3r::Sp3RSpace& base_space,
      const sp3r::U3Subspace& bra_subspace,
      const sp3r::U3Subspace& ket_subspace,
      u3::UCoefCache& u_coef_cache,
      const HalfInt& S_bra = HalfInt(0,2),
      const HalfInt& S_ket = HalfInt(0,2),
      const HalfInt& J_bra = HalfInt(0,2),
      const HalfInt& J_ket = HalfInt(0,2)
    ) const
    {
      // Check if operator is dependent on spin S and total angular momentum J
      if (type_obs_ == 1)
      {
        // No S and J dependence
        return obs_func_(omega,rho,base_space,
          bra_subspace,ket_subspace,u_coef_cache
        );
      }
      else if (type_obs_ == 2)
      {
        // S and J dependence
        return obs_func_sj_(omega,rho,base_space,
          bra_subspace,ket_subspace,u_coef_cache,
          S_bra,S_ket,J_bra,J_ket
        );
      }
      else
      {
        return basis::OperatorBlock<double>::Zero(
          bra_subspace.dimension(),
          ket_subspace.dimension()
        );
      }
    }

  private:
    const sp3r::OperatorSpaceSp3RSJ& op_space_;
    std::function<basis::OperatorBlock<double>(
      const u3::U3&,
      int,
      const sp3r::Sp3RSpace&,
      const sp3r::U3Subspace&,
      const sp3r::U3Subspace&,
      u3::UCoefCache&
      )> obs_func_;
      std::function<basis::OperatorBlock<double>(
      const u3::U3&,
      int,
      const sp3r::Sp3RSpace&,
      const sp3r::U3Subspace&,
      const sp3r::U3Subspace&,
      u3::UCoefCache&,
      const HalfInt&,
      const HalfInt&,
      const HalfInt&,
      const HalfInt&
      )> obs_func_sj_;
    int type_obs_;

  };
  

  // Observable is linear combination of operators and coefficients
  // Obs = Observable = sum_i {coef_i * op_i}
  // Store observable as vector of pairs: {operator,coefficient}
  typedef std::vector<std::pair<sp3r::OperatorInfoSp3RSJ, double>> ObservableSp3RSJ;


  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Hypersectors: U3Subspaces connected by operators with good J0, S0, SU(3)>SO(3)
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  class Sp3RSJHypersectors
      : public basis::BaseHypersectors<Sp3RSJSpace,OperatorSpaceSp3RSJ>
  {
  public:
    // Default constructor
    Sp3RSJHypersectors() = default;
    
    // Constructor
    Sp3RSJHypersectors(
        const Sp3RSJSpace& space, const OperatorSpaceSp3RSJ& operator_space, const bool& su3_generator = false
      );
  };



}

#endif
