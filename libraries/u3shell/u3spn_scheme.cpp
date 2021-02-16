/****************************************************************
  u3spn_scheme.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT
****************************************************************/

#include "u3shell/u3spn_scheme.h"

#include <sstream>

#include "am/am.h"
#include "fmt/format.h"

namespace u3shell {

  std::string U3SPN::Str() const
  {
    std::ostringstream ss;

    ss << omegaS_.Str()
       << "x" << "(" << Sp_ << "x" << Sn_ << ")";

    return ss.str();
  }

  SubspaceU3SPN::SubspaceU3SPN (const u3shell::U3SPN& omegaSPN, int dimension)
  {
    // set labels
    labels_ = omegaSPN;

    // save dimension
    dimension_ = dimension;
  }

  std::string SubspaceU3SPN::LabelStr() const
  {

    return fmt::format(
                       "[{}]",
                       U3SPN().Str()
                       );
  }

  SpaceU3SPN::SpaceU3SPN(const std::map<u3shell::U3SPN,int>& subspace_dimensions)
  {
    for (auto& omegaSPN_dimension : subspace_dimensions)
      {

        // define aliases for key and value
        const u3shell::U3SPN& omegaSPN = omegaSPN_dimension.first;
        const int& dimension = omegaSPN_dimension.second;

        // construct subspace
        SubspaceU3SPN subspace(omegaSPN,dimension);

        // push subspace if nonempty
        assert(subspace.size()!=0);
        PushSubspace(subspace);
      }
  }

  SectorsU3SPN::SectorsU3SPN(const SpaceU3SPN& space, const OperatorLabelsU3S& operator_labels,
                             bool spin_scalar)
    : BaseSectors(space)
  {
    for (int bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index)
      for (int ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index)
        {
          // retrieve subspaces
          const SubspaceU3SPN& bra_subspace = space.GetSubspace(bra_subspace_index);
          const SubspaceU3SPN& ket_subspace = space.GetSubspace(ket_subspace_index);

          // verify selection rules
          bool allowed = true;
          // U(1)
          allowed &= (ket_subspace.N() + operator_labels.N0() - bra_subspace.N() == 0);
          // spin
          allowed &= am::AllowedTriangle(ket_subspace.S(),operator_labels.S0(),bra_subspace.S());
          // proton and neutron spin
          if (spin_scalar)
            {
              assert(operator_labels.S0()==0);
              allowed &= (ket_subspace.Sp()==bra_subspace.Sp()) && (ket_subspace.Sn()==bra_subspace.Sn());
              }
          // find SU(3) multiplicity and check SU(3) selection
          int multiplicity = 0;
          if (allowed)
            {
              multiplicity = u3::OuterMultiplicity(ket_subspace.SU3(),operator_labels.x0(),bra_subspace.SU3());
              allowed &= (multiplicity > 0);
            }

          // push sectors (tagged by multiplicity)
          if (allowed)
            for (int multiplicity_index = 1; multiplicity_index <= multiplicity; ++multiplicity_index)
              {
                PushSector(bra_subspace_index,ket_subspace_index,multiplicity_index);
              }
        }
  }

  SectorsU3SPN::SectorsU3SPN(
    const SpaceU3SPN& space_bra, const SpaceU3SPN& space_ket, 
    const OperatorLabelsU3S& operator_labels,bool spin_scalar
  )
  {
    for (int bra_subspace_index=0; bra_subspace_index<space_bra.size(); ++bra_subspace_index)
      for (int ket_subspace_index=0; ket_subspace_index<space_ket.size(); ++ket_subspace_index)
        {
          // retrieve subspaces
          const SubspaceU3SPN& bra_subspace = space_bra.GetSubspace(bra_subspace_index);
          const SubspaceU3SPN& ket_subspace = space_ket.GetSubspace(ket_subspace_index);

          // verify selection rules
          bool allowed = true;
          // U(1)
          allowed &= (ket_subspace.N() + operator_labels.N0() - bra_subspace.N() == 0);
          // spin
          allowed &= am::AllowedTriangle(ket_subspace.S(),operator_labels.S0(),bra_subspace.S());
          // proton and neutron spin
          if (spin_scalar)
            {
              assert(operator_labels.S0()==0);
              allowed &= (ket_subspace.Sp()==bra_subspace.Sp()) && (ket_subspace.Sn()==bra_subspace.Sn());
              }
          // find SU(3) multiplicity and check SU(3) selection
          int multiplicity = 0;
          if (allowed)
            {
              multiplicity = u3::OuterMultiplicity(ket_subspace.SU3(),operator_labels.x0(),bra_subspace.SU3());
              allowed &= (multiplicity > 0);
            }

          // push sectors (tagged by multiplicity)
          if (allowed)
            for (int multiplicity_index = 1; multiplicity_index <= multiplicity; ++multiplicity_index)
              {
                PushSector(SectorType(bra_subspace_index,ket_subspace_index,bra_subspace,ket_subspace,multiplicity_index));
              }
        }
  }



  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
