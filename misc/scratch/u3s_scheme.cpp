/****************************************************************
  u3s_scheme.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "u3shell/u3s_scheme.h"

#include <sstream>

#include "am/am.h"
#include "fmt/format.h"

namespace u3shell {

  SubspaceU3S::SubspaceU3S (const u3::U3S& omegaS, int dimension)
  {
    // set labels
    labels_ = omegaS;

    // save dimension
    dimension_ = dimension;
  }

  std::string SubspaceU3S::Str() const
  {

    return fmt::format(
                       "[{}]",
                       U3S().Str()
                       );
  }

  SpaceU3S::SpaceU3S(const std::map<u3::U3S,int>& subspace_dimensions)
  {
    for (auto& omegaS_dimension : subspace_dimensions)
      {

        // define aliases for key and value
        const u3::U3S& omegaS = omegaS_dimension.first;
        const int& dimension = omegaS_dimension.second;

        // construct subspace
        SubspaceU3S subspace(omegaS,dimension);

        // push subspace if nonempty
        assert(subspace.size()!=0);
        PushSubspace(subspace);
      }
  }
  
  SectorsU3S::SectorsU3S(SpaceU3S& space, const OperatorLabelsU3S& operator_labels)
  {
    for (int bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index)
      for (int ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index)
        {
          // retrieve subspaces
          const SubspaceU3S& bra_subspace = space.GetSubspace(bra_subspace_index);
          const SubspaceU3S& ket_subspace = space.GetSubspace(ket_subspace_index);
          // std::cout << bra_subspace.Str() << ket_subspace.Str() << std::endl;

          // verify selection rules
          bool allowed = true;
          // U(1)
          allowed &= (ket_subspace.N() + operator_labels.N0() - bra_subspace.N() == 0);
          // spin
          allowed &= am::AllowedTriangle(ket_subspace.S(),operator_labels.S0(),bra_subspace.S());
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
