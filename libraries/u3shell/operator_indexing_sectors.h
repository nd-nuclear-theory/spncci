/****************************************************************
  operator_indexing_sectors.h

  Enumeration of operator sectors
                                  
  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  3/23/22 (aem): Created.
****************************************************************/

#ifndef OPERATOR_INDEXING_SECTORS_H_
#define OPERATOR_INDEXING_SECTORS_H_

#include <unordered_set>
#include "basis/basis.h"
#include "basis/degenerate.h"

namespace u3shell::relative
{
////////////////////////////////////////////////////////////////
/// Selection rules:
///  parity_bar==(exchange_symm_bar + 1)%2.
///  S0 x L0 -> J0
///
/// OperatorSpace
///  -> OperatorSubspace [S0,exhcange_sym_bar]
///    -> OperatorStates [i]
///
/// spatial Space
///  OperatorSpace
///    ->OperatorParitySpace [parity_bar]
///      -> OperatorN0Space [N0]
///        -> OperatorL0Space [L0]
///           -> OperatorSU3Subspace (kappa0)[omega0]
///               ->OperatorStates [Nbar]  Nbarp fixed by N0.
///
////////////////////////////////////////////////////////////////

template<typename tOperatorSpatialSpaceType, typename tOperatorSpatialSubspaceType,
 typename tOperatorSpinSpaceType,typename tOperatorSpinSubspaceType>
class OperatorU3SpinSectors
  : public basis::BaseSectors<
      tOperatorSpatialSpaceType,
      tOperatorSpinSpaceType,
      basis::BaseSector<tOperatorSpatialSubspaceType,tOperatorSpinSubspaceType>
      >
{
private:
      // Constructor
    using BaseSectorsType = basis::BaseSectors<
      tOperatorSpatialSpaceType,
      tOperatorSpinSpaceType,
      basis::BaseSector<tOperatorSpatialSubspaceType,tOperatorSpinSubspaceType>
      >;

  public:

    // Bring bra_space and ket_space into current name space
    using BaseSectorsType::bra_space;
    using BaseSectorsType::ket_space;

    // Default constructor
    OperatorU3SpinSectors() = default;

    OperatorU3SpinSectors(
      std::shared_ptr<const tOperatorSpatialSpaceType> spatial_space_ptr,
      std::shared_ptr<const tOperatorSpinSpaceType> spin_space_ptr,
      const unsigned int J0
    )
    : BaseSectorsType{std::move(spatial_space_ptr),std::move(spin_space_ptr)},J0_(J0)
    {
      for (const auto& spatial_subspace : bra_space())
        for (int i_spin = 0; i_spin < ket_space().size(); ++i_spin) {
          const auto& spin_subspace = ket_space().GetSubspace(i_spin);
          if (spatial_subspace.parity_bar() != (spin_subspace.exchange_symm_bar() + 1) % 2)
            continue;

          const auto& S0 = spin_subspace.S0();
          for (const auto& N0subspace : spatial_subspace)
            for (int i_spatial = 0; i_spatial < N0subspace.size(); ++i_spatial)
            {
              const auto& L0subspace = N0subspace.GetSubspace(i_spatial);
              if (am::AllowedTriangle(S0, L0subspace.L0(), J0))
              {
                auto spatial_ptr = N0subspace.GetSubspacePtr(i_spatial);
                auto spin_ptr = ket_space().GetSubspacePtr(i_spin);

                BaseSectorsType::EmplaceSector(i_spatial, i_spin,spatial_ptr,spin_ptr,1);
              }
            }
       }
    }

    inline unsigned int J0() const {return J0_;}


private:
  unsigned int J0_;
};

}//relative namespace



#endif
