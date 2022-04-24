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

template<typename tOperatorSpatialSubspaceType,typename tOperatorSpinSubspaceType>
class OperatorU3SpinSector
  : public basis::BaseSector<tOperatorSpatialSubspaceType,tOperatorSpinSubspaceType>
{
private:
  using BaseSectorType =
     basis::BaseSector<tOperatorSpatialSubspaceType, tOperatorSpinSubspaceType>;

public:
 using BaseSectorType::bra_subspace;
 using BaseSectorType::bra_subspace_index;
 using BaseSectorType::ket_subspace;
 using BaseSectorType::ket_subspace_index;

 using KeyType = std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>;

 OperatorU3SpinSector() = default;

 OperatorU3SpinSector(
     const std::size_t spatial_subspace_index,
     const std::size_t spin_subspace_index,
     std::shared_ptr<const tOperatorSpatialSubspaceType> spatial_subspace_ptr,
     std::shared_ptr<const tOperatorSpinSubspaceType> spin_subspace_ptr,
     const std::size_t parity_space_index,
     const std::size_t N0_space_index
   )
     : BaseSectorType{spatial_subspace_index, spin_subspace_index, std::move(spatial_subspace_ptr), std::move(spin_subspace_ptr)},
       parity_space_index_{parity_space_index},
       N0_space_index_{N0_space_index}
 {}

 inline std::size_t parity_space_index() const { return parity_space_index_; }
 inline std::size_t N0_space_index() const { return N0_space_index_; }

 inline std::size_t element_offset(
     const std::size_t spin_index,
     const std::size_t x0_kappa0_index,
     const std::size_t Nbar_index
   ) const
  {
   return spin_index * bra_subspace().dimension() + x0_kappa0_index + Nbar_index;
  }


  inline std::size_t element_offset(
      const std::size_t spin_index,
      const u3::SU3& x0,
      const unsigned int kappa0,
      const unsigned int Nbar
    ) const
  {
    std::size_t x0_index = bra_subspace().LookUpSubspaceIndex(x0);
    std::size_t x0_kappa0_index =
        bra_subspace().GetSubspaceOffset(x0_index, kappa0);
    std::size_t Nbar_index =
        bra_subspace().GetSubspace(x0_index).LookUpStateIndex(Nbar);
    return element_offset(spin_index, x0_kappa0_index, Nbar_index);
  }

  inline KeyType Key() const
  {
    return {
        parity_space_index(),
        N0_space_index(),
        bra_subspace_index(),
        ket_subspace_index()
      };
  }

private:
  std::size_t parity_space_index_;
  std::size_t N0_space_index_;

};



template<typename tOperatorSpatialSpaceType, typename tOperatorSpatialSubspaceType,
 typename tOperatorSpinSpaceType,typename tOperatorSpinSubspaceType>
class OperatorU3SpinSectors
  : public basis::BaseSectors<
      tOperatorSpatialSpaceType,
      tOperatorSpinSpaceType,
      OperatorU3SpinSector<tOperatorSpatialSubspaceType,tOperatorSpinSubspaceType>
      >
{
private:
      // Constructor
    using BaseSectorsType = basis::BaseSectors<
      tOperatorSpatialSpaceType,
      tOperatorSpinSpaceType,
      OperatorU3SpinSector<tOperatorSpatialSubspaceType,tOperatorSpinSubspaceType>
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
      //By parity/exchange_symm
      //By S0
      //By N0
      //By L0
      for (uint8_t  parity : {0,1})
      {
        std::size_t parity_space_index = bra_space().LookUpSubspaceIndex(parity);
        if(parity_space_index == basis::kNone)
        {
          continue;
        }

        const auto& parity_space = bra_space().GetSubspace(parity_space_index);

        for (int i_spin = 0; i_spin < ket_space().size(); ++i_spin)
        {
          const auto& spin_subspace = ket_space().GetSubspace(i_spin);
          if((parity_space.parity_bar()+spin_subspace.exchange_symm_bar())%2 !=1)
            continue;

          const auto& S0 = spin_subspace.S0();
          for(int N0_space_index=0; N0_space_index<parity_space.size(); ++N0_space_index)
          {
            const auto& N0subspace = parity_space.GetSubspace(N0_space_index);
            for (int i_spatial = 0; i_spatial < N0subspace.size(); ++i_spatial)
            {
              const auto& L0subspace = N0subspace.GetSubspace(i_spatial);
              if (am::AllowedTriangle(S0, L0subspace.L0(), J0))
              {
                auto spatial_ptr = N0subspace.GetSubspacePtr(i_spatial);
                auto spin_ptr = ket_space().GetSubspacePtr(i_spin);

                BaseSectorsType::EmplaceSector(
                  i_spatial, i_spin,
                  spatial_ptr,spin_ptr,
                  parity_space_index,
                  N0_space_index
                  );
              }
            }
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
