/****************************************************************
recurrence_seeds.h

    Populating recurrence seed matricex 


  Anna E. McCoy[1] and Patrick J. Fasano[2,3]
  [1] Institute for Nuclear Theory
  [2] University of Notre Dame
  [3] Lawrence Berkeley National Laboratory

  SPDX-License-Identifier: MIT

  09/24/21 (aem) : Created.
****************************************************************/

#ifndef RECURRENCE_SEEDS_H_
#define RECURRENCE_SEEDS_H_

#include <array>
#include <functional>  // for std::hash
#include <map>
#include <tuple>
#include <utility>

#include "am/am.h"
#include "am/halfint.h"
#include "basis/basis.h"
#include "basis/degenerate.h"
#include "u3shell/operator_indexing_spin.h"
#include "lgi/lgi.h"
#include "sp3rlib/multiplicity_tagged.h"
#include "sp3rlib/u3.h"
#include "spncci_basis/recurrence_indexing.h"
#include "lgi/recurrence_lgi.h"
#include "u3ncsm/u3ncsm_interface.h"

namespace spncci::seeds
{

basis::OperatorBlocks<double>
GetRecurrenceSeedsFromFile(
  const lgi::MultiplicityTaggedLGIVector& lgi_vector,
  const std::vector<int>& lgi_full_space_index_lookup,
  const spncci::spatial::RecurrenceSpace<u3shell::spatial::OneCoordType>& spatial_recurrence_space,
  const spncci::spin::RecurrenceSpace<lgi::LGI, u3shell::spin::twobody::OperatorLabelsST>& spin_recurrence_space
);


basis::OperatorBlocks<double>
GetRecurrenceSeedsFromFile(
  const nuclide::NuclideType& nuclide,
  const spncci::spatial::RecurrenceSpace<u3shell::spatial::OneCoordType>& spatial_recurrence_space,
  const spncci::spin::RecurrenceSpace<lgi::LGI, u3shell::spin::twobody::OperatorLabelsST>& spin_recurrence_space,
  const std::string& seed_filename_template
)
{
  const auto& [Z,N] = nuclide;
  bool intrinsic = true;
  HalfInt Nsigma0 = nuclide::Nsigma0ForNuclide(nuclide,intrinsic);
  for(auto recurrence_sp3r_space_index=0; recurrence_sp3r_space_index < spatial_recurrence_space.size(); ++recurrence_sp3r_space_index)
    {
      // Get RecurrenceSp3RSpace and extract labels
      const auto& recurrence_sp3r_space = spatial_recurrence_space.GetSubspace(recurrence_sp3r_space_index);
      const auto& sigma_ket = recurrence_sp3r_space.sigma_ket();
      const auto& sigma_bra = recurrence_sp3r_space.sigma_bra();
      const auto parity_bar = recurrence_sp3r_space.parity_bar();

      // Look up corresponding RecurrenceLGISpace
      const auto exchange_symm_bar = (parity_bar+1)%2;
      const auto& spin_recurrence_lgi_space = spin_recurrence_space.LookUpSubspace({sigma_ket,sigma_bra,exchange_symm_bar});

      // using default template for seed filename
      std::string seed_filename = spncci::seeds::seed_filename(
        Z,N,Nsigma0,sigma_bra,sigma_ket,parity_bar
        );

      auto rows = recurrence_sp3r_space.dimension();
      auto cols = spin_recurrence_lgi_space.dimension();
      // TODO: allocate and read seeds

    }

}





}


#endif  // RECURRENCE_SEEDS_H_
