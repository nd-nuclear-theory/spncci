/****************************************************************
  seed_gen.h

  Computes seeds for spncci recurrence using lsu3shell

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  2/24/22 (aem): Created.
****************************************************************/

#ifndef SEED_GEN_H_
#define SEED_GEN_H_

#include <unordered_map>
#include "am/halfint.h"
#include "fmt/format.h"

#include "lgi/lgi.h"
#include "lgi/recurrence_lgi.h"
#include "spncci_basis/recurrence_indexing.h"
#include "u3shell/relative_operator.h"
#include "utilities/nuclide.h"

#include "LSU3/ncsmSU3xSU2Basis.h"

namespace spncci{
namespace seeds{

  inline std::string seed_filename(
    const int Z, const int N,
    const HalfInt& Nsigma0,
    const u3::U3& sigma_bra,
    const u3::U3& sigma_ket,
    const int parity_bar
    )
    {
      int Nex_bra(sigma_bra.N()-Nsigma0);
      int Nex_ket(sigma_ket.N()-Nsigma0);

      std::string filename
        = fmt::format("seeds_Z{:02d}_N{:02d}_Nex{:02d}_lm{:02d}_mu{:02d}_Nex{:02d}_lm{:02d}_mu{:02d}_gbar{:1d}.dat",
              Z,N,Nex_bra,sigma_bra.SU3().lambda(),sigma_bra.SU3().mu(),
              Nex_ket,sigma_ket.SU3().lambda(),sigma_ket.SU3().mu(),parity_bar
            );

      return filename;
    }

    std::unordered_map<u3shell::U3SPN,basis::OperatorBlock<double>>
    GetLGIExpansions(
      const nuclide::NuclideType& nuclide,
      const spncci::spin::RecurrenceLGISpace<lgi::LGI,spncci::spin::UnitTensorLabelsST>& lgi_space,
      const HalfInt& Nsigma0
      );


    basis::OperatorBlock<double> GenerateRecurrenceSeedBlock(
      const nuclide::NuclideType& nuclide,
      const HalfInt& Nsigma0,
      const int N1v,
      const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& unit_tensor_labels,
      const spncci::spin::RecurrenceSpace<lgi::LGI, spncci::spin::UnitTensorLabelsST>& spin_recurrence_space,
      const spncci::spatial::RecurrenceSpace& spatial_recurrence_space,
      const CBaseSU3Irreps& baseSU3Irreps,
      const std::string& operator_dir,
      const int recurrence_sp3r_space_index
      );



}} // namespace
#endif
