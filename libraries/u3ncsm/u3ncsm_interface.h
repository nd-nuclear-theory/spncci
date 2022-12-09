/****************************************************************
  seed_gen.h

  Computes seeds for spncci recurrence using lsu3shell

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  2/24/22 (aem): Created.
****************************************************************/

#ifndef U3NCSM_INTERFACE_H_
#define U3NCSM_INTERFACE_H_

#include <unordered_map>

#include "am/halfint.h"
#include "fmt/format.h"
#include "lgi/lgi.h"
// #include "lgi/recurrence_lgi.h"
#include "spncci_basis/recurrence_indexing.h"
#include "u3shell/relative_operator.h"
#include "utilities/nuclide.h"

namespace spncci::seeds
{

inline std::string seed_filename(
    const int Z,
    const int N,
    const HalfInt& Nsigma0,
    const u3::U3& sigma_bra,
    const u3::U3& sigma_ket,
    const int parity_bar
  )
{
  int Nex_bra(sigma_bra.N() - Nsigma0);
  int Nex_ket(sigma_ket.N() - Nsigma0);

  std::string filename = fmt::format(
      "seeds_Z{:02d}_N{:02d}_Nex{:02d}_lm{:02d}_mu{:02d}_Nex{:02d}_lm{:02d}_mu{"
      ":02d}_gbar{:1d}.dat",
      Z,
      N,
      Nex_bra,
      sigma_bra.SU3().lambda(),
      sigma_bra.SU3().mu(),
      Nex_ket,
      sigma_ket.SU3().lambda(),
      sigma_ket.SU3().mu(),
      parity_bar
    );

  return filename;
}


}  // namespace spncci::seeds

namespace lgi
{

inline std::string lgi_expansion_filename(const int Z, const int N, const lgi::LGI& lgi)
{
  const auto& [Nex, sigma, Sp, Sn, S] = lgi.Key();
  std::string filename = fmt::format(
      "lgi_expansion_Z{:02d}_N{:02d}_Nex{:02d}_lm{:02d}_mu{:02d}_2Sp{:02}_2Sn{:"
      "02}_2S{:02}.dat",
      Z,
      N,
      Nex,
      sigma.SU3().lambda(),
      sigma.SU3().mu(),
      TwiceValue(Sp),
      TwiceValue(Sn),
      TwiceValue(S)
    );
  return filename;
}

std::unordered_map<u3shell::U3SPN, basis::OperatorBlock<double>> GetLGIExpansions(
    const nuclide::NuclideType& nuclide,
    const spncci::spin::RecurrenceLGISpace<lgi::LGI, u3shell::spin::twobody::OperatorLabelsST>&
        lgi_space,
    const HalfInt& Nsigma0
  );
/// Reads in LGI expansions from file.

}  // namespace lgi


#endif
