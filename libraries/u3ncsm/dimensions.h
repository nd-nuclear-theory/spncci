/****************************************************************
  dimensions.h

  Functions used to get dimensions of lgi subspaces directly
  from the lsu3shell basis

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  11/1/21 (aem): Created.
****************************************************************/

#ifndef LGI_DIMENSION_H_
#define LGI_DIMENSION_H_

#include <vector>

#include "am/halfint.h"
#include "lgi/lgi.h"
#include "utilities/nuclide.h"

namespace lsu3shell
{
std::map<u3shell::U3SPN, unsigned int> LSU3ShellBasisDimensions(
    const nuclide::NuclideType& nuclide, const HalfInt& Nsigma0, const int Nmax
  );

std::map<u3shell::U3SPN, unsigned int> LSU3ShellCMFBasisDimensions(
    const HalfInt& Nsigma0,
    const int& Nmax,
    const std::map<u3shell::U3SPN, unsigned int>& u3spn_dimensions
  );


// renaming lsu3shell_cmf_basis_dimensions
std::map<u3shell::U3SPN, unsigned int> LSU3ShellCMFBasisDimensions(
    const nuclide::NuclideType& nuclide, const HalfInt& Nsigma0, const int& Nmax
  );
// Overload which generates u3spn_dimensions before getting cmf dimensions.

// ////TEMP
//  inline std::map<u3shell::U3SPN,unsigned int>
//  generate_lsu3shell_basis_dimensions(
//    const nuclide::NuclideType& nuclide,
//    const HalfInt& Nsigma0,
//    const int& Nmax
//  )
//  {
//    return LSU3ShellBasisDimensions(nuclide, Nsigma0,Nmax);
//  }

// print_lsu3shell_basis_info
void PrintLSU3ShellBasisInfo(
    const nuclide::NuclideType& nuclide, const int Nmax, const int Nstep
  );
// Function to print out lsu3shell basis.
// Nex x_proton S_proton x_neutron S_neutron x S dim

}  // namespace lsu3shell

namespace lgi
{
// renaming get_lgi_vector to

MultiplicityTaggedLGIVector GetLGIVector(
    const std::map<u3shell::U3SPN, unsigned int>& lsu3shell_cmf_dimensions,
    const HalfInt& Nsigma0,
    const unsigned int Nmax
  );

inline MultiplicityTaggedLGIVector GetLGIVector(
    const nuclide::NuclideType& nuclide,
    const HalfInt& Nsigma0,
    const unsigned int Nmax
  )
{
  auto lsu3shell_cmf_dimensions =
      lsu3shell::LSU3ShellCMFBasisDimensions(nuclide, Nsigma0, Nmax);

  return GetLGIVector(lsu3shell_cmf_dimensions, Nsigma0, Nmax);
}


}  // namespace lgi
#endif
