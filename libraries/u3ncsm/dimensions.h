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
#include "lgi/lgi.h"
#include "utilities/nuclide.h"
#include "am/halfint.h"

namespace lsu3shell
{
  std::map<u3shell::U3SPN,unsigned int> 
  lsu3shell_basis_dimensions(
    const nuclide::NuclideType& nuclide, 
    const HalfInt& Nsigma0,
    const int& Nmax
  );

  std::map<u3shell::U3SPN, unsigned int>
  lsu3shell_cmf_basis_dimensions(
    const HalfInt& Nsigma0,
    const int& Nmax, 
    const std::map<u3shell::U3SPN, unsigned int>& u3spn_dimensions
    );


  std::map<u3shell::U3SPN, unsigned int>
  lsu3shell_cmf_basis_dimensions(
    const nuclide::NuclideType& nuclide, 
    const HalfInt& Nsigma0,
    const int& Nmax
  );
  // Overload which generates u3spn_dimensions before getting cmf dimensions.
 ////TEMP
  inline std::map<u3shell::U3SPN,unsigned int> 
  generate_lsu3shell_basis_dimensions(
    const nuclide::NuclideType& nuclide, 
    const HalfInt& Nsigma0,
    const int& Nmax
  )
  {
    return lsu3shell_basis_dimensions(nuclide, Nsigma0,Nmax);
  }

  void print_lsu3shell_basis_info(
    const nuclide::NuclideType& nuclide,
    const int Nmax,
    const int Nstep
  );
  // Function to print out lsu3shell basis.
  // Nex x_proton S_proton x_neutron S_neutron x S dim

}//namespace lsu3shell

namespace lgi
{
  MultiplicityTaggedLGIVector get_lgi_vector(
      const nuclide::NuclideType& nuclide, 
      const HalfInt& Nsigma0,
      const unsigned int& Nmax
    );



}
#endif
