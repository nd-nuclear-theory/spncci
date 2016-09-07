/****************************************************************
  lsu3shell_basis.h

  Input of lsu3shell basis listing.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  8/1/16 (aem,mac): Created.
  9/7/16 (mac): Split from lsu3shell_interface.

****************************************************************/

#ifndef LSU3SHELL_BASIS_H_
#define LSU3SHELL_BASIS_H_

#include "u3shell/u3spn_scheme.h"

namespace lsu3shell
{

  // struct LSU3BasisGroup
  // // Data on a "multiplicity group" within the lsu3shell basis.
  // //
  // // Much of the input info is merely for diagnostic purposes on the
  // // lsu3shell configuration, which does not actually enter into our
  // // subsequent calculations.
  // //
  // // Fields:
  // //   omegaSPN (u3shell::U3SPN) -- U(3) and spin labels
  // //   Nexp, Nexn, Nex (int) -- configuration excitation quanta
  // //   dim -- group size
  // //   start_index -- starting index within U3SPN subspace
  // {
  // }

  typedef std::tuple<u3shell::U3SPN,int,int> LSU3BasisGroup;
  typedef std::vector<LSU3BasisGroup> LSU3BasisTable;

  void ReadLSU3Basis(
      HalfInt Nsigma_0, 
      const std::string& filename, 
      LSU3BasisTable& lsu3_basis_table, 
      std::map<u3shell::U3SPN,int>& subspace_dimensions
    );


}
#endif
