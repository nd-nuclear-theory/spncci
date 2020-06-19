/****************************************************************
  lsu3shell_basis.h

  Input of lsu3shell basis listing.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  8/1/16 (aem,mac): Created.
  9/7/16 (mac): Split from lsu3shell_interface.
  6/26/17 (mac): Add sanity checks on input basis labels.
  9/5/17 (aem) : Renamed LSU3Basis to LSU3ShellBasis
****************************************************************/

#ifndef LSU3SHELL_BASIS_H_
#define LSU3SHELL_BASIS_H_

#include "u3shell/u3spn_scheme.h"

namespace lsu3shell
{

  ////////////////////////////////////////////////////////////////
  // data structures for lsu3shell basis labels
  ////////////////////////////////////////////////////////////////

  struct LSU3ShellBasisGroupLabels
  // Data on the "multiplicity group" underlying an lsu3shell basis state.
  //
  // Fields:
  //   omegaSPN (u3shell::U3SPN) : U(3) and spin labels
  //   ip, in (int) : configuration indices
  //   Np, Nn, Nex (int) : configuration excitation quanta
  {

    LSU3ShellBasisGroupLabels(
        const u3shell::U3SPN& omegaSPN_,
        int ip_, int in_, int Np_, int Nn_, int Nex_
      )
    : omegaSPN(omegaSPN_), ip(ip_), in(in_), Np(Np_), Nn(Nn_), Nex(Nex_)
    {}

    u3shell::U3SPN omegaSPN;
    int ip, in, Np, Nn, Nex;
  };

  typedef std::vector<std::vector<LSU3ShellBasisGroupLabels>> U3SPNBasisLSU3Labels;
  // Container to hold lsu3shell basis provenance info (entries for each U3SPN
  // subspace, for each basis state).

  struct LSU3ShellBasisGroupData
    : LSU3ShellBasisGroupLabels
  // Indexing information for lsu3shell multiplicity group.
  //
  // TODO: This inheritance LSU3ShellBasisGroupData : LSU3ShellBasisGroupLabels
  // violates the "is a type of" rule for inheritance.  Consider
  // restructuring/abolishing.
  //
  // Fields:
  //   dim : group size
  //   start_index : starting index within U3SPN subspace
  {
    LSU3ShellBasisGroupData (const LSU3ShellBasisGroupLabels& lsu3_basis_group_labels_, int dim_, int start_index_)
      : LSU3ShellBasisGroupLabels(lsu3_basis_group_labels_), dim(dim_), start_index(start_index_)
    {}

    int dim;
    int start_index;
  };

  typedef std::vector<LSU3ShellBasisGroupData> LSU3ShellBasisTable;
  // Container to hold basis group information for each basis group.

  ////////////////////////////////////////////////////////////////
  // lsu3shell basis input
  ////////////////////////////////////////////////////////////////
  
  void ReadLSU3ShellBasis(
      HalfInt Nsigma_0, 
      const std::string& filename, 
      LSU3ShellBasisTable& lsu3_basis_table,
      U3SPNBasisLSU3Labels& basis_provenance,
      u3shell::SpaceU3SPN& space
    );
  // Read lsu3shell basis multiplicity group table and set up
  // resulting mapping to U3SPN basis.
  //
  // Arguments:
  //   Nsigma_0 (HalfInt) : Base excitation to use in calculating
  //     total oscillator quanta
  //   filename (std::string) : Input filename
  //   ...
  //   space (output) : lsu3shell basis space 
  //   basis_provenance (output) : look up table between 

}
#endif
