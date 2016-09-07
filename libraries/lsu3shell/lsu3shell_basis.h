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

  struct LSU3BasisGroupLabels
  // Data on the "multiplicity group" underlying an lsu3shell basis state.
  //
  // Fields:
  //   omegaSPN (u3shell::U3SPN) : U(3) and spin labels
  //   Np, Nn, Nex (int) : configuration excitation quanta
  {

    LSU3BasisGroupLabels(
        const u3shell::U3SPN& omegaSPN_,
        int Np_, int Nn_, int Nex_
      )
    : omegaSPN(omegaSPN_), Np(Np_), Nn(Nn_), Nex(Nex_)
    {}

    u3shell::U3SPN omegaSPN;
    int Np, Nn, Nex;
  };

  typedef std::vector<std::vector<LSU3BasisGroupLabels>> U3SPNBasisProvenance;
  // Container to hold lsu3shell basis provenance info (entries for each U3SPN
  // subspace, for each basis state).

  // typedef std::tuple<LSU3BasisGroupLabels,int,int> LSU3BasisTuple;
  // // Information on single lsu3shell multiplicity group.
  // //
  // // Fields:
  // //   omegaSPN (u3shell::U3SPN) -- U(3) and spin labels
  // //   Nexp, Nexn, Nex (int) -- configuration excitation quanta
  // //   dim -- group size
  // //   start_index -- starting index within U3SPN subspace

  struct LSU3BasisGroupData
    : LSU3BasisGroupLabels
  // Indexing information for lsu3shell multiplicity group.
  //
  // Fields:
  //   dim : group size
  //   start_index : starting index within U3SPN subspace
  {
    LSU3BasisGroupData (const LSU3BasisGroupLabels& lsu3_basis_group_labels_, int dim_, int start_index_)
      : LSU3BasisGroupLabels(lsu3_basis_group_labels_), dim(dim_), start_index(start_index_)
    {}

    int dim;
    int start_index;
  };

  typedef std::vector<LSU3BasisGroupData> LSU3BasisTable;
  // Container to hold basis group information for each basis group.

  
  void ReadLSU3Basis(
      HalfInt Nsigma_0, 
      const std::string& filename, 
      LSU3BasisTable& lsu3_basis_table,
      U3SPNBasisProvenance& basis_provenance,
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


}
#endif
