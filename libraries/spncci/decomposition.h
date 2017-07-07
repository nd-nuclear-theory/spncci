/****************************************************************
  decomposition.h

  Code to generate eigenfunction decompositions.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  6/22/17 (mac): Created.
****************************************************************/

#ifndef SPNCCI_SPNCCI_DECOMPOSITION_H_
#define SPNCCI_SPNCCI_DECOMPOSITION_H_

#include "spncci/branching.h"

namespace spncci
{

  void CalculateNexDecompositions(
      const spncci::SpaceSpJ& spj_space,
      const std::vector<spncci::Matrix>& eigenvectors,
      std::vector<spncci::Matrix>& Nex_decompositions,
      HalfInt Nsigma0,
      int Nmax
    );
  // Calculate decompositions of eigenstates w.r.t. Nex.
  //
  // This calculation is redundant to the decomposition
  // w.r.t. BabySpNCCI subspaces, since each BabySpNCCI subspace has
  // good Nex.  The Nex decomposition is just a coarser binning of the
  // BabySpNCCI decompsition.  However, we keep the calculation of the
  // Nex decomposition as a separate stand-alone routine since this
  // decomposition provides a basic, immediate diagnostic, also
  // provided by M-scheme NCCI codes.
  //
  // Both the eigenvectors and the decompositions are stored as a
  // vector of matrices:
  //
  //   J_subspace_index -> matrix(basis_index,eigenvector_index)
  //
  // That is, each J subspace gets its own matrix, and each column of
  // the matrix stores the wave function amplitudes or decomposition
  // probabilities for a single eigenvector.
  //
  // Arguments:
  //   spj_space (input): final branched basis
  //   eigenvectors (input): eigenvectors in final branched basis
  //   Nex_decompositions (output): decompositions of eigenvectors
  //   Nsigma0 (input): base U(1) quantum number for calculation
  //     (used to convert N_omega to Nex)
  //   Nmax (input): highest Nex in basis (used to size the
  //     decomposition matrix)

  void CalculateBabySpNCCIDecompositions(
      const spncci::SpaceSpJ& spj_space,
      const std::vector<spncci::Matrix>& eigenvectors,
      std::vector<spncci::Matrix>& baby_spncci_decompositions,
      int baby_spncci_space_size
    );
  // Calculate decompositions of eigenstates w.r.t. BabySpNCCI
  // subspace.
  //
  // Both the eigenvectors and the decompositions are stored as a
  // vector of matrices:
  //
  //   J_subspace_index -> matrix(basis_index,eigenvector_index)
  //
  // That is, each J subspace gets its own matrix, and each column of
  // the matrix stores the wave function amplitudes or decomposition
  // probabilities for a single eigenvector.
  //
  // Arguments:
  //   spj_space (input): final branched basis
  //   eigenvectors (input): eigenvectors in final branched basis
  //   baby_spncci_decompositions (output): decompositions of eigenvectors
  //   baby_spncci_space_size (input): number of BabySpNCCI subspaces
  //     (used to size the decomposition matrix)

}  // namespace

#endif
