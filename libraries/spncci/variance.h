/****************************************************************
  variance.h

  Compute variance of operators in spncci basis
                                  
  Anna E. McCoy
  TRIUMF

  SPDX-License-Identifier: MIT

  6/17/19 (aem): Created.
****************************************************************/

#ifndef SPNCCI_VARIANCE_H_
#define SPNCCI_VARIANCE_H_

#include "spncci/spncci_basis.h"
#include "spncci/parameters.h"
// Functions used to calculate variance of operator matrix when the 
// space is broken up into a reference space H and the remainder of the full space V


namespace spncci
{

  void GetVarianceBlock(
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const u3shell::ObservableSpaceU3S& observable_space,
    const HalfInt& J0,
    const spncci::SpaceSpBasis& spbasis_bra, //For a given J
    const spncci::SpaceSpBasis& spbasis_ket, //For a given J
    const std::vector<spncci::LGIPair>& lgi_pairs, //Defines tiles to get computed 
    int observable_index, int hw_index,
    spncci::OperatorBlock& operator_matrix
  );
  // For give lgi pairs, compute operator tiles and accumulate in operator matrix


  void CalculateVariance(
    const spncci::OperatorBlock& eigenvectors,
    const spncci::OperatorBlock& operator_matrix,
    spncci::OperatorBlock& variance_block
  );
  // Given the operator matrix with rows corresponding to H space and columns corresponding to V spaces
  // and eigenvectors in H space, compute variance given by <psi|H12*H21|psi>


  void GetEigensystemH(
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const u3shell::ObservableSpaceU3S& observable_space,
    int hw_index, const spncci::RunParameters& run_parameters,
    const std::set<int>& irrep_families_H,
    const std::vector<spncci::SpaceSpBasis>& spbasis_H_byJ,
    std::vector<spncci::Vector>& eigenvalues,  // eigenvalues by J subspace
    std::vector<spncci::Matrix>& eigenvectors  // eigenvectors by J subspace
  );
  // Calculate and diagonalize Hamiltonian matrix in H space

  void GetVariances(
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const u3shell::ObservableSpaceU3S& observable_space,
    int observable_index, int hw_index, const HalfInt& J0,
    const std::vector<std::pair<int,int>>& sectors_J,
    const spncci::RunParameters& run_parameters,
    const std::set<int>& irrep_families_H,
    const std::set<int>& irrep_families_V,
    const std::vector<spncci::SpaceSpBasis>& spbasis_H_byJ,
    std::vector<spncci::Matrix>& eigenvectors,  // eigenvectors by J subspace
    std::vector<std::vector<double>>& variances
  );
  // Given H and V spaces and eigenvectors obtained in the H space, get variances for given
  // observable for each valueof J.  


void GetVariancesForIrrepFamilies(
    const std::vector<int>& irrep_families,
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const u3shell::ObservableSpaceU3S& observable_space,
    int observable_index, int hw_index,
    const spncci::RunParameters& run_parameters,
    const std::set<int>& irrep_families_H,
    std::vector<std::vector<std::vector<double>>>& variance_table,
    std::vector<int>& list_irrep_families_V
  );
  // Iterate over all irrep families not in H and compute the variance in the new space
  // formed by H+irrep_family.  Varience stored in variance table index by index of 
  // irrep_family in list_irrep_families_V then by J index.
  //
  // inputs: irrep_families, baby_spncci_space, observable_space, observable_index, hw_index,
  //          run_paramters, irrep_families_H
  //
  // output: variance_table, list_irrep_families_V


void SortIrrepFamiliesByVariance(
    const std::vector<std::vector<std::vector<double>>>& variances,
    const std::vector<int>& irrep_families_V,
    int J_index, int eigenvalue_index,
    std::vector<int>& irrep_families_by_variance,
    double variance_threshold=1e-8
  );
  // Orders irrep families from largest to smallest variance,
  //  input: variances, irrep_families_V, J_index, eigenvalue_index
  //  output: irrep_families_by_variance
void SortIrrepFamiliesByNex(
  const lgi::MultiplicityTaggedLGIVector& lgi_families,
  std::vector<std::vector<int>>& irrep_families_by_Nex,
  int Nmax
  );

void TestingVariances(
    const spncci::RunParameters& run_parameters, 
    int hw_index, int J_index, int eigenvalue_index,
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const std::vector<u3shell::ObservableSpaceU3S>& observable_spaces,
    const std::vector<int>& irrep_families,
    const std::set<int>& reference_H,
    const lgi::MultiplicityTaggedLGIVector& lgi_families,
    double variance_threshold
  );


}

#endif