/****************************************************************
  decomposition.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/decomposition.h"


namespace spncci
{
  // Note on decomposition algorithm:
  //
  // Accumulation of probabilities is done for all the
  // eigenvectors in the subspace in parallel, by means of
  // "column-wise partial reduction operations".  See Eigen
  // documentation:
  //
  //   man/eigen/eigen-doc-3.2.8/group__TutorialReductionsVisitorsBroadcasting.html
  //
  // Note that column-wise partial reductions return a row vector.
  //
  // We wish to accumulate the probability contributions from a
  // group of rows (in the eigenvector matrix) corresponding to
  // the "substates" of a single basis "state", in the
  // terminology of basis/degenerate.h.  (These are marked by
  // the vertical side bar in the diagram below.)  We thus
  // calculate the .columnwise().squaredNorm() of these rows of
  // the eigenvector matrix.  The result is a row vector of the
  // summed squares of the amplitudes from these rows.  We must
  // accumulate these onto the appropriate row of probabilities
  // in the decomposition matrix.
  //
  // Source (eigenvectors_J):
  //   
  //                     (v1)_0     (v2)_0     (v3)_0
  //                     (v1)_1     (v2)_1     (v3)_1
  //                     ...
  //   offset ->       | (v1)_n     (v2)_n     (v3)_n
  //                /  | (v1)_n+1   (v2)_n+1   (v3)_n+1
  //   (degeneracy)-   | ...      
  //                \  | (v1)_n+d-1 (v2)_n+d-1 (v3)_n+d-1
  //                     ...
  //
  // Result of column-wise partial reduction:
  //
  //                     sum1       sum2       sum3
  //
  // Target (xxx_decompositions_J):
  //
  //                     P(v1;0)   P(v2;0)  P(v3;0)
  //                     P(v1;1)   P(v2;1)  P(v3;1)
  //                     ...
  //   accumulate to ->  P(v1;m)   P(v2;m)  P(v3;m)
  //                     ...

  void CalculateNexDecompositions(
      const spncci::SpaceSpJ& spj_space,
      const std::vector<spncci::Matrix>& eigenvectors,
      std::vector<spncci::Matrix>& Nex_decompositions,
      HalfInt Nsigma0,
      int Nmax
    )
  {
    for (int spj_subspace_index=0; spj_subspace_index<spj_space.size(); ++spj_subspace_index)
      // for each J subspace
      {
        // set up aliases for current J subspace
        const SubspaceSpJ& spj_subspace = spj_space.GetSubspace(spj_subspace_index);
        const spncci::Matrix& eigenvectors_J = eigenvectors[spj_subspace_index];
        spncci::Matrix& Nex_decompositions_J = Nex_decompositions[spj_subspace_index];

        // initialize decomposition matrix
        const int decomposition_size = Nmax+1;
        const int num_eigenvectors = eigenvectors_J.cols();
        Nex_decompositions_J = spncci::Matrix::Zero(decomposition_size,num_eigenvectors);

        // accumulate probability
        for (int spj_state_index=0; spj_state_index<spj_subspace.size(); ++spj_state_index)
          // for each (composite) state
          {
            // retrieve basis state information
            StateSpJ spj_state(spj_subspace,spj_state_index);
            int offset = spj_state.offset();
            int degeneracy = spj_state.degeneracy();
            int Nex = int(spj_state.N()-Nsigma0);
            assert((0<=Nex)&&(Nex<=Nmax));

            // accumulate probability from this (composite) state
            Nex_decompositions_J.row(Nex) += eigenvectors_J.block(offset,0,degeneracy,num_eigenvectors).colwise().squaredNorm();
          }
            
          
      }
  }

  void CalculateBabySpNCCIDecompositions(
      const spncci::SpaceSpJ& spj_space,
      const std::vector<spncci::Matrix>& eigenvectors,
      std::vector<spncci::Matrix>& baby_spncci_decompositions,
      int baby_spncci_space_size
    )
  {
    for (int spj_subspace_index=0; spj_subspace_index<spj_space.size(); ++spj_subspace_index)
      // for each J subspace
      {
        // set up aliases for current J subspace
        const SubspaceSpJ& spj_subspace = spj_space.GetSubspace(spj_subspace_index);
        const spncci::Matrix& eigenvectors_J = eigenvectors[spj_subspace_index];
        spncci::Matrix& baby_spncci_decompositions_J = baby_spncci_decompositions[spj_subspace_index];

        // initialize decomposition matrix
        const int decomposition_size = baby_spncci_space_size;
        const int num_eigenvectors = eigenvectors_J.cols();
        baby_spncci_decompositions_J = spncci::Matrix::Zero(decomposition_size,num_eigenvectors);

        // accumulate probability
        for (int spj_state_index=0; spj_state_index<spj_subspace.size(); ++spj_state_index)
          // for each (composite) state
          {
            // retrieve basis state information
            StateSpJ spj_state(spj_subspace,spj_state_index);
            int offset = spj_state.offset();
            int degeneracy = spj_state.degeneracy();
            int baby_spncci_subspace_index = spj_state.baby_spncci_subspace_index();
            assert((0<=baby_spncci_subspace_index)&&(baby_spncci_subspace_index<baby_spncci_space_size));

            // accumulate probability from this (composite) state
            baby_spncci_decompositions_J.row(baby_spncci_subspace_index) += eigenvectors_J.block(offset,0,degeneracy,num_eigenvectors).colwise().squaredNorm();
          }
            
          
      }
  }


}  // namespace
