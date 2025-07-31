/****************************************************************
  decomposition.cpp

  Generates decompositions for sp3rlib.

  Colin V. Coane
  
  SPDX-License-Identifier: MIT 

  07/08/24 (cvc): Created.

****************************************************************/


#include "sp3rlib/decomposition.h"

namespace sp3r
{
  //
  void CalculateNexDecompositions(
    const std::vector<sp3r::Sp3RSJSpace>& spaces_sp3rsj,
    const std::vector<sp3r::Matrix>& eigenvectors,
    std::vector<sp3r::Matrix>& Nex_decompositions,
    HalfInt Nsigma0,
    int Nmax
  )
  {
    // loop through Sp3RSJ spaces
    for (int sp3rsj_space_index=0; sp3rsj_space_index<spaces_sp3rsj.size(); sp3rsj_space_index++)
    {
      int offset = 0;
      // alias for current S-J space
      const sp3r::Sp3RSJSpace sp3rsj_space = spaces_sp3rsj[sp3rsj_space_index];

      const sp3r::Matrix& eigenvectors_SJ = eigenvectors[sp3rsj_space_index];
      sp3r::Matrix& Nex_decompositions_SJ = Nex_decompositions[sp3rsj_space_index];

      // Initialize decomposition matrix
      const int decomposition_size = Nmax + 1;
      const int num_eigenvectors = eigenvectors_SJ.cols();
      Nex_decompositions_SJ = sp3r::Matrix::Zero(decomposition_size,num_eigenvectors);

      // Accumulate probability
      // Loop over subspaces
      for (int sp3rsj_subspace_index=0; sp3rsj_subspace_index<sp3rsj_space.size(); sp3rsj_subspace_index++)
      {
        // Retrieve subspace information
        const sp3r::U3Subspace& u3_subspace = sp3rsj_space.GetSubspace(sp3rsj_subspace_index);
        int upsilon_max = u3_subspace.upsilon_max();
        // Retrieve basis state information
        for (int sp3rsj_state_index=0; sp3rsj_state_index<u3_subspace.size(); sp3rsj_state_index++)
        {
          // Retrieve state information
          //int state_offset = u3_subspace.GetStateOffset(sp3rsj_state_index);
          sp3r::SO3State so3_state(u3_subspace,sp3rsj_state_index);
          int kappa_max = so3_state.kappa_max();
          int state_degeneracy = kappa_max*upsilon_max;

        }
      }
    }
  }


      // accumulate probability
      for (int spj_subspace_index=0; spj_subspace_index<spj_space.size(); ++spj_subspace_index)
        // for each (composite) state
        {
          // retrieve basis state information          
          const spncci::SubspaceSpBasis& spj_subspace=spj_space.GetSubspace(spj_subspace_index);
          for(int spj_state_index=0; spj_state_index<spj_subspace.size(); ++spj_state_index)
            {
              StateSpBasis spj_state(spj_subspace,spj_state_index);    
              int degeneracy = spj_state.degeneracy();
              int Nex = int(spj_state.omega().N()-Nsigma0);
              assert((0<=Nex)&&(Nex<=Nmax));
 
              // accumulate probability from this (composite) state
               Nex_decompositions_J.row(Nex)+=eigenvectors_J.block(offset,0,degeneracy,num_eigenvectors).colwise().squaredNorm();
              offset+=degeneracy;
            }
          

        }


}

