/****************************************************************
  sp3rsj_tiling.cpp

  Colin V. Coane
  
  SPDX-License-Identifier: MIT 

  02/05/24 (cvc): Created.
  05/23/24 (cvc): Updated.
  06/03/24 (cvc): Added constructor for Sp3RSJ hyperblocks.
  06/03/24 (cvc): Refactored to .h/.cpp format.
  06/11/24 (cvc): Fixed looping over hypersectors.
  07/05/24 (cvc): Observables changed to linear combinations of operators.
  07/08/24 (cvc): ConstructSp3RHyperBlocks() computes Sp(3,R)xSU(2)>SU(2) RMEs.

****************************************************************/

#include "sp3rlib/sp3rsj_tiling.h"

namespace sp3r
{
  
// Construct hyperblocks
void ConstructSp3RHyperBlocks(
  const sp3r::Sp3RSpace& sp3r_space,
  const sp3r::Sp3RSJSpace& bra_space,
  const sp3r::Sp3RSJSpace& ket_space,
  const sp3r::OperatorSpaceSp3RSJ& op_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  const int bra_index,
  const int ket_index,
  HalfInt& J_bra, HalfInt& J_ket, HalfInt& J_op,
  HalfInt& S_bra, HalfInt& S_ket, HalfInt& S_op,
  u3::UCoefCache& u_coef_cache,
  u3::WCoefCache& w_coef_cache,
  const sp3r::OperatorInfoSp3RSJ& operator_info,
  sp3r::Sp3RSJHypersectors& hypersectors,
  basis::OperatorHyperblocks<double>& operator_hyperblocks
)
{
  // Get sigma
  const u3::U3 sigma = sp3r_space.sigma();
  // Get rho for coupling of operators: TODO FIX THIS IF NEEDED????
  int rho0 = 1;
  // Loop over operator space
  for (int op_index=0; op_index<op_space.size(); op_index++)
  {
    // Get operator subspace
    const OperatorU3Subspace& op_subspace = op_space.GetSubspace(op_index);
    const u3::U3& omega_bra = bra_subspace.U3();
    const u3::U3& omega_ket = ket_subspace.U3();
    const u3::U3& omega_op = op_subspace.U3();
    unsigned int multiplicity = u3::OuterMultiplicity(omega_ket, omega_op, omega_bra);
    // ASK ABOUT MULTIPLICITY INDEX - DO WE NEED IT???
    // Loop over multiplicities
    for (int mult_idx = 1; mult_idx <= multiplicity; mult_idx++)
    {
      // Look up corresponding hypersector index
      int hypersector_index = hypersectors.LookUpHypersectorIndex(bra_index,ket_index,op_index,mult_idx);

      // If not in hypersectors, continue
      if (hypersector_index==-1)
        continue;

      // Get observable info
      // Populate hyperblock with operator content
      operator_hyperblocks[hypersector_index][0] = operator_info.ComputeRMESector(omega_op,mult_idx,sp3r_space,bra_subspace,ket_subspace,u_coef_cache,S_bra,S_ket,J_bra,J_ket);
    }
    // End looping over multiplicity
  }
  // End looping over operator space
}


// Get operator tile
void GetSp3ROperatorTile(
  const sp3r::Sp3RSpace& sp3r_space,
  const sp3r::Sp3RSJSpace& bra_space,
  const sp3r::Sp3RSJSpace& ket_space,
  const sp3r::OperatorSpaceSp3RSJ& op_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  const int bra_index,
  const int ket_index,
  HalfInt& J_bra, HalfInt& J_ket, HalfInt& J_op,
  HalfInt& S_bra, HalfInt& S_ket, HalfInt& S_op,
  u3::UCoefCache& u_coef_cache,
  u3::WCoefCache& w_coef_cache,
  const sp3r::OperatorInfoSp3RSJ& operator_info,
  sp3r::Sp3RSJHypersectors& hypersectors,
  basis::OperatorHyperblocks<double>& operator_hyperblocks,
  basis::OperatorBlock<double>& tile
)
{
  // Subspave labels and info
  const u3::U3& omega_bra = bra_subspace.U3();
  const u3::U3& omega_ket = ket_subspace.U3();

  // Get tile dimensions
  // Tile dimension = upsilon_max*(total states in U3 irrep)
  const int tile_dimension_bra = bra_subspace.dimension()*bra_subspace.total_states();
  const int tile_dimension_ket = ket_subspace.dimension()*ket_subspace.total_states();

  // Construct tile
  tile = basis::OperatorBlock<double>::Zero(tile_dimension_bra,tile_dimension_ket);

  // Branch hypersectors to J and add to tile
  // Loop over operator space
  for (int op_index=0; op_index<op_space.size(); op_index++)
  {
    // Get operator subspace
    const OperatorU3Subspace& op_subspace = op_space.GetSubspace(op_index);
    const u3::U3& omega_op = op_subspace.U3();
    unsigned int multiplicity = u3::OuterMultiplicity(omega_ket, omega_op, omega_bra);
    // Look up corresponding hypersector index
    // Loop over multiplicities
    for (int mult_idx = 1; mult_idx <= multiplicity; mult_idx++)
    {
      int hypersector_index = hypersectors.LookUpHypersectorIndex(bra_index,ket_index,op_index,mult_idx);

      // If not in hypersectors, continue
      if (hypersector_index==-1)
        continue;
      
      // Get unbranched block
      basis::OperatorBlock<double> block = operator_hyperblocks[hypersector_index][0];
      int block_rows=block.rows();
      int block_cols=block.cols();

      // Get list of L_bra, L_ket values for this sector
      MultiplicityTagged<unsigned int>::vector L_bra_list = u3::BranchingSO3(omega_bra.SU3());
      MultiplicityTagged<unsigned int>::vector L_ket_list = u3::BranchingSO3(omega_ket.SU3());
      // Get allowed L_op values and multiplicities
      MultiplicityTagged<unsigned int>::vector L_op_list = op_subspace.GetOpStateLabels();
      
      // Loop over bra states
      for (const auto& [L_bra, kappa_max_bra] : L_bra_list)
      {
        // Check if bra state exists
        if (not bra_subspace.ContainsState(L_bra))
          continue;

        for (const auto& [L_ket, kappa_max_ket] : L_ket_list)
        {
          // Check if ket state exists
          if (not ket_subspace.ContainsState(L_ket))
            continue;
          
          // Get L_bra, L_ket indices
          int state_index_bra = bra_subspace.LookUpStateIndex(L_bra);
          int state_index_ket = ket_subspace.LookUpStateIndex(L_ket);
          // Loop over operator states
          for (const auto& [L_op, kappa_max_op] : L_op_list)
          {
            // Get LSJ 9J Coefficient
            double LSJcoef = am::Unitary9J(L_ket,S_ket,J_ket,L_op,S_op,J_op,L_bra,S_bra,J_bra);

            // Get offsets in subspace
            int offset_bra = bra_subspace.GetStateOffset(state_index_bra);
            int offset_ket = ket_subspace.GetStateOffset(state_index_ket);

            // Get global indices in matrix
            // index = offset*(block dimension)
            int index_bra = offset_bra*block_rows;
            // Loop over kappa for bra state
            for (int kappa_bra = 1; kappa_bra <= kappa_max_bra; kappa_bra++)
            {
              // Get global indices in matrix
              // index = offset*(block dimension)
              int index_ket = offset_ket*block_cols;
              // Loop over kappa for ket state
              for (int kappa_ket = 1; kappa_ket <= kappa_max_ket; kappa_ket++)
              {
                // Loop over kappa for operator state
                for (int kappa_op = 1; kappa_op <= kappa_max_op; kappa_op++)
                {
                  // W coefficient
                  double Wcoef = u3::WCached(
                    w_coef_cache,omega_ket.SU3(),kappa_ket,L_ket,
                    omega_op.SU3(),kappa_op,L_op,
                    omega_bra.SU3(),kappa_bra,L_bra,mult_idx
                  );

                  // Accumulate branched block in tile
                  tile.block(index_bra,index_ket,block_rows,block_cols) += Wcoef*LSJcoef*block;
                }
                // End looping over kappa_op
                // Increment starting index for ket
                index_ket += block_cols;
              }
              // end looping over kappa_ket
              // Increment starting index for bra
              index_bra += block_rows;
            }
            // End looping over kappa_bra
          }
          // End looping over L_op
        }
        // End looping over L_ket
      }
      // End looping over L_bra
    }
    // End looping over multiplicities
  }
  // End looping over operator subspaces
}


// Construct Operator Matrix
void ConstructSp3ROperatorMatrix(
  const sp3r::Sp3RSpace& sp3r_space,
  const sp3r::Sp3RSJSpace& bra_space,
  const sp3r::Sp3RSJSpace& ket_space,
  const sp3r::ObservableSp3RSJ& observable_list,
  basis::OperatorBlock<double>& operator_matrix
)
{
  // Extract info from base space and bra/ket spaces (dims, labels, etc.)

  // Ensure same Sp(3,R) irrep for bra_space, ket_space, sp3r_space
  assert(bra_space.sigma() == ket_space.sigma());
  assert(sp3r_space.sigma() == ket_space.sigma());
  assert(sp3r_space.sigma() == bra_space.sigma());
  assert(bra_space.Nn_max() == ket_space.Nn_max());
  assert(sp3r_space.Nn_max() == ket_space.Nn_max());
  assert(sp3r_space.Nn_max() == bra_space.Nn_max());
  

  // Get dimensions of operator matrix
  int basis_size_bra = bra_space.GetBranchedDimension();
  int basis_size_ket = ket_space.GetBranchedDimension();

  // Get S, J for bra, ket, operator spaces
  HalfInt J_bra = bra_space.J();
  HalfInt J_ket = ket_space.J();
  
  HalfInt S_bra = bra_space.S();
  HalfInt S_ket = ket_space.S();
  
  // Initialize full operator matrix
  operator_matrix = basis::OperatorBlock<double>::Zero(basis_size_bra,basis_size_ket);

  // Coefficient caches
  u3::UCoefCache u_coef_cache;
  u3::WCoefCache w_coef_cache;

  // Loop over operators which contribute to observable
  for(const auto& op_pair : observable_list)
  {
    // Get operator and coefficient
    const auto& [operator_info,op_coef] = op_pair;

    const sp3r::OperatorSpaceSp3RSJ op_space = operator_info.OperatorSpaceSp3RSJ();
    HalfInt J_op = op_space.J0();
    HalfInt S_op = op_space.S0();

    // Construct hypersectors
    sp3r::Sp3RSJHypersectors hypersectors(bra_space,op_space);
    basis::OperatorHyperblocks<double> operator_hyperblocks;
    basis::SetHyperoperatorToZero(hypersectors,operator_hyperblocks);

    // Loop over sectors (bra/ket subspaces) and generate operator tiles for each sector
    for (std::size_t bra_index = 0; bra_index < bra_space.size(); bra_index++)
      for (std::size_t ket_index = 0; ket_index < ket_space.size(); ket_index++)
      {
        // Get subspaces and subspace offsets in matrix
        const sp3r::U3Subspace& bra_subspace = bra_space.GetSubspace(bra_index);
        const sp3r::U3Subspace& ket_subspace = ket_space.GetSubspace(ket_index);

        const int bra_subspace_offset = bra_space.GetBranchedSubspaceOffset(bra_index);
        const int ket_subspace_offset = ket_space.GetBranchedSubspaceOffset(ket_index);

        const u3::U3& omega_bra = bra_subspace.U3();
        const u3::U3& omega_ket = ket_subspace.U3();

        // Get dimension of tile in operator matrix
        // Tile dimension = upsilon_max*(total states in U3 irrep)
        const int tile_dimension_bra = bra_subspace.dimension()*bra_subspace.total_states();
        const int tile_dimension_ket = ket_subspace.dimension()*ket_subspace.total_states();

        // Construct hyperblocks corresponding to current bra-ket sector and hypersectors
        sp3r::ConstructSp3RHyperBlocks(
          sp3r_space,bra_space,ket_space,op_space,bra_subspace,ket_subspace,
          bra_index,ket_index,J_bra,J_ket,J_op,S_bra,S_ket,S_op,
          u_coef_cache,w_coef_cache,operator_info,
          hypersectors,operator_hyperblocks
        );
        
        // Construct operator tile using computed hypersectors and hyperblocks for bra/ket subspace
        // <bra|Op|ket> = \sum_{Op labels} <bra|Op^(tensor)|ket>
        // <bra|Op^(tensor)|ket> = \sum (LSJ 9j symbol)*(SU(3) Wcoef)*<bra||Op^(tensor)||ket>

        // Initialize tile for sector
        basis::OperatorBlock<double> tile;
        // Construct operator tile
        sp3r::GetSp3ROperatorTile(
          sp3r_space,bra_space,ket_space,op_space,bra_subspace,ket_subspace,
          bra_index,ket_index,J_bra,J_ket,J_op,S_bra,S_ket,S_op,
          u_coef_cache,w_coef_cache,operator_info,
          hypersectors,operator_hyperblocks,tile
        );

        // Place tile in matrix
        operator_matrix.block(bra_subspace_offset,ket_subspace_offset,tile_dimension_bra,tile_dimension_ket) += op_coef*tile;

      }
      // End looping all bra/ket sectors

  }
  // End looping all operators
}


}