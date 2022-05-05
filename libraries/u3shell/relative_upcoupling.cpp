/****************************************************************
  recurrence_upcoupling.cpp
                                  
  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT
****************************************************************/
#include "u3shell/relative_upcoupling.h"

#include <unordered_set>
#include "am/am.h"
#include "am/wigner_gsl.h"
#include "fmt/format.h"
#include "u3shell/operator_indexing_sectors.h"
#include "utilities/utilities.h"
// extern double zero_threshold;

namespace u3shell::relative
{

  basis::OperatorBlock<double> UpcoupleLST(
    const std::tuple<unsigned int,unsigned int, unsigned int, unsigned int>& labels_op,
    const std::tuple<unsigned int,unsigned int, unsigned int>& labels_bra,
    const std::tuple<unsigned int,unsigned int, unsigned int>& labels_ket,
    const basis::RelativeSectorsLSJT& sectors_lsjt,
    basis::OperatorBlocks<double>& blocks_lsjt
    )
  {
    const auto&[L0,S0,T0,J0] = labels_op;
    const auto&[Lbar,Sbar,Tbar] = labels_ket;
    const auto&[Lbarp,Sbarp,Tbarp] = labels_bra;

    const auto& bra_space_lsjt = sectors_lsjt.bra_space();
    const auto& ket_space_lsjt = sectors_lsjt.ket_space();

    const int Jbar_min = abs(static_cast<int>(Lbar) - static_cast<int>(Sbar));
    const int Jbarp_min = abs(static_cast<int>(Lbarp) - static_cast<int>(Sbarp));

    std::size_t bra_index = bra_space_lsjt.LookUpSubspaceIndex({Lbarp, Sbarp, Jbarp_min, Tbarp, Lbarp%2});
    std::size_t ket_index = ket_space_lsjt.LookUpSubspaceIndex({Lbar, Sbar, Jbar_min, Tbar, Lbar%2});

    basis::OperatorBlock<double> temp;
    if(bra_index==basis::kNone || ket_index == basis::kNone)
      // return zero sized block
      return temp;

    // Otherwise initialize block as zero matrix.
    std::size_t dimp = bra_space_lsjt.GetSubspace(bra_index).dimension();
    std::size_t dim = ket_space_lsjt.GetSubspace(ket_index).dimension();
    temp = basis::OperatorBlock<double>::Zero(dimp, dim);

    // Summing over Jbar and Jbarp
    for (int Jbar = Jbar_min; Jbar <= Lbar + Sbar; ++Jbar)
      for (int Jbarp = Jbarp_min; Jbarp <= Lbarp + Sbarp; ++Jbarp)
      {
        if (!am::AllowedTriangle(Jbar, J0, Jbarp))
          continue;

          size_t bra_index = bra_space_lsjt.LookUpSubspaceIndex({Lbarp, Sbarp, Jbarp, Tbarp, Lbarp%2});
          size_t ket_index = ket_space_lsjt.LookUpSubspaceIndex({Lbar, Sbar, Jbar, Tbar, Lbar%2});

          bool get_adoint_sector = bra_index>ket_index;

          size_t block_index;
          if(get_adoint_sector)
            block_index = sectors_lsjt.LookUpSectorIndex(ket_index,bra_index);
          else
            block_index = sectors_lsjt.LookUpSectorIndex(bra_index, ket_index);

          auto& block = blocks_lsjt[block_index];

          // If the sector is diagonal, then we need to fill in lower triangle of the sector
          // because only the upper triangle was saved to file.
          if (bra_index==ket_index)
          {
            int nmax = block.cols() - 1;
            for (int n = 0; n <= nmax; ++n)
              for (int np = 0; np < n; ++np)
              {
                block(n, np) = block(np, n);
              }
          }

          // accumulating blocks
          double so3_coef =
              am::Unitary9J(Lbar, Sbar, Jbar, L0, S0, J0, Lbarp, Sbarp, Jbarp)
              * (am::dim(Jbarp) * am::dim(S0) * am::dim(L0))
              / (am::dim(J0) * am::dim(Sbarp) * am::dim(Lbarp));

          if(get_adoint_sector)
          {

            assert(block.rows() == dim);
            assert(block.cols() == dimp);

            // It is assumed here that the operator is a spherical tensor,
            // i.e., (V^{J_0}_{M_0})^\dagger = (-1)^{M_0}V^{J_0}_{\tilde{M}_0}. Thus
            // \braket{J'||V^{J_0}||J}=(-1)^{J'-J} \hat{J}/\hat{J'}\braket{J||V^{J_0}||J'}
            // SU_J(2) and SU_T(2).
            double conjugation_factor
              = ParitySign(Jbar-Jbarp+Tbar-Tbarp)*(Hat(Jbar)*Hat(Tbar))/(Hat(Jbarp)*Hat(Tbarp));
            temp += so3_coef * conjugation_factor*block.transpose();

          }
          else
          {
            assert(block.rows() == dimp);
            assert(block.cols() == dim);

            temp += so3_coef * block;
          }
      }
    return temp;
  }


  std::vector<double> UpcoupleU3ST(
      const int Nbar_max,
      const uint8_t g0,
      const u3shell::relative::OperatorSectors& sectors_u3st,
      const std::array<basis::RelativeSectorsLSJT,3>& component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& component_blocks
    )
  {
    //Note: Currently only considering parity conserving operators...
    assert(g0==0);
    const auto& J0 = sectors_u3st.J0();
    const size_t num_elements = sectors_u3st.num_elements();
    std::vector<double> rme_array(num_elements,0.0);

    for (size_t index_u3st = 0; index_u3st < sectors_u3st.size(); ++index_u3st)
    {
      const auto& sector_u3st = sectors_u3st.GetSector(index_u3st);
      const auto sector_u3st_offset = sectors_u3st.GetSectorOffset(index_u3st);
      const u3shell::spatial::onecoord::OperatorL0Space& spatial_subspace =
          sector_u3st.bra_subspace();
      const u3shell::spin::twobody::OperatorSubspace& spin_subspace =
          sector_u3st.ket_subspace();
      const auto& S0 = spin_subspace.S0();
      const auto& exchange_symm_bar = spin_subspace.exchange_symm_bar();
      const auto& L0 = spatial_subspace.L0();
      const int parity_bar = (exchange_symm_bar + 1) % 2;

      for (size_t index_spin_state = 0; index_spin_state < spin_subspace.size();
           ++index_spin_state)
      {
        const auto& [S00, T0, Sbar, Sbarp, Tbar, Tbarp] =
            spin_subspace.GetState(index_spin_state).labels();

        assert(S00 == S0);
        const basis::RelativeSectorsLSJT& sectors_lsjt = component_sectors[T0];

        if(sectors_lsjt.size()==0)
          continue;

        basis::OperatorBlocks<double>& blocks = component_blocks[T0];

        for (int Lbar = parity_bar; Lbar <= Nbar_max; Lbar += 2)
          for (int Lbarp = parity_bar; Lbarp <= Nbar_max; Lbarp += 2)
          {
            if (!am::AllowedTriangle(Lbar, L0, Lbarp))
              continue;

            basis::OperatorBlock<double> block_lst;
            const std::tuple<unsigned int,unsigned int, unsigned int> bra_label{Lbarp,Sbarp,Tbarp};
            const std::tuple<unsigned int,unsigned int, unsigned int> ket_label{Lbar,Sbar,Tbar};

            // Upcouple to LST
            block_lst = UpcoupleLST({L0,S0,T0,J0},bra_label,ket_label,sectors_lsjt,blocks);
            if(block_lst.rows()==0)
              continue;

            // Upcoupling to U3ST
            // RMEs are stored in array:
            //  -> by L0x{S0,exchange_bar} sector
            //    -> by spin_state (T0,Sbar,Sbarp,Tbar,Tbarp)
            //      -> by x0, kappa0, Nbar (recall each sector has fixed N0 so Nbarp is fixed.)
            for (size_t index_x0 = 0; index_x0 < spatial_subspace.size(); ++index_x0)
            {
              const auto& x0_subspace = spatial_subspace.GetSubspace(index_x0);
              const auto& x0 = x0_subspace.x0();
              const auto N0 = x0_subspace.N0();
              for (int kappa0 = 1;
                   kappa0 <= spatial_subspace.GetSubspaceDegeneracy(index_x0);
                   ++kappa0)
              {
                const auto x0kappa0_offset =
                    spatial_subspace.GetSubspaceOffset(index_x0, kappa0);
                for (int index_spatial_state=0; index_spatial_state<x0_subspace.size(); ++index_spatial_state)
                {
                  const auto& spatial_state = x0_subspace.GetState(index_spatial_state);
                  const auto Nbar = spatial_state.Nbar();
                  if(Nbar<Lbar) continue;
                  const auto Nbarp = spatial_state.Nbarp();
                  if(Nbarp<Lbarp) continue;

                  int n = int((static_cast<int>(Nbar) - static_cast<int>(Lbar)) / 2);
                  int np = int((static_cast<int>(Nbarp) - static_cast<int>(Lbarp)) / 2);

                  size_t array_index =
                      sector_u3st_offset +
                      sector_u3st.element_offset(
                        index_spin_state,x0kappa0_offset,index_spatial_state
                      );

                  rme_array[array_index] +=
                      u3::W({Nbar,0u}, 1, Lbar, x0, kappa0, L0, {Nbarp,0u}, 1, Lbarp, 1)
                      * (u3::dim(x0) * am::dim(Lbarp))
                      / (1.*u3::dim({Nbarp,0u}) * am::dim(L0))
                      * ParitySign(n + np)
                      * block_lst(np, n);
                }
              }
            }
          }
      }
    }

    return rme_array;
  }


  std::vector<double> UpcoupleU3ST(
      const unsigned int Nbar_max,
      const u3shell::relative::OperatorSectors& sectors_u3st,
      const std::string& input_filename,
      const unsigned int input_Nbar_max,
      const unsigned int input_Jmax
    )
  {
    bool verbose = true;
    assert(utils::FileExists(input_filename, verbose));
    basis::RelativeSpaceLSJT relative_space_lsjt(input_Nbar_max, input_Jmax);
    std::array<basis::RelativeSectorsLSJT, 3> sectors_lsjt;
    std::array<basis::OperatorBlocks<double>, 3> blocks_lsjt;
    basis::RelativeOperatorParametersLSJT op_labels_lsjt;

    basis::ReadRelativeOperatorLSJT(
        input_filename,
        relative_space_lsjt,
        op_labels_lsjt,
        sectors_lsjt,
        blocks_lsjt,
        true
      );

    std::vector<double> rmes = UpcoupleU3ST(
        Nbar_max, op_labels_lsjt.g0, sectors_u3st, sectors_lsjt, blocks_lsjt
      );

    return UpcoupleU3ST(
        Nbar_max, op_labels_lsjt.g0, sectors_u3st, sectors_lsjt, blocks_lsjt
      );
  }


  std::pair<u3shell::relative::OperatorParameters,std::vector<double>> UpcoupleU3ST(
      const unsigned int Nbar_max,
      const std::string& input_filename,
      const unsigned int input_Nbar_max,
      const unsigned int input_Jmax
    )
  {
    bool verbose = true;
    assert(utils::FileExists(input_filename, verbose));
    basis::RelativeSpaceLSJT relative_space_lsjt(input_Nbar_max, input_Jmax);
    std::array<basis::RelativeSectorsLSJT, 3> sectors_lsjt;
    std::array<basis::OperatorBlocks<double>, 3> blocks_lsjt;
    basis::RelativeOperatorParametersLSJT op_labels_lsjt;

    basis::ReadRelativeOperatorLSJT(
        input_filename,
        relative_space_lsjt,
        op_labels_lsjt,
        sectors_lsjt,
        blocks_lsjt,
        true
      );

    std::set<uint8_t> Allowed_T0_values;
    for(uint8_t T0 = op_labels_lsjt.T0_min; T0<=op_labels_lsjt.T0_max; ++T0)
      Allowed_T0_values.insert(T0);

    u3shell::relative::OperatorParameters parameters(
        Nbar_max, op_labels_lsjt.J0, {}, {}, {}, Allowed_T0_values
      );

    u3shell::relative::OperatorSectors sectors_u3st(
        std::make_shared<const u3shell::spatial::onecoord::OperatorSpace>(parameters),
        std::make_shared<const u3shell::spin::twobody::OperatorSpace>(parameters),
        parameters.J0
      );

    std::vector<double> rmes = UpcoupleU3ST(
        Nbar_max, op_labels_lsjt.g0, sectors_u3st, sectors_lsjt, blocks_lsjt
      );

    return {parameters,rmes};
  }


  }  // namespace u3shell::relative
