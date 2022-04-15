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

// extern double zero_threshold;

namespace u3shell::relative
{
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
        // fmt::print("{} {} [{} {}, {} {}]\n",S0,T0,Sbarp,Tbarp,Sbar,Tbar);
        assert(S00 == S0);
        const basis::RelativeSectorsLSJT& sectors_lsjt = component_sectors[T0];
        basis::OperatorBlocks<double>& blocks = component_blocks[T0];

        const auto& bra_space_lsjt = sectors_lsjt.bra_space();
        const auto& ket_space_lsjt = sectors_lsjt.ket_space();

        for (int Lbar = parity_bar; Lbar <= Nbar_max; Lbar += 2)
          for (int Lbarp = parity_bar; Lbarp <= Nbar_max; Lbarp += 2)
          {
            if (!am::AllowedTriangle(Lbar, L0, Lbarp))
              continue;

            // from lsjt_scheme: indexing of states within sector
            // for (int N = L; N <= Nmax; N +=2)
            //   PushStateLabels(StateLabelsType(N));
            auto dimp = int((Nbar_max - Lbarp - parity_bar) / 2) + 1;
            auto dim = int((Nbar_max - Lbar - parity_bar) / 2) + 1;
            basis::OperatorBlock<double> temp =
                basis::OperatorBlock<double>::Zero(dimp, dim);

            // Summing over Jbar and Jbarp
            for (int Jbar = abs(Lbar - static_cast<int>(Sbar));
                 Jbar <= Lbar + Sbar;
                 ++Jbar)
              for (int Jbarp = abs(Lbarp - static_cast<int>(Sbarp));
                   Jbarp <= Lbarp + Sbarp;
                   ++Jbarp)
              {
                if (!am::AllowedTriangle(Jbar, J0, Jbarp))
                  continue;

                size_t bra_index = bra_space_lsjt.LookUpSubspaceIndex(
                    {Lbarp, Sbarp, Jbarp, Tbarp, parity_bar}
                  );
                size_t ket_index = ket_space_lsjt.LookUpSubspaceIndex(
                    {Lbar, Sbar, Jbar, Tbar, parity_bar}
                  );
                size_t block_index =
                    sectors_lsjt.LookUpSectorIndex(bra_index, ket_index);
                auto& block = blocks[block_index];

                assert(block.rows() == dimp);
                assert(block.cols() == dim);

                // If the sector is diagonal, then fill in lower triangle of the sector
                if (sectors_lsjt.GetSector(block_index).IsDiagonal())
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

                temp += so3_coef * block;
              }

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
                for (unsigned int Nbar = Lbar; Nbar <= Nbar_max; Nbar += 2)
                {
                  auto index_spatial_state = x0_subspace.LookUpStateIndex({Nbar});
                  const auto Nbarp = Nbar + N0;

                  int n = int((static_cast<int>(Nbar) - Lbar) / 2);
                  int np = int((static_cast<int>(Nbarp) - Lbarp) / 2);

                  size_t array_index =
                      sector_u3st_offset +
                      sector_u3st.element_offset(
                        index_spin_state,x0kappa0_offset,index_spatial_state
                      );

                  // fmt::print("index: {}  ({} {} ; {} {} {}|| {} {})  {}x{}\n",
                  //   array_index,Nbar,Lbar,x0,kappa0,L0,Nbarp,Lbarp,np,n
                  //   );

                  rme_array[array_index] +=
                      u3::W({Nbar,0u}, 1, Lbar, x0, kappa0, L0, {Nbarp,0u}, 1, Lbarp, 1)
                      * (u3::dim(x0) * am::dim(Lbarp))
                      / (1.*u3::dim({Nbar,0u}) * am::dim(L0)) * ParitySign(n + np)
                      * temp(np, n);
                }
              }
            }
          }
      }
    }

    return rme_array;
  }



  }  // namespace u3shell::relative
