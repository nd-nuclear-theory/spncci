/****************************************************************
  recurrence_upcoupling_test.cpp

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  4/2/22 (aem): Created.
****************************************************************/
// #include <fstream>
#include <cassert>
#include "fmt/format.h"
#include "u3shell/operator_indexing_spin.h"
#include "u3shell/operator_indexing_spatial.h"
#include "u3shell/operator_indexing_sectors.h"
#include "basis/lsjt_scheme.h"
#include "basis/lsjt_operator.h"
#include "u3shell/relative_upcoupling.h"
#include "sp3rlib/u3coef.h"
#include "utilities/utilities.h"


int main(int argc, char **argv)
{
  u3::U3CoefInit(100);

    unsigned int J0=0;
    // Nbar_max=4
    int Nmax=2, N1v=1;
    int Jmax=Nmax+2*N1v+1;

    std::unordered_set<u3::U3> Allowed_w0_values = {{0,{0u,0u}}};
    std::set<unsigned int> Allowed_L0_values =  {0};
    std::set<uint8_t> Allowed_S0_values = {0};
    std::set<uint8_t> Allowed_T0_values = {0};

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Setting up u3st sectors
    u3shell::relative::OperatorParameters
      identity_parameters(
        N1v,Nmax,J0,
        Allowed_w0_values,
        Allowed_L0_values,
        Allowed_S0_values,
        Allowed_T0_values
        );

    auto spatial_ptr = std::make_shared<const u3shell::spatial::onecoord::OperatorSpace>(identity_parameters);
    auto spin_ptr = std::make_shared<const u3shell::spin::twobody::OperatorSpace>(identity_parameters);
    u3shell::relative::OperatorSectors sectors_u3st(spatial_ptr,spin_ptr,identity_parameters.J0);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Getting lsjt sectors and rmes
    std::string inputfile
    = fmt::format("{}/spncci/data/relative_interactions/relative_identity_Nmax04.dat",
      utils::get_spncci_project_root_dir()
      );

    fmt::print("{}\n",inputfile);
    bool verbose=true;
    assert(utils::FileExists(inputfile, verbose));

    basis::RelativeSpaceLSJT relative_space_lsjt(Nmax+2*N1v, Jmax);

    std::array<basis::RelativeSectorsLSJT,3> isospin_component_sectors_lsjt;
    std::array<basis::OperatorBlocks<double>,3> isospin_component_matrices_lsjt;
    basis::RelativeOperatorParametersLSJT op_labels_lsjt;

    basis::ReadRelativeOperatorLSJT(
      inputfile,relative_space_lsjt,op_labels_lsjt,
      isospin_component_sectors_lsjt, isospin_component_matrices_lsjt, true
      );

    fmt::print("upcoupling\n");
    uint8_t g0=0;
    std::vector<double> rme_array
      = u3shell::relative::UpcoupleU3ST(
        identity_parameters.Nbar_max,
        g0,
        sectors_u3st,
        isospin_component_sectors_lsjt,
        isospin_component_matrices_lsjt
      );


    const auto& spatial_space = sectors_u3st.bra_space();
    const auto& spin_space = sectors_u3st.ket_space();
    fmt::print("num sectors {}\n",sectors_u3st.size());
    fmt::print("num elements {}\n",sectors_u3st.num_elements());
    assert(sectors_u3st.num_elements()==rme_array.size());

    for(int sector_index = 0; sector_index<sectors_u3st.size(); ++sector_index)
      {
        const auto& sector = sectors_u3st.GetSector(sector_index);
        const auto& parity_space = spatial_space.GetSubspace(sector.parity_space_index());
        const auto parity_bar = parity_space.parity_bar();
        const auto& N0_space = parity_space.GetSubspace(sector.N0_space_index());
        const auto N0 = N0_space.N0();
        const auto& L0_space = sector.bra_subspace();
        const auto& S0_space = sector.ket_subspace();
        const auto L0 = L0_space.L0();
        const auto& [S0, exchange_symm_bar] = S0_space.labels();
        fmt::print("parity: {}  N0: {}  L0: {} S0: {}\n",parity_bar,N0,L0,S0);
        const auto sector_offset = sectors_u3st.GetSectorOffset(sector_index);
        for(int x0_index=0; x0_index<L0_space.size(); ++x0_index)
        {
          const auto& x0_subspace = L0_space.GetSubspace(x0_index);
          const auto& x0 = x0_subspace.x0();
          const auto kappa0_max = L0_space.GetSubspaceDegeneracy(x0_index);
          fmt::print("x0: {}  kappa0_max: {}\n",x0,kappa0_max);
          for(unsigned int kappa0=1; kappa0<=kappa0_max; ++kappa0)
          {
            const auto x0_kappa0_index = L0_space.GetSubspaceOffset(x0_index,kappa0);
            for (int Nbar_index = 0; Nbar_index < x0_subspace.size(); ++Nbar_index)
            {
              const auto& state = x0_subspace.GetState(Nbar_index);
              const auto Nbar = state.Nbar();
              const auto Nbarp = state.Nbarp();

                for(std::size_t spin_index=0; spin_index<S0_space.size(); ++spin_index)
                {
                  fmt::print("{} {}\n",Nbarp,Nbar);
                const auto element_offset
                  = sector.element_offset(spin_index,x0_kappa0_index,Nbar_index);
                const auto element_offset2 = sector.element_offset(spin_index,x0,kappa0,Nbar);
                assert(element_offset==element_offset2);
                std::cout<<rme_array[sector_offset+element_offset]<<std::endl;
                }


            }
          }

        }
      }


}
