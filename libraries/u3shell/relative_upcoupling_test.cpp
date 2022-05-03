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
#include "u3shell/operator_parameters.h"
#include "basis/lsjt_scheme.h"
#include "basis/lsjt_operator.h"
#include "u3shell/relative_upcoupling.h"
#include "sp3rlib/u3coef.h"
#include "utilities/utilities.h"
#include "u3shell/relative_operator.h"
#include "u3shell/upcoupling.h"
#include "boost/math/constants/constants.hpp"

namespace u3shell
{
  // Define test function for upcoupled interaction. 
  double JISP16RME(
    const TensorLabelsU3ST& tensor_labels,
    const StateLabelsNST& bra,
    const StateLabelsNST& ket
  )
  {
    const auto& [x0,S0,T0] = tensor_labels;
    const auto& [Nbar,Sbar,Tbar] = ket;
    const auto& [Nbarp,Sbarp,Tbarp] = bra;
    
  }
}


int main(int argc, char **argv)
{
  u3::U3CoefInit(39);

    // unsigned int J0=0;
    // Nbar_max=4
    // int Nmax=2, N1v=1;
    // int Jmax=Nmax+2*N1v+1;

    // std::unordered_set<u3::U3> Allowed_w0_values = {{0,{0u,0u}}};
    // std::set<unsigned int> Allowed_L0_values =  {0};
    // std::set<uint8_t> Allowed_S0_values = {0};
    // std::set<uint8_t> Allowed_T0_values = {0};
    //Nbar_max,J0,Allowed_w0_values,Allowed_L0_values,Allowed_S0_values,Allowed_T0_values
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    std::vector<std::tuple<u3shell::relative::OperatorParameters,
      u3shell::relative::RMEFunction,
      int,
      std::string>
    >
    operator_data = {
      {
        u3shell::relative::IdentityParameters(4),
        u3shell::relative::IdentityRME,
        4,
        "relative_identity_Nmax04.dat",
      },
      {
        {10,2,{},{2},{0},{0}},
        u3shell::relative::QuadrupoleRME,
        20,
        "relative_quadrupole_T0_Nmax20.dat"
      },

      {
        {10,2,{},{2},{0},{1}},
        u3shell::relative::QuadrupoleRME,
        20,
        "relative_quadrupole_T1_Nmax20.dat"
      },
      {
        // HamiltonianParameters(const unsigned int Nbar_max)
        // Lambda function for testing that uses old upcoupling
      }


    };

    for(const auto& [parameters,rme_function,Nmax_lsjt,filename] : operator_data)
    {
      std::string inputfile = fmt::format(
          "{}/spncci/data/relative_interactions/{}",
          utils::get_spncci_project_root_dir(),
          filename
        );

      fmt::print("{}\n", inputfile);
      bool verbose = true;
      assert(utils::FileExists(inputfile, verbose));

      auto spatial_ptr =
          std::make_shared<const u3shell::spatial::onecoord::OperatorSpace>(
              parameters
            );

      auto spin_ptr =
          std::make_shared<const u3shell::spin::twobody::OperatorSpace>(parameters);

      u3shell::relative::OperatorSectors sectors_u3st(
          spatial_ptr, spin_ptr, parameters.J0
        );
      // std::cout<<sectors_u3st.DebugStr()<<std::endl;

      ////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Getting lsjt sectors and rmes
      int Jmax = Nmax_lsjt + 1;
      basis::RelativeSpaceLSJT relative_space_lsjt(Nmax_lsjt, Jmax);
      std::array<basis::RelativeSectorsLSJT, 3> isospin_component_sectors_lsjt;
      std::array<basis::OperatorBlocks<double>, 3> isospin_component_matrices_lsjt;
      basis::RelativeOperatorParametersLSJT op_labels_lsjt;
      std::cout<<"Reading operator "<<std::endl;

      basis::ReadRelativeOperatorLSJT(
          inputfile,
          relative_space_lsjt,
          op_labels_lsjt,
          isospin_component_sectors_lsjt,
          isospin_component_matrices_lsjt,
          true
        );


      fmt::print("upcoupling\n");
      uint8_t g0 = 0;
      std::vector<double> rme_array = u3shell::relative::UpcoupleU3ST(
          parameters.Nbar_max,
          g0,
          sectors_u3st,
          isospin_component_sectors_lsjt,
          isospin_component_matrices_lsjt
      );


      fmt::print("Constructing relative operator\n");
      u3shell::relative::RelativeOperator relative_operator(parameters,rme_function);

      if(true)
      {
        const auto& spatial_space = sectors_u3st.bra_space();
        const auto& spin_space = sectors_u3st.ket_space();
        fmt::print("num sectors {}\n",sectors_u3st.size());
        fmt::print("num elements {}\n",sectors_u3st.num_elements());
        assert(sectors_u3st.num_elements()==rme_array.size());

        for(int sector_index = 0; sector_index<sectors_u3st.size(); ++sector_index)
        {
          const auto& sector = sectors_u3st.GetSector(sector_index);
          const auto& parity_space =
              spatial_space.GetSubspace(sector.parity_space_index());
          const auto parity_bar = parity_space.parity_bar();
          const auto& N0_space =
              parity_space.GetSubspace(sector.N0_space_index());
          const auto N0 = N0_space.N0();
          const auto& L0_space = sector.bra_subspace();
          const auto& S0_space = sector.ket_subspace();
          const auto L0 = L0_space.L0();
          const auto& [exchange_symm_bar,S0] = S0_space.labels();
          
          const auto sector_offset = sectors_u3st.GetSectorOffset(sector_index);
          for (int x0_index = 0; x0_index < L0_space.size(); ++x0_index)
          {
            const auto& x0_subspace = L0_space.GetSubspace(x0_index);
            const auto& x0 = x0_subspace.x0();
            const auto kappa0_max = L0_space.GetSubspaceDegeneracy(x0_index);
            
            for (unsigned int kappa0 = 1; kappa0 <= kappa0_max; ++kappa0)
            {
              const auto x0_kappa0_index =
                  L0_space.GetSubspaceOffset(x0_index, kappa0);
              for (int Nbar_index = 0; Nbar_index < x0_subspace.size();
                   ++Nbar_index)
              {
                const auto& state = x0_subspace.GetState(Nbar_index);
                const auto Nbar = state.Nbar();
                const auto Nbarp = state.Nbarp();
                
                for (std::size_t spin_index = 0; spin_index < S0_space.size();
                     ++spin_index)
                {
                  const auto spin_state = S0_space.GetState(spin_index);
                  const auto& [S00, T0, Sbar, Sbarp, Tbar, Tbarp] =
                      spin_state.labels();
                 const auto element_offset = sector.element_offset(
                      spin_index, x0_kappa0_index, Nbar_index
                    );
                  const auto element_offset2 =
                      sector.element_offset(spin_index, x0, kappa0, Nbar);
                  assert(element_offset == element_offset2);
                  std::set<u3::SU3> qx0_values = {{2,0},{1,1},{0,2}};
                  
                  const auto& operator_rmes = relative_operator.rmes();
                  int array_index = sector_offset + element_offset;
                  if(fabs(rme_array[array_index]-operator_rmes[array_index])>1e-10)
                  {
                    fmt::print("parity: {}  N0: {}  L0: {} S0: {}\n", parity_bar, N0, L0, S0);
                    fmt::print("spins: {} [{} {}] [{} {}]\n", T0, Sbar, Sbarp, Tbar, Tbarp);
                    fmt::print("x0: {}  kappa0_max: {}\n", x0, kappa0_max);
                    fmt::print("Nbar: {} {}\n", Nbarp, Nbar);

                    std::cout
                        << rme_array[sector_offset + element_offset] << std::endl;
                    std::cout
                        <<
                        operator_rmes[array_index]
                        << std::endl;
                    std::cout<<
                    sqrt(5./(16*boost::math::constants::pi<double>()))*u3shell::Qrel({Nbarp, Sbarp, Tbarp}, {Nbar, Sbar, Tbar})
                    <<std::endl;
                  }
                }
              }
            }
          }
        }
      }
    }
    }
