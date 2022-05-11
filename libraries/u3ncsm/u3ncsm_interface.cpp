/****************************************************************
  seed_gen.cpp

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT
****************************************************************/
#include "u3ncsm/u3ncsm_interface.h"

#include "mcutils/eigen.h"
#include "mcutils/parsing.h"
#include "utilities/utilities.h"


namespace spncci::seeds{

    std::unordered_map<u3shell::U3SPN,basis::OperatorBlock<double>>
    GetLGIExpansions(
      const nuclide::NuclideType& nuclide,
      const spncci::spin::RecurrenceLGISpace<lgi::LGI,u3shell::spin::twobody::OperatorLabelsST>& lgi_space,
      const HalfInt& Nsigma0
      )
    {
      const auto&[Z,N] = nuclide;
      std::unordered_map<u3shell::U3SPN,basis::OperatorBlock<double>> lgi_expansions;
      std::unordered_set<lgi::LGI> temp_list;
      const auto& [sigma_ket,sigma_bra,exchange_symm_bar] = lgi_space.labels();
      int Nex_bra = int(sigma_bra.N()-Nsigma0);
      int Nex_ket = int(sigma_ket.N()-Nsigma0);
      for(const auto& spin_space : lgi_space)
        {
          const auto& [S_ket, S_bra] = spin_space.labels();
          for(const auto& spin_subspace : spin_space)
            {
              const auto& ket_upstream_labels = spin_subspace.ket_upstream_labels();
              const auto& bra_upstream_labels = spin_subspace.bra_upstream_labels();
              const auto [Sp_ket, Sn_ket] = ket_upstream_labels;
              const auto [Sp_bra, Sn_bra] = bra_upstream_labels;

              temp_list.insert({{sigma_bra,Sp_bra,Sn_bra,S_bra},Nex_bra});
              temp_list.insert({{sigma_ket,Sp_ket,Sn_ket,S_ket},Nex_ket});
            }
        }

      for(const auto& lgi : temp_list)
        {
          std::string filename = lgi::lgi_expansion_filename(Z,N,lgi);

          // TODO: Make construction optional if lgi expansion doesn't exist
          bool exit_on_nonexist=true, warn_on_overwrite=false;
          mcutils::FileExistCheck(filename,exit_on_nonexist,warn_on_overwrite);
          lgi_expansions[lgi.u3spn()] = utils::ReadOperatorBlockBinary(filename);
        }

      return lgi_expansions;
    }

}// namespaces
