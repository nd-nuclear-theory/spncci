/****************************************************************
  seed_gen.cpp

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT
****************************************************************/
#include "u3ncsm/seed_gen.h"

#include "LookUpContainers/CWig9lmLookUpTable.h"
#include "SU3ME/CInteractionPN.h"
#include "SU3ME/InteractionPPNN.h"
#include "UNU3SU3/UNU3SU3Basics.h"
#include "mcutils/eigen.h"
#include "mcutils/parsing.h"
// #include "spncci/recurrence_seeds.h"
#include "u3ncsm/dimensions.h"
#include "u3ncsm/lgi_gen.h"
#include "u3ncsm/su3rme.h"
#include "utilities/utilities.h"


namespace spncci::seeds
{

basis::OperatorBlock<double> GenerateRecurrenceSeedBlock(
    const nuclide::NuclideType& nuclide,
    const HalfInt& Nsigma0,
    const int N1v,
    const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& unit_tensor_labels,
    const spncci::spin::RecurrenceSpace<lgi::LGI, u3shell::spin::twobody::OperatorLabelsST>&
        spin_recurrence_space,
    const spncci::spatial::RecurrenceSpace<u3shell::spatial::OneCoordType>& spatial_recurrence_space,
    const CBaseSU3Irreps& baseSU3Irreps,
    const std::string& operator_dir,
    const int recurrence_sp3r_space_index
  )
{
  std::cout<<"spatial "<<spatial_recurrence_space.size()<<std::endl;
  basis::OperatorBlock<double> recurrence_seed_block;

  MPI_Comm individual_comm = MPI_COMM_SELF;
  const auto& [Z, N] = nuclide;

  // Get relevant spin::RecurrenceLGISpace space and extract
  // sigma, sigma' and operator parity information
  ////////////////////////////////////////////////////////////////////////////////////
  // spatial::RecurrenceSpace() []
  // ->spatial::RecurrenceSp3RSpace() [sigma,sigma',parity_bar]
  //   ->spatial::RecurrenceNnsumSpace() [Nsum]
  //     ->spatial::RecurrenceU3Space() [omega,omega']->(upsilon x upsilon')
  //       ->spatial::RecurrenceOperatorSubspace() [x0] ->rho0_max
  //         ->spatial::RecurrenceOperatorState() [Nbar,Nbar']
  //
  // Get spatial operator subspace for given sigma_ket,sigma_bra
  // Note: for seeds we only care about the Nsum=0 subspace (with index=0) which
  // has a single spatial::RecurrenceU3Space with omega=sigma and omega' = sigma'
  // which also has index=0
  ///////////////////////////////////////////////////////////////////////////////////
  const auto& spatial_recurrence_sp3r_space
      = spatial_recurrence_space.GetSubspace(recurrence_sp3r_space_index);
  const auto& [sigma_ket, sigma_bra, parity_bar]
      = spatial_recurrence_sp3r_space.labels();
  const auto& spatial_recurrence_u3s_space
      = spatial_recurrence_sp3r_space.GetSubspace(0).GetSubspace(0);

  int Nex_ket(sigma_ket.N() - Nsigma0);
  int Nex_bra(sigma_bra.N() - Nsigma0);

  int exchange_symm_bar = (parity_bar + 1) % 2;
  int spin_recurrence_space_index = spin_recurrence_space.LookUpSubspaceIndex(
      {sigma_ket, sigma_bra, exchange_symm_bar}
    );
  const auto& spin_recurrence_lgi_space
      = spin_recurrence_space.GetSubspace(spin_recurrence_space_index);

  // // If there is no corresponding
  if (spin_recurrence_space_index == basis::kNone)
  {
    std::cout << "No corresponding spin space for " << sigma_ket.Str() << " "
              << sigma_bra.Str() << " " << parity_bar << std::endl;
    return recurrence_seed_block;
  }

  // Allocate seed block
  recurrence_seed_block = basis::OperatorBlock<double>::Zero(
      spatial_recurrence_u3s_space.dimension(),
      spin_recurrence_lgi_space.dimension()
    );

  // Get LGI expansions
  auto lgi_expansions =
      lgi::GetLGIExpansions(nuclide, spin_recurrence_lgi_space, Nsigma0);

  // Loop over all unit tensors
  // TODO: openMP parallelize over operator index.
  for (int operator_index = 0; operator_index < unit_tensor_labels.size();
       ++operator_index)
  {
    // Extract operator labels and apply spatial selection rules
    const auto& [operator_labels, relative_bra, relative_ket] =
        unit_tensor_labels[operator_index].Key();
    const auto& [dN0, x0, S0, T0, g0] = operator_labels;
    if (dN0 + sigma_ket.N() != sigma_bra.N())
      continue;

    int spatial_operator_index =
        spatial_recurrence_u3s_space.LookUpSubspaceIndex(x0);
    if (spatial_operator_index == basis::kNone)
      continue;
    const auto& spatial_recurrence_operator_subspace =
        spatial_recurrence_u3s_space.GetSubspace(spatial_operator_index);
    // std::cout<<"spatial_operator_index "<<spatial_operator_index<<std::endl;

    const auto& [Nbar, Sbar, Tbar] = relative_ket;
    const auto& [Nbarp, Sbarp, Tbarp] = relative_bra;
    int spatial_recurrence_operator_state_index =
        spatial_recurrence_operator_subspace.LookUpStateIndex({{Nbar}});
    if (spatial_recurrence_operator_state_index == basis::kNone)
      continue;
    // std::cout<<"spatial_recurrence_operator_state_index
    // "<<spatial_recurrence_operator_state_index<<std::endl; Look up spatial
    // and spin operator subspace indices

    int rhot_max = u3::OuterMultiplicity(sigma_ket.SU3(), x0, sigma_bra.SU3());
    // assert(rhot_max==spatial_recurrence_operator_subspace.degeneracy());
    // if(rhot_max<1) continue;

    // Definite labels for use in generate seeds
    auto spin_tensor_labels = u3shell::spin::twobody::OperatorLabelsST(
        int(S0), int(T0), int(Sbar), int(Sbarp), int(Tbar), int(Tbarp)
      );
    SU3xSU2::LABELS w0(1, x0.lambda(), x0.mu(), TwiceValue(S0));

    // Read in operator from file
    std::string operator_base_name =
        fmt::format("{}/relative_unit_{:06}", operator_dir, operator_index);
    std::ofstream interaction_log_file("/dev/null");
    bool log_is_on = false;
    bool generate_missing_rme = true;

    // The constructor for interactionPPNN must always be called before the
    // constructor for interactionPN otherwise a global variable doesn't get
    // correctly allocated...F#$%*$
    CInteractionPPNN interactionPPNN(
        baseSU3Irreps, log_is_on, interaction_log_file
      );
    CInteractionPN interactionPN(
        baseSU3Irreps, generate_missing_rme, log_is_on, interaction_log_file
      );

    interactionPPNN.LoadTwoBodyOperator(operator_base_name + ".PPNN");
    interactionPN.AddOperator(operator_base_name + ".PN");
    interactionPPNN.TransformTensorStrengthsIntoPP_NN_structure();
    ////////////////////////////////////////////////////////////////////////////////////
    // Iterating through spin space
    ////////////////////////////////////////////////////////////////////////////////////
    // spin::RecurrenceSpace() []
    //   spin::RecurrenceLGISpace() [sigma,sigma',exchange_symm_bar]
    //   ->spin::RecurrenceSpinSpace() [S,S']
    //     -> spin::RecurrenceSpinSubspace() [Sp,Sn,Sp',Sn']/[T,T']->(gamma,gamma')
    //       ->spin::RecurrenceOperatorState() [operator_index]
    ////////////////////////////////////////////////////////////////////////////////////
    for (int spin_space_index = 0;
         spin_space_index < spin_recurrence_lgi_space.size();
         ++spin_space_index)
    {
      const auto& spin_recurrence_spin_space =
          spin_recurrence_lgi_space.GetSubspace(spin_space_index);
      const auto& [S_ket, S_bra] = spin_recurrence_spin_space.labels();
      if (!am::AllowedTriangle(S_ket, S0, S_bra))
        continue;

      for (int spin_subspace_index = 0;
           spin_subspace_index < spin_recurrence_spin_space.size();
           ++spin_subspace_index)
      {
        const auto& spin_recurrence_subspace =
            spin_recurrence_spin_space.GetSubspace(spin_subspace_index);
        const auto& ket_upstream_labels =
            spin_recurrence_subspace.ket_upstream_labels();
        const auto& bra_upstream_labels =
            spin_recurrence_subspace.bra_upstream_labels();
        if (!lgi::UnitTensorAllowed(
                ket_upstream_labels, spin_tensor_labels, bra_upstream_labels
              ))
          continue;

        const auto [Sp_ket, Sn_ket] = ket_upstream_labels;
        const auto [Sp_bra, Sn_bra] = bra_upstream_labels;

        // Look up spin operator index
        int spin_operator_index =
            spin_recurrence_subspace.LookUpStateIndex(spin_tensor_labels);
        // std::cout<<"spin_operator_index "<<spin_operator_index<<std::endl;
        int jdiag = 0;
        int ndiag = 1;

        // Define bra and ket basis
        proton_neutron::ModelSpaceMapType model_space_map_ket;
        model_space_map_ket[Nex_ket][{TwiceValue(Sp_ket), TwiceValue(Sn_ket)}]
                           [TwiceValue(S_ket)]
                               .emplace_back(
                                   1, sigma_ket.SU3().lambda(), sigma_ket.SU3().mu()
                                 );
        proton_neutron::ModelSpace ket_ncsmModelSpace(Z, N, model_space_map_ket);
        lsu3::CncsmSU3xSU2Basis ket_ncsm_basis;
        ket_ncsm_basis.ConstructBasis(
            ket_ncsmModelSpace, jdiag, ndiag, individual_comm
          );

        // Construct bra basis from model space
        proton_neutron::ModelSpaceMapType model_space_map_bra;
        model_space_map_bra[Nex_bra][{TwiceValue(Sp_bra), TwiceValue(Sn_bra)}]
                           [TwiceValue(S_bra)]
                               .emplace_back(
                                   1, sigma_bra.SU3().lambda(), sigma_bra.SU3().mu()
                                 );
        proton_neutron::ModelSpace bra_ncsmModelSpace(Z, N, model_space_map_bra);
        lsu3::CncsmSU3xSU2Basis bra_ncsm_basis;
        bra_ncsm_basis.ConstructBasis(
            bra_ncsmModelSpace, jdiag, ndiag, individual_comm
          );

        // std::cout<<"Calculate unit tensor RMEs in lsu3shell basis"<<std::endl;
        basis::OperatorBlocks<double> su3_unit_tensor_rmes = lsu3shell::CalculateRME(
            interactionPPNN, interactionPN, bra_ncsm_basis, ket_ncsm_basis, dN0, w0, rhot_max
          );

        // // Read in lgi expansion from file
        // std::string lgi_expansion_filename_bra
        // =lgi::lgi_expansion_filename(Z,N,{{sigma_bra,Sp_bra,Sn_bra,S_bra},Nex_bra});
        // std::string lgi_expansion_filename_ket
        // =lgi::lgi_expansion_filename(Z,N,{{sigma_ket,Sp_ket,Sn_ket,S_ket},Nex_ket});
        // basis::OperatorBlock<double> lgi_expansion_bra =
        // ReadOperatorBlockBinary(lgi_expansion_filename_bra); basis::OperatorBlock<double>
        // lgi_expansion_ket = ReadOperatorBlockBinary(lgi_expansion_filename_ket);

        // testing
        basis::OperatorBlock<double> lgi_expansion_bra =
            lgi_expansions[{sigma_bra, Sp_bra, Sn_bra, S_bra}];
        basis::OperatorBlock<double> lgi_expansion_ket =
            lgi_expansions[{sigma_ket, Sp_ket, Sn_ket, S_ket}];
        // assert(mcutils::IsZero(lgi_expansion_ket_test-lgi_expansion_ket,1e-8));
        // assert(mcutils::IsZero(lgi_expansion_bra_test-lgi_expansion_bra,1e-8));

        // Loop over outer multiplicity and transform each block from lsu3shell
        // basis to spncci basis. Save seeds into recurrence seed matrix
        for (int irhot = 0; irhot < rhot_max; ++irhot)
        {
          int rho0 = irhot + 1;
          // Row index in recurrence seed matrix
          int target_index_row =
              spatial_recurrence_u3s_space.GetSubspaceOffset(
                  spatial_operator_index, rho0
                )
              + spatial_recurrence_operator_state_index;


          basis::OperatorBlock<double> block =
              lgi_expansion_bra.transpose() * su3_unit_tensor_rmes[irhot]
              * lgi_expansion_ket;

          // // TEMP: write block to file for comparison with old seeds.
          // std::string temp_filename
          //   =
          //   fmt::format("temp/seeds_Z{:02d}_N{:02d}_Nex{:02d}_lm{:02d}_mu{:02d}_2Sp{:02d}_2Sn{:02d}_2S{:02d}_Nex{:02d}_lm{:02d}_mu{:02d}_2Sp{:02d}_2Sn{:02d}_2S{:02d}_rho{}_i{}.dat",
          //         Z,N,Nex_bra,sigma_bra.SU3().lambda(),sigma_bra.SU3().mu(),TwiceValue(Sp_bra),TwiceValue(Sn_bra),TwiceValue(S_bra),
          //         Nex_ket,sigma_ket.SU3().lambda(),sigma_ket.SU3().mu(),TwiceValue(Sp_ket),TwiceValue(Sn_ket),TwiceValue(S_ket),
          //         rho0,operator_index
          //       );

          // utils::WriteOperatorBlockBinary(block,temp_filename);

          for (int i = 0; i < block.rows(); i++)
            for (int j = 0; j < block.cols(); j++)
            {
              int gamma_bra = i + 1;
              int gamma_ket = j + 1;

              // Column index in recurrence seed matrix
              int target_index_col =
                  spin_recurrence_lgi_space.GetSubspaceOffset(spin_space_index)
                  + spin_recurrence_spin_space.GetSubspaceOffset(
                      spin_subspace_index, gamma_ket, gamma_bra
                    )
                  + spin_operator_index;

              recurrence_seed_block(target_index_row, target_index_col) =
                  block(i, j);
            }
        }
      }
    }
  }
  return recurrence_seed_block;
}

}  // namespace spncci::seeds
