/****************************************************************
  lgi_gen.cpp

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT
****************************************************************/
#include "u3ncsm/lgi_gen.h"

#include <omp.h>

#include <fstream>
#include <iostream>

#include "LSU3/ncsmSU3xSU2Basis.h"
#include "LookUpContainers/CWig9lmLookUpTable.h"
#include "SU3ME/CInteractionPN.h"
#include "SU3ME/InteractionPPNN.h"
#include "UNU3SU3/UNU3SU3Basics.h"
#include "am/am.h"
#include "lsu3shell/lsu3shell_basis.h"
#include "u3ncsm/dimensions.h"
#include "u3ncsm/su3rme.h"
#include "utilities/null_solver.h"

namespace lgi
{
// proton_neutron::ModelSpaceMapType generate_bra_model_space_Brel(
//     const lgi::LGI& lgi, const lgi::MultiplicityTaggedLGIVector& lgi_vector
//   )
// {
//   const auto& [Nex, sigma, Sp, Sn, S] = lgi.Key();
//   // bra model space consists subspaces connected by B[-2(0,2)]
//   proton_neutron::ModelSpaceMapType model_space_map_bra;
//   for (const auto& [lgi_p, dim_p] : lgi_vector)
//   {
//     const auto& [Nex_p, sigma_p, Sp_p, Sn_p, S_p] = lgi_p.Key();

//     // Apply selection rule
//     if (Nex - 2 != Nex_p)
//     {
//       continue;
//     }
//     if (Sp != Sp_p)
//     {
//       continue;
//     }
//     if (Sn != Sn_p)
//     {
//       continue;
//     }
//     if (S != S_p)
//     {
//       continue;
//     }
//     if (u3::OuterMultiplicity(sigma.SU3(), u3::SU3(0, 2), sigma_p.SU3()) == 0)
//     {
//       continue;
//     }

//     // If passed all selection rules, add to bra spaces
//     SU3::LABELS w(1, sigma_p.SU3().lambda(), sigma_p.SU3().mu());
//     model_space_map_bra[Nex_p][{TwiceValue(Sp_p), TwiceValue(Sn_p)}][TwiceValue(S_p)]
//         .push_back(w);
//   }
//   return model_space_map_bra;
// }


proton_neutron::ModelSpaceMapType
GenerateBraModelSpaceSp3RLoweringOperator(
    const HalfInt& Nsigma0,
    const proton_neutron::ModelSpaceMapType& model_space_map_ket,
    const std::map<u3shell::U3SPN, unsigned int>& u3ncsm_basis_map
  )
{
  proton_neutron::ModelSpaceMapType model_space_map_bra;

  for (const auto& [Nex, SPN_map] : model_space_map_ket)
  {
    if (int(Nex - 2) < 0)
      continue;

    for (const auto& [SSPN, S_map] : SPN_map)
    {
      const auto& [SSp, SSn] = SSPN;
      for (const auto& [SS, su3_vector] : S_map)
      {
        std::set<u3::SU3> bra_su3_irreps;
        for (const auto& x : su3_vector)
        {
          const auto product_vector =
              u3::KroneckerProduct({x.lm, x.mu}, {0u, 2u});

          for (const auto& [xp, mult] : product_vector)
          {
            u3shell::U3SPN bra_irrep(
                {Nsigma0 + Nex - 2, xp}, {SSp, 2}, {SSn, 2}, {SS, 2}
              );
            if (u3ncsm_basis_map.count(
                    {{Nsigma0 + Nex - 2, xp}, {SSp, 2}, {SSn, 2}, {SS, 2}}
                  ))
            {
              bra_su3_irreps.insert(xp);
            }
          }
        }

        for (const auto& xp : bra_su3_irreps)
        {
          // SU3::LABELS w(1,xp.lambda(),xp.mu());
          model_space_map_bra[Nex-2][SSPN][SS].push_back({1, xp.lambda(), xp.mu()});
        }
      }
    }
  }

  return model_space_map_bra;
}


basis::OperatorBlock<double> NcmMatrix(
    const int A,
    const double Ntot,
    const CBaseSU3Irreps& baseSU3Irreps,
    const lsu3::CncsmSU3xSU2Basis& bra_ncsm_basis,
    const lsu3::CncsmSU3xSU2Basis& ket_ncsm_basis,
    const std::string& nrel_operator_base_name
  )
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Calculate Ncm matrix
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
{
  //////////////////////////////////////////////////////////////////////////////////////////
  // Load Nrel operator used to calculate Ncm
  //////////////////////////////////////////////////////////////////////////////////////////
  // std::string nrel_operator_base_name = fmt::format("{}/Nrel",operator_dir);

  std::ofstream interaction_log_file("/dev/null");
  bool log_is_on = false;
  bool generate_missing_rme = true;

  CInteractionPPNN interactionPPNN_Nrel(
      baseSU3Irreps, log_is_on, interaction_log_file
    );
  CInteractionPN interactionPN_Nrel(
      baseSU3Irreps, generate_missing_rme, log_is_on, interaction_log_file
    );
  interactionPPNN_Nrel.LoadTwoBodyOperator(nrel_operator_base_name + ".PPNN");
  interactionPN_Nrel.AddOperator(nrel_operator_base_name + ".PN");
  interactionPPNN_Nrel.TransformTensorStrengthsIntoPP_NN_structure();


  // w0(rho0,lambda0,mu0,twice_S0)
  SU3xSU2::LABELS w0(1, 0, 0, 0);
  int dN0 = 0;
  int rhot_max = 1;

  // There is only one block in the operator_blocks returned by CalculateRME because rhot_max=1;
  basis::OperatorBlock<double> Ncm =
      -2. / A
      * lsu3shell::CalculateRME(
          interactionPPNN_Nrel,
          interactionPN_Nrel,
          bra_ncsm_basis,
          ket_ncsm_basis,
          dN0,
          w0,
          rhot_max
        )[0];

  int dim = Ncm.rows();
  Ncm += Ntot * Eigen::MatrixXd::Identity(dim, dim);

  return Ncm;
}


basis::OperatorBlock<double> BMatrix(
    const CBaseSU3Irreps& baseSU3Irreps,
    const lsu3::CncsmSU3xSU2Basis& bra_ncsm_basis,
    const lsu3::CncsmSU3xSU2Basis& ket_ncsm_basis,
    const std::string& brel_operator_base_name
  )
//////////////////////////////////////////////////////////////////////////////////////////
// Calculate Brel matrix
/////////////////////////////////////////////////////////////////////////////////////
{
  std::ofstream interaction_log_file("/dev/null");
  bool log_is_on = false;
  bool generate_missing_rme = true;

  CInteractionPPNN interactionPPNN_Brel(
      baseSU3Irreps, log_is_on, interaction_log_file
    );
  CInteractionPN interactionPN_Brel(
      baseSU3Irreps, generate_missing_rme, log_is_on, interaction_log_file
    );

  interactionPPNN_Brel.LoadTwoBodyOperator(brel_operator_base_name + ".PPNN");
  interactionPN_Brel.AddOperator(brel_operator_base_name + ".PN");
  interactionPPNN_Brel.TransformTensorStrengthsIntoPP_NN_structure();

  int dN0 = -2;
  // w0(rho0,lambda0,mu0,twice_S0)
  SU3xSU2::LABELS w0(1, 0, 2, 0);
  int rhot_max = 1;

  return lsu3shell::CalculateRME(
      interactionPPNN_Brel,
      interactionPN_Brel,
      bra_ncsm_basis,
      ket_ncsm_basis,
      dN0,
      w0,
      rhot_max
    )[0];
}


basis::OperatorBlock<double> GenerateLGIExpansion(
    const nuclide::NuclideType& nuclide,
    const HalfInt& Nsigma0,
    const unsigned int N0,
    const CBaseSU3Irreps& baseSU3Irreps,
    const std::map<u3shell::U3SPN, unsigned int>& u3ncsm_basis_map,
    const MultiplicityTagged<lgi::LGI>& tagged_lgi,
    const std::string& operator_dir,
    const MPI_Comm& individual_comm,
    double zero_threshold
  )
{
  // Notes: probably want to pull out operator load into main program and pass
  // in a argument If so, need to understand how Nop selection rules applied
  // internally. Nrel and Brel would have to be combined into one operator
  // because of funky pointer modification barried deep in lsu3shell code but
  // will have two different Nop which may be problematic. Aslo want to switch
  // to using B+Ncm rather than Bintr+Ncm.
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Calculate rmes of B and Ncm for each U3SpSnS subspace corresponding to an
  // LGI Then solve for simultaneous null space of B and Ncm
  // Extract label information
  const auto& [lgi, gamma_max] = tagged_lgi;
  // const auto& [lgi, gamma_max] = lgi_vector[lgi_index];
  const auto& [Nex, sigma, Sp, Sn, S] = lgi.Key();

  // If Nex==0, then all of the irreps in the subspace are lgi and the basis
  // for the null spaces is given by the identity matrix
  if (Nex == 0)
  {
    return Eigen::MatrixXd::Identity(gamma_max, gamma_max);
  }

  // With communicator split, ndiag is always 1
  int jdiag = 0;
  int ndiag = 1;
  const auto& [Z, N] = nuclide;

  // Model space consists of irreps with quantum numbers Nex (lambda,mu) Sp Sn S
  // map used for defining model_space_map organized into {N : {SpSn : {S : <(lambda,mu)>}}}
  proton_neutron::ModelSpaceMapType model_space_map_ket;
  model_space_map_ket[Nex][{TwiceValue(Sp), TwiceValue(Sn)}][TwiceValue(S)]
      .emplace_back(1, sigma.SU3().lambda(), sigma.SU3().mu());
  proton_neutron::ModelSpace ket_ncsmModelSpace(Z, N, model_space_map_ket);

  // Construct ket basis from model space
  // num in communicator must be equal to ndiag
  lsu3::CncsmSU3xSU2Basis ket_ncsm_basis;
  ket_ncsm_basis.ConstructBasis(ket_ncsmModelSpace, jdiag, ndiag, individual_comm);
  int dim_ket = lsu3shell::get_num_U3PNSPN_irreps(ket_ncsm_basis);

  // Nrel has quantum numbers 0(0,0)0 0, so bra and ket model space are the same
  // But we'll need a different bra model space for Brel
  // proton_neutron::ModelSpaceMapType model_space_map_bra =
  //     generate_bra_model_space_Brel(lgi, lgi_vector);

  // HalfInt Nsigma0 = nuclide::Nsigma0ForNuclide(nuclide);
  proton_neutron::ModelSpaceMapType model_space_map_bra
    = GenerateBraModelSpaceSp3RLoweringOperator(Nsigma0,model_space_map_ket,u3ncsm_basis_map);


  // If model_space_map_bra is not empty, then construct the basis for calculating B
  bool compute_B_matrix = false;
  lsu3::CncsmSU3xSU2Basis bra_ncsm_basis;
  int dim_B_bra = 0;
  if (model_space_map_bra.size() != 0)
  {
    proton_neutron::ModelSpace bra_ncsmModelSpace(Z, N, model_space_map_bra);
    bra_ncsm_basis.ConstructBasis(
        bra_ncsmModelSpace, jdiag, ndiag, individual_comm
      );
    dim_B_bra = lsu3shell::get_num_U3PNSPN_irreps(bra_ncsm_basis);
    compute_B_matrix = true;
  }

  int dim_bra = dim_ket + dim_B_bra;

  // Initialize matrix
  basis::OperatorBlock<double> BNcm_matrix =
      Eigen::MatrixXd::Zero(dim_bra, dim_ket);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Calculate Ncm matrix
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::string nrel_operator_base_name = fmt::format("{}/Nrel", operator_dir);
  BNcm_matrix.block(0, 0, dim_ket, dim_ket) = lgi::NcmMatrix(
      Z + N, double(N0 + Nex), baseSU3Irreps, ket_ncsm_basis, ket_ncsm_basis, nrel_operator_base_name
    );
  //////////////////////////////////////////////////////////////////////////////////////////
  // Calculate Brel matrix
  /////////////////////////////////////////////////////////////////////////////////////
  if (compute_B_matrix)
  {
    std::string brel_operator_base_name = fmt::format("{}/Brel", operator_dir);
    BNcm_matrix.block(dim_ket, 0, dim_B_bra, dim_ket) = BMatrix(
        baseSU3Irreps, bra_ncsm_basis, ket_ncsm_basis, brel_operator_base_name
      );
  }

  //////////////////////////////////////////////////////////////////////////////////////////
  // Get null vectors
  /////////////////////////////////////////////////////////////////////////////////////
  basis::OperatorBlock<double> lgi_expansion =
      utils::FindNullSpaceSVD(BNcm_matrix, gamma_max, zero_threshold);

  // TODO: Make optional
  if (true)
  {
    assert(mcutils::IsZero(BNcm_matrix * lgi_expansion, zero_threshold));
    assert(mcutils::IsZero(
        BNcm_matrix.block(dim_ket, 0, dim_B_bra, dim_ket) * lgi_expansion
      ));
    assert(mcutils::IsZero(
        BNcm_matrix.block(0, 0, dim_ket, dim_ket) * lgi_expansion, zero_threshold
      ));
    assert(not mcutils::IsZero(lgi_expansion, zero_threshold));
    assert(mcutils::IsZero(
        lgi_expansion.transpose() * lgi_expansion
            - Eigen::MatrixXd::Identity(gamma_max, gamma_max),
        zero_threshold
      ));
  }


  return lgi_expansion;
}
}  // namespace lgi
