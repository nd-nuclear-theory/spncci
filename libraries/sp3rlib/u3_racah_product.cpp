/****************************************************************
  u3_racah_product.cpp

  Colin V. Coane
  
  SPDX-License-Identifier: MIT 

  06/03/24 (cvc): Created. Split from u3_racah_product.cpp.
  06/03/24 (cvc): Refactored to .h/.cpp file format.
  07/02/24 (cvc): Functions for specific tensors moved to u3_tensor.h.

****************************************************************/

#include "sp3rlib/u3_racah_product.h"

namespace sp3r_operator
{


// Couple Sp(3,R)>U(3) operators [T^wT x S^wS]^w0
//
// Function inputs currently must have the form below:
// Tfunc(wT, rhoT, base_space, bra_space, ket_space, u_coef_cache) = <bra ||T^wT|| ket>_(rhoT)
// Sfunc(wS, rhoS, base_space, bra_space, ket_space, u_coef_cache) = <bra ||S^wS|| ket>_(rhoS)
//
basis::OperatorBlock<double> CoupledU3Operators(
  std::function<basis::OperatorBlock<double>(
    const u3::U3&,
    int,
    const sp3r::Sp3RSpace&,
    const sp3r::U3Subspace&,
    const sp3r::U3Subspace&,
    u3::UCoefCache&
    )> Tfunc,
  std::function<basis::OperatorBlock<double>(
    const u3::U3&,
    int,
    const sp3r::Sp3RSpace&,
    const sp3r::U3Subspace&,
    const sp3r::U3Subspace&,
    u3::UCoefCache&
    )> Sfunc,
  const u3::U3& omega_T,
  const u3::U3& omega_S,
  const u3::U3& omega0,
  int rho0,
  int outer_multiplicity,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
  )
{
  // Get sigma
  const u3::U3 sigma = base_space.sigma();

  // Subspace labels
  const auto& omega_bra = bra_subspace.omega();
  const auto& omega_ket = ket_subspace.omega();

  // Get labels to sum over
  // max allowed rho0_prime for coupling wS x wT -> w0
  // same as max allowed rho0 (as wT x wS isomorphic to wS x wT)
  int rho0prime_max = u3::OuterMultiplicity(omega_S,omega_T,omega0);

  // Create matrix of Sp(3,R)>U(3) branched rmes
  basis::OperatorBlock<double> TxS = basis::OperatorBlock<double>::Zero(
      bra_subspace.upsilon_max(),
      ket_subspace.upsilon_max()
      );

  //auto [Twrap, Swrap] = wrap_functions(Tfunc,Sfunc,omega_bar_vec,omega_T,omega_S,
  //                        omega_bra,base_space,bra_subspace,ket_subspace,u_coef_cache);

  // Calculate matrix elements
  // Iterate over omega = omega_ket, omega_prime = omega_bra
  // Only loop if w_ket x w0 -> w_bra AND wT x wS -> rho0 w0 is allowed
  if ((u3::OuterMultiplicity(omega_T,omega_S,omega0)>=rho0))
    if (u3::OuterMultiplicity(omega_ket,omega0,omega_bra)>0)
    {
      //
      for (int upsilon_bra = 1; upsilon_bra <= bra_subspace.upsilon_max(); upsilon_bra++)
      {
        for (int upsilon_ket = 1; upsilon_ket <= ket_subspace.upsilon_max(); upsilon_ket++)
        {
          // get rows and columns
          int row = upsilon_bra - 1;
          int col = upsilon_ket - 1;

          // sum over matrix elements of uncoupled operators
          double sum = 0;

          // Loop over all possible subspaces
          for (const auto& barred_subspace : base_space)
          {
            const auto& omega_bar = barred_subspace.omega();
            // get outer multiplicities corresponding to barred subspaces
            int rho_T_max = u3::OuterMultiplicity(omega_bar,omega_T,omega_bra);
            int rho_S_max = u3::OuterMultiplicity(omega_ket,omega_S,omega_bar);
            // loop over outer multiplicities, including only allowed subspaces
            for (int rho_S = 1; rho_S <= rho_S_max; rho_S++)
              for (int rho_T = 1; rho_T <= rho_T_max; rho_T++)
              {
                // compute and sum matrix elements of uncoupled operators
                basis::OperatorBlock<double> Tmat = Tfunc(omega_T,rho_T,base_space,bra_subspace,barred_subspace,u_coef_cache);
                basis::OperatorBlock<double> Smat = Sfunc(omega_S,rho_S,base_space,barred_subspace,ket_subspace,u_coef_cache);
                // matrix dot product to sum over barred subspace labels
                basis::OperatorBlock<double> TS_dot = Tmat.row(row)*Smat.col(col);
                // loop over rho0prime
                for (int rho0prime = 1; rho0prime <= rho0prime_max; rho0prime++)
                {
                  sum += u3::Phi(omega_T.SU3(),omega_S.SU3(),omega0.SU3(),rho0,rho0prime)
                    * u3::U(omega_ket.SU3(),omega_S.SU3(),omega_bra.SU3(),
                    omega_T.SU3(),omega_bar.SU3(),rho_S,rho_T,
                    omega0.SU3(),rho0prime,outer_multiplicity)
                    * (TS_dot(0,0));
                }
              }
          }
          TxS(row,col) = sum;
        }
      }
    }
  // return computed sector
  return TxS;
}


}