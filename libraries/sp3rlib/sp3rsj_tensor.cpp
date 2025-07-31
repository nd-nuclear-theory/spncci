/****************************************************************
  sp3rsj_tensor.cpp

  Colin V. Coane
  
  SPDX-License-Identifier: MIT 

  07/08/24 (cvc): Created.

****************************************************************/

#include "sp3rlib/sp3rsj_tensor.h"

namespace sp3r_operator
{

// Wrapper for J^2
basis::OperatorBlock<double> JdotJSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache,
  const HalfInt& S_bra,
  const HalfInt& S_ket,
  const HalfInt& J_bra,
  const HalfInt& J_ket
)
{
  // J value of ket subspace as a double
  double Jval = (J_ket.TwiceValue()/2.0);
  // Return J(J+1)*Identity
  if ((J_bra == J_ket) && (S_bra == S_ket))
  {
    //std::cout<<fmt::format("Return I_{}x{}",bra_subspace.dimension(),ket_subspace.dimension())<<std::endl;
    return (Jval)*(Jval+1.0)*IdentitySp3R(omega,rho,base_space,bra_subspace,ket_subspace,u_coef_cache);
  }
  else
  {
    //std::cout<<fmt::format("Return 0_{}x{}",bra_subspace.dimension(),ket_subspace.dimension())<<std::endl;
    return (0.0)*IdentitySp3R(omega,rho,base_space,bra_subspace,ket_subspace,u_coef_cache);
  }
}

// Wrapper for S^2
basis::OperatorBlock<double> SdotSSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache,
  const HalfInt& S_bra,
  const HalfInt& S_ket,
  const HalfInt& J_bra,
  const HalfInt& J_ket
)
{
  // S value of ket subspace as a double
  double Sval = (S_ket.TwiceValue()/2.0);
  // Return S(S+1)*Identity
  if ((J_bra == J_ket) && (S_bra == S_ket))
  {
    //std::cout<<fmt::format("Return I_{}x{}",bra_subspace.dimension(),ket_subspace.dimension())<<std::endl;
    return Sval*(Sval+1.0)*IdentitySp3R(omega,rho,base_space,bra_subspace,ket_subspace,u_coef_cache);
  }
  else
  {
    //std::cout<<fmt::format("Return 0_{}x{}",bra_subspace.dimension(),ket_subspace.dimension())<<std::endl;
    return (0.0)*IdentitySp3R(omega,rho,base_space,bra_subspace,ket_subspace,u_coef_cache);
  }
}

// Wrapper for LdotS = (1/2)*(J^2 - S^2 - L^2)
basis::OperatorBlock<double> LdotSSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache,
  const HalfInt& S_bra,
  const HalfInt& S_ket,
  const HalfInt& J_bra,
  const HalfInt& J_ket
)
{
  // operator matrix
  basis::OperatorBlock<double> LdotS_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );
  
  LdotS_matrix += (1.0/2.0)*JdotJSp3R(omega,rho,base_space,bra_subspace,ket_subspace,u_coef_cache,S_bra,S_ket,J_bra,J_ket);
  LdotS_matrix -= (1.0/2.0)*SdotSSp3R(omega,rho,base_space,bra_subspace,ket_subspace,u_coef_cache,S_bra,S_ket,J_bra,J_ket);
  LdotS_matrix -= (1.0/2.0)*LdotLSp3R(omega,rho,base_space,bra_subspace,ket_subspace,u_coef_cache);

  return LdotS_matrix;
}

}
