/****************************************************************
  u3_tensor.cpp

  Colin V. Coane
  
  SPDX-License-Identifier: MIT 

  07/02/24 (cvc): Created. Functions moved from u3_racah_product.h.
  07/02/24 (cvc): Fixed issue with multiplicity of SU3GeneratorSp3R().
  07/03/24 (cvc): Q dot Q implemented for mass quadrupole.

****************************************************************/

#include "sp3rlib/u3_tensor.h"

namespace sp3r_operator
{

// Wrapper for Identity
basis::OperatorBlock<double> IdentitySp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  const u3::U3 sigma = base_space.sigma();
  const u3::U3 omega_Id(0,{0u,0u});
  // Subspace labels
  const auto& omega_bra = bra_subspace.omega();
  const auto& omega_ket = ket_subspace.omega();

  // Return Identity matrix elements for given sector
  //std::cout<<fmt::format("Return Identity")<<std::endl;
  // Check for coupling to (0,0) ???
  //if (omega_bra == omega_ket)
  if ((omega_bra == omega_ket) && (omega.SU3() == omega_Id.SU3()))
  {
    //std::cout<<fmt::format("Return I_{}x{}",bra_subspace.dimension(),ket_subspace.dimension())<<std::endl;
    return basis::OperatorBlock<double>::Identity(
      bra_subspace.dimension(),
      ket_subspace.dimension()
    );
  }
  else
  {
    //std::cout<<fmt::format("Return 0_{}x{}",bra_subspace.dimension(),ket_subspace.dimension())<<std::endl;
    return basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
    );
  }
}


// Wrapper for symplectic raising and lowering operators
// laddering between states in w', w subspaces of Sp(3,R) irrep defined by base_space
// Raising Operator, A^(2,0)
basis::OperatorBlock<double> RaisingOperatorSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  const u3::U3 sigma = base_space.sigma();
  // Ask about doing check vs just returning regardless of omega
  u3::U3 omega_A(2,{2u,0u});
  if (omega.SU3() == omega_A.SU3())
  {
    // Return raising operator A^(2,0) for given sector
    return sp3r::Sp3rRaisingOperator(sigma,bra_subspace,ket_subspace,u_coef_cache);
  }
  else
  {
    return basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
    );
  }
}


// Lowering Operator, B^(0,2)
basis::OperatorBlock<double> LoweringOperatorSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  const u3::U3 sigma = base_space.sigma();
  // Ask about doing check vs just returning regardless of omega
  u3::U3 omega_B(-2,{0u,2u});
  if (omega.SU3() == omega_B.SU3())
  {
    // Return lowering operator B^(0,2) for given sector
    return sp3r::Sp3rLoweringOperator(sigma,bra_subspace,ket_subspace,u_coef_cache);
  }
  else
  {
    return basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
    );
  }
}

// SU(3) Generators, C^(1,1)
basis::OperatorBlock<double> SU3GeneratorSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  const u3::U3 sigma = base_space.sigma();
  // Ask about doing check vs just returning regardless of omega
  u3::U3 omega_C(0,{1u,1u});
  if ((omega.SU3() == omega_C.SU3()) && (rho == 1))
  {
    // Return SU(3) generator C^(1,1) for given sector
    return sp3r::SU3Generator(sigma,bra_subspace,ket_subspace);
  }
  else
  {
    return basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
    );
  }
}


// Harmonic Oscillator Hamiltonian, H_0^(0,0)
basis::OperatorBlock<double> HarmonicOscillatorSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  const u3::U3 sigma = base_space.sigma();
  // Subspace labels
  const auto& omega_bra = bra_subspace.omega();
  const auto& omega_ket = ket_subspace.omega();
  // Ask about doing check vs just returning regardless of omega
  u3::U3 omega_H(0,{0u,0u});
  if ((omega.SU3() == omega_H.SU3()) && (omega_bra == omega_ket))
  {
    // Return Harmonic Oscillator Quanta for given sector
    return (omega_bra.N().TwiceValue()/2.0)*basis::OperatorBlock<double>::Identity(
      bra_subspace.dimension(),
      ket_subspace.dimension()
    );
  }
  else
  {
    return basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
    );
  }
}

// Mass Quadrupole Operator Q_2 = sqrt(3)*[C^(1,1) + A^(2,0) + B^(0,2)]_2
basis::OperatorBlock<double> MassQuadrupoleSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  // Prefactor
  // pi = 2*acos(0.0)
  double prefactor = std::sqrt(3.0);

  // operator matrix
  basis::OperatorBlock<double> Q_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );

  Q_matrix += prefactor*RaisingOperatorSp3R(omega,rho,base_space,bra_subspace,ket_subspace,u_coef_cache);
  Q_matrix += prefactor*LoweringOperatorSp3R(omega,rho,base_space,bra_subspace,ket_subspace,u_coef_cache);
  Q_matrix += prefactor*SU3GeneratorSp3R(omega,rho,base_space,bra_subspace,ket_subspace,u_coef_cache);

  return Q_matrix;
}


// SU(3) Casimir
basis::OperatorBlock<double> SU3CasimirSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  const u3::U3 sigma = base_space.sigma();
  // Subspace labels
  const auto& omega_bra = bra_subspace.omega();
  const auto& omega_ket = ket_subspace.omega();

  // Return matrix elements for given sector
  if (omega_bra == omega_ket)
  {
    double casimir_rme = u3::Casimir2(omega_bra.SU3());
    return basis::OperatorBlock<double>::Identity(
      bra_subspace.dimension(),
      ket_subspace.dimension()
    )*casimir_rme;
  }
  else
  {
    return  basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
    );
  }
}


///////////////////////////////////////////
//////////// COUPLED OPERATORS ////////////
///////////////////////////////////////////


// Coupled SU(3) Generators [C^(1,1) x C^(1,1)]
basis::OperatorBlock<double> CtimesCSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  // U(3) Labels
  u3::U3 omega_C(0,{1u,1u});
  int rho0 = 1;

  // operator matrix
  basis::OperatorBlock<double> CxC_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );
  

  CxC_matrix += CoupledU3Operators(&SU3GeneratorSp3R,&SU3GeneratorSp3R,
                    omega_C,omega_C,omega,rho0,rho,base_space,
                    bra_subspace,ket_subspace,u_coef_cache);
  
  return CxC_matrix;
}

// Coupled Raising Operators [A^(2,0) x A^(2,0)]
basis::OperatorBlock<double> AtimesASp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  // U(3) Labels
  u3::U3 omega_A(2,{2u,0u});
  int rho0 = 1;

  // operator matrix
  basis::OperatorBlock<double> AxA_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );
  

  AxA_matrix += CoupledU3Operators(&RaisingOperatorSp3R,&RaisingOperatorSp3R,
                    omega_A,omega_A,omega,rho0,rho,base_space,
                    bra_subspace,ket_subspace,u_coef_cache);
  
  return AxA_matrix;
}

// Coupled Lowering Operators [B^(0,2) x B^(0,2)]
basis::OperatorBlock<double> BtimesBSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  // U(3) Labels
  u3::U3 omega_B(-2,{0u,2u});
  int rho0 = 1;

  // operator matrix
  basis::OperatorBlock<double> BxB_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );
  

  BxB_matrix += CoupledU3Operators(&LoweringOperatorSp3R,&LoweringOperatorSp3R,
                    omega_B,omega_B,omega,rho0,rho,base_space,
                    bra_subspace,ket_subspace,u_coef_cache);
  
  return BxB_matrix;
}

// Coupled Symplectic Raising/Lowering Operators [A^(2,0) x B^(0,2)]
basis::OperatorBlock<double> AtimesBSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  // U(3) Labels
  u3::U3 omega_A(2,{2u,0u});
  u3::U3 omega_B(-2,{0u,2u});
  int rho0 = 1;

  // operator matrix
  basis::OperatorBlock<double> AxB_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );
  

  AxB_matrix += CoupledU3Operators(&RaisingOperatorSp3R,&LoweringOperatorSp3R,
                    omega_A,omega_B,omega,rho0,rho,base_space,
                    bra_subspace,ket_subspace,u_coef_cache);
  
  return AxB_matrix;
}


// Coupled Symplectic Lowering/Raising Operators [B^(0,2) x A^(2,0)]
basis::OperatorBlock<double> BtimesASp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  // U(3) Labels
  u3::U3 omega_A(2,{2u,0u});
  u3::U3 omega_B(-2,{0u,2u});
  int rho0 = 1;

  // operator matrix
  basis::OperatorBlock<double> BxA_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );
  

  BxA_matrix += CoupledU3Operators(&LoweringOperatorSp3R,&RaisingOperatorSp3R,
                    omega_B,omega_A,omega,rho0,rho,base_space,
                    bra_subspace,ket_subspace,u_coef_cache);
  
  return BxA_matrix;
}


// Coupled Generator/Raising Operators [C^(1,1) x A^(2,0)]
basis::OperatorBlock<double> CtimesASp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  // U(3) Labels
  u3::U3 omega_A(2,{2u,0u});
  u3::U3 omega_C(0,{1u,1u});
  int rho0 = 1;

  // operator matrix
  basis::OperatorBlock<double> CxA_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );
  

  CxA_matrix += CoupledU3Operators(&SU3GeneratorSp3R,&RaisingOperatorSp3R,
                    omega_C,omega_A,omega,rho0,rho,base_space,
                    bra_subspace,ket_subspace,u_coef_cache);
  
  return CxA_matrix;
}

// Coupled Raising/Generator Operators [A^(2,0) x C^(1,1)]
basis::OperatorBlock<double> AtimesCSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  // U(3) Labels
  u3::U3 omega_A(2,{2u,0u});
  u3::U3 omega_C(0,{1u,1u});
  int rho0 = 1;

  // operator matrix
  basis::OperatorBlock<double> AxC_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );
  

  AxC_matrix += CoupledU3Operators(&RaisingOperatorSp3R,&SU3GeneratorSp3R,
                    omega_A,omega_C,omega,rho0,rho,base_space,
                    bra_subspace,ket_subspace,u_coef_cache);
  
  return AxC_matrix;
}


// Coupled Generator/Lowering Operators [C^(1,1) x B^(0,2)]
basis::OperatorBlock<double> CtimesBSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  // U(3) Labels
  u3::U3 omega_B(-2,{0u,2u});
  u3::U3 omega_C(0,{1u,1u});
  int rho0 = 1;

  // operator matrix
  basis::OperatorBlock<double> CxB_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );
  

  CxB_matrix += CoupledU3Operators(&SU3GeneratorSp3R,&LoweringOperatorSp3R,
                    omega_C,omega_B,omega,rho0,rho,base_space,
                    bra_subspace,ket_subspace,u_coef_cache);
  
  return CxB_matrix;
}

// Coupled Lowering/Generator Operators [B^(0,2) x C^(1,1)]
basis::OperatorBlock<double> BtimesCSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  // U(3) Labels
  u3::U3 omega_B(-2,{0u,2u});
  u3::U3 omega_C(0,{1u,1u});
  int rho0 = 1;

  // operator matrix
  basis::OperatorBlock<double> BxC_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );
  

  BxC_matrix += CoupledU3Operators(&LoweringOperatorSp3R,&SU3GeneratorSp3R,
                    omega_B,omega_C,omega,rho0,rho,base_space,
                    bra_subspace,ket_subspace,u_coef_cache);
  
  return BxC_matrix;
}



// Angular momentum operator L^2 = -sqrt(3)*[C_1 x C_1]_00
// L^2 = alpha*[C^(1,1) x C^(1,1)]^(0,0) + beta*[C^(1,1) x C^(1,1)]^(2,2)
basis::OperatorBlock<double> LdotLSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  // U(3) Labels
  u3::U3 omega_0(0,{0u,0u});
  u3::U3 omega_2(0,{2u,2u});
  int rho0 = 1;

  // Coefficient from dot product
  double dot_coef = (-1.0)*std::sqrt(3.0);
  // Coefficients from recoupling to good SU(3)
  double alpha_00 = (-1.0)*std::sqrt(3.0/8.0);
  double alpha_22 = (-1.0)*std::sqrt(5.0/8.0);

  // operator matrix
  basis::OperatorBlock<double> LxL_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );


  // return [C^(1,1) x C^(1,1)]^omega components of L^2
  if ((omega.SU3() == omega_0.SU3()))
  {
    LxL_matrix += dot_coef*alpha_00*CtimesCSp3R(omega,rho,base_space,bra_subspace,ket_subspace,u_coef_cache);
  }
  else if ((omega.SU3() == omega_2.SU3()))
  {
    LxL_matrix += dot_coef*alpha_22*CtimesCSp3R(omega,rho,base_space,bra_subspace,ket_subspace,u_coef_cache);
  }
  
  // return operator
  return LxL_matrix;
}



// Elliott Quadrupole operator QdotQ = 3*sqrt(5)*[C_2 x C_2]_00
basis::OperatorBlock<double> QdotQElliottSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  // U(3) Labels
  u3::U3 omega_0(0,{0u,0u});
  u3::U3 omega_2(0,{2u,2u});
  int rho0 = 1;

  // Coefs
  // Coefficient from dot product
  double dot_coef = (3.0)*std::sqrt(5.0);
  // Coefficients from recoupling to good SU(3)
  double alpha_00 = (1.0)*std::sqrt(5.0/8.0);
  double alpha_22 = (-1.0)*std::sqrt(3.0/8.0);

  // operator matrix
  basis::OperatorBlock<double> QxQElliott_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );


  // return [C^(1,1) x C^(1,1)]^omega components of Elliott QdotQ
  if ((omega.SU3() == omega_0.SU3()))
  {
    QxQElliott_matrix += dot_coef*alpha_00*CtimesCSp3R(omega,rho,base_space,bra_subspace,ket_subspace,u_coef_cache);
  }
  else if ((omega.SU3() == omega_2.SU3()))
  {
    QxQElliott_matrix += dot_coef*alpha_22*CtimesCSp3R(omega,rho,base_space,bra_subspace,ket_subspace,u_coef_cache);
  }
  
  // return operator
  return QxQElliott_matrix;
}



// Sp(3,R) Mass Quadrupole operator QdotQ = 3*sqrt(5)*[(C_2 + A_2 + B_2) x (C_2 + A_2 + B_2)]_00
basis::OperatorBlock<double> QdotQSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  // U(3) Labels
  u3::U3 omega_20(2,{2u,0u});
  u3::U3 omega_02(-2,{0u,2u});
  u3::U3 omega_00(0,{0u,0u});
  u3::U3 omega_22(0,{2u,2u});
  u3::U3 omega_40(4,{4u,0u});
  u3::U3 omega_04(-4,{0u,4u});
  int rho0 = 1;

  // Coefs
  // Coefficient from dot product
  double dot_coef = (3.0)*std::sqrt(5.0);

  // operator matrix
  basis::OperatorBlock<double> QxQ_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );


  // return [C^(1,1) x C^(1,1)]^omega components of L^2
  if ((omega.SU3() == omega_00.SU3()))
  {
    QxQ_matrix += dot_coef*(1.0)*std::sqrt(5.0/8.0)*CtimesCSp3R(omega,rho,base_space,bra_subspace,ket_subspace,u_coef_cache);
    QxQ_matrix += dot_coef*(2.0)*std::sqrt(5.0/6.0)*AtimesBSp3R(omega,rho,base_space,bra_subspace,ket_subspace,u_coef_cache);
    QxQ_matrix += dot_coef*std::sqrt(20.0/9.0)*HarmonicOscillatorSp3R(omega,rho,base_space,bra_subspace,ket_subspace,u_coef_cache);
  }
  else if ((omega.SU3() == omega_22.SU3()))
  {
    QxQ_matrix += dot_coef*(-1.0)*std::sqrt(3.0/8.0)*CtimesCSp3R(omega,rho,base_space,bra_subspace,ket_subspace,u_coef_cache);
    QxQ_matrix += dot_coef*(2.0)*std::sqrt(1.0/6.0)*AtimesBSp3R(omega,rho,base_space,bra_subspace,ket_subspace,u_coef_cache);
  }
  else if ((omega.SU3() == omega_02.SU3()))
  {
    QxQ_matrix += dot_coef*std::sqrt(5.0/9.0)*AtimesASp3R(omega,rho,base_space,bra_subspace,ket_subspace,u_coef_cache);
    QxQ_matrix += dot_coef*(2.0)*CtimesBSp3R(omega,rho,base_space,bra_subspace,ket_subspace,u_coef_cache);
    QxQ_matrix += dot_coef*std::sqrt(40.0/3.0)*LoweringOperatorSp3R(omega,rho,base_space,bra_subspace,ket_subspace,u_coef_cache);
  }
  else if ((omega.SU3() == omega_20.SU3()))
  {
    QxQ_matrix += dot_coef*std::sqrt(5.0/9.0)*BtimesBSp3R(omega,rho,base_space,bra_subspace,ket_subspace,u_coef_cache);
    QxQ_matrix += dot_coef*(2.0)*AtimesCSp3R(omega,rho,base_space,bra_subspace,ket_subspace,u_coef_cache);
    QxQ_matrix += dot_coef*std::sqrt(40.0/3.0)*RaisingOperatorSp3R(omega,rho,base_space,bra_subspace,ket_subspace,u_coef_cache);
  }
  else if ((omega.SU3() == omega_40.SU3()))
  {
    QxQ_matrix += dot_coef*(2.0/3.0)*AtimesASp3R(omega,rho,base_space,bra_subspace,ket_subspace,u_coef_cache);
  }
  else if ((omega.SU3() == omega_04.SU3()))
  {
    QxQ_matrix += dot_coef*(2.0/3.0)*BtimesBSp3R(omega,rho,base_space,bra_subspace,ket_subspace,u_coef_cache);
  }
  
  // return operator
  return QxQ_matrix;
}


}

