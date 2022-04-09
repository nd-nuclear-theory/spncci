/****************************************************************
 sp3r_operators.cpp

  Anna E. McCoy
  University of Notre Dame

  SPDX-License-Identifier: MIT
****************************************************************/
#include "sp3rlib/sp3r_operator.h"



namespace sp3r
{
  basis::OperatorBlock<double> Sp3rRaisingOperator(
    const u3::U3& sigma,
    const sp3r::U3Subspace& bra_subspace,
    const sp3r::U3Subspace& ket_subspace,
    u3::UCoefCache& u_coef_cache
    )
  {
    const auto& omega_bra = bra_subspace.omega();
    const auto& omega_ket = ket_subspace.omega();
    if(!u3::OuterMultiplicity(omega_ket,{2,{2u,0u}},omega_bra))
      return basis::OperatorBlock<double>::Zero(bra_subspace.upsilon_max(),ket_subspace.upsilon_max());

    const auto& bra_nonorthogonal_basis = bra_subspace.nonorthogonal_basis();
    const auto& ket_nonorthogonal_basis = ket_subspace.nonorthogonal_basis();

    basis::OperatorBlock<double> A_boson_matrix =
      u3boson::U3BosonCreationOperator(sigma,bra_nonorthogonal_basis,ket_nonorthogonal_basis,u_coef_cache);

    return bra_subspace.K_matrix()*A_boson_matrix*ket_subspace.Kinv_matrix();
  }


  basis::OperatorBlock<double> Sp3rLoweringOperator(
    const u3::U3& sigma,
    const sp3r::U3Subspace& bra_subspace,
    const sp3r::U3Subspace& ket_subspace,
    u3::UCoefCache& u_coef_cache
    )
  {
    const auto& x_bra=bra_subspace.omega().SU3();
    const auto& x_ket=ket_subspace.omega().SU3();

    double conjugation_factor
      =ParitySign(u3::ConjugationGrade(x_bra)+u3::ConjugationGrade(x_ket))
        *std::sqrt(1.*u3::dim(x_ket)/u3::dim(x_bra));

    return conjugation_factor*Sp3rRaisingOperator(sigma,ket_subspace,bra_subspace,u_coef_cache).transpose();

  }



  basis::OperatorBlock<double> SU3Generator(
    const u3::U3& sigma,
    const sp3r::U3Subspace& bra_subspace,
    const sp3r::U3Subspace& ket_subspace
    )
  {
    const auto& omega_bra = bra_subspace.omega();
    const auto& omega_ket = ket_subspace.omega();
    if(omega_bra==omega_ket)
      {
        //TODO generalize phase
        int phase = omega_bra.SU3().mu()==0?1:-1;
        double rme=phase*std::sqrt(2.0*Casimir2(omega_bra.SU3()));
        return rme*basis::OperatorBlock<double>::Identity(bra_subspace.upsilon_max(),ket_subspace.upsilon_max());
      }
    else
      return basis::OperatorBlock<double>::Zero(bra_subspace.upsilon_max(),ket_subspace.upsilon_max());
  }


  Eigen::MatrixXd  U3Operator(
      const sp3r::Sp3RSpace& sp3r_space, 
      const u3::U3& omegap, 
      const u3::U3& omega
    )
  {
    // Look up  relevant subspaces in Sp(3,R) irrep
    int subspace_index_bra=sp3r_space.LookUpSubspaceIndex(omegap);
    int subspace_index_ket=sp3r_space.LookUpSubspaceIndex(omega);

    const sp3r::U3Subspace& subspace_bra =  sp3r_space.GetSubspace(subspace_index_bra);
    const sp3r::U3Subspace& subspace_ket =  sp3r_space.GetSubspace(subspace_index_ket);

    //Construct boson matrix 
    int vp_max=subspace_bra.size();
    int v_max=subspace_ket.size();

    if(omegap==omega)
      {
        double Crme=sqrt(2.0*Casimir2(omega.SU3()));
        return Crme*Eigen::MatrixXd::Identity(v_max,v_max);
      }

    else
      return Eigen::MatrixXd::Zero(vp_max,v_max);
  }


}
