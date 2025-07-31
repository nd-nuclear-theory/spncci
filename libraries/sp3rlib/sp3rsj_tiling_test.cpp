/****************************************************************
  sp3rsj_tiling_test.cpp

  Colin V. Coane
  
  SPDX-License-Identifier: MIT 

  06/05/24 (cvc): Created.
  07/03/24 (cvc): Modified tests.
  07/08/24 (cvc): Modified tests.

****************************************************************/

#include "sp3rlib/sp3r.h"
#include "sp3rlib/sp3rsj.h"
#include "sp3rlib/sp3r_operator.h"
#include "sp3rlib/u3_racah_product.h"
#include "sp3rlib/sp3rsj_tiling.h"
#include "sp3rlib/u3_tensor.h"
#include "sp3rlib/sp3rsj_tensor.h"

#include <iostream>
#include "sp3rlib/u3coef.h"
#include "sp3rlib/u3.h"
#include "fmt/format.h"
#include "am/halfint.h"
#include "am/am.h"


// main
int main(int argc, char **argv)
{
  u3::U3CoefInit(39);

  ////////////////////////////////////////////////////////////////
  // Sp(3,R) x SU(2) operator construction test
  ////////////////////////////////////////////////////////////////



  u3::U3 sigma(16,{2u,1u});
  //u3::U3 sigma(5,{2u,0u});
  //u3::U3 sigma({29u,2u},{2u,1u});

  // List of U(3) operator labels
  std::vector<u3::U3>
  op_list = {
      u3::U3(2,{2u,0u}),
      u3::U3(-2,{0u,2u}),
      u3::U3(0,{0u,0u}),
      //u3::U3(0,{1u,1u}),
      u3::U3(0,{2u,2u}),
      u3::U3(4,{4u,0u}),
      u3::U3(-4,{0u,4u})
    };

  // Full QxQxQ space labels
  // CHECK IF CORRECT
  std::vector<u3::U3> QxQ_su3_list = {
      u3::U3(2,{2u,0u}),
      u3::U3(-2,{0u,2u}),
      u3::U3(0,{0u,0u}),
      u3::U3(0,{1u,1u}),
      u3::U3(0,{2u,2u}),
      u3::U3(0,{0u,3u}),
      u3::U3(0,{3u,0u}),
      u3::U3(2,{0u,1u}),
      u3::U3(2,{1u,2u}),
      u3::U3(2,{3u,1u}),
      u3::U3(-2,{1u,0u}),
      u3::U3(-2,{1u,3u}),
      u3::U3(-2,{2u,1u}),
      u3::U3(4,{4u,0u}),
      u3::U3(-4,{0u,4u})
    };
  

  int Nn_max=2;
  //int Nn_max=4;
  u3::UCoefCache u_coef_cache;
  int rho0 = 1;
  int outer_multiplicity = 1;
  u3::U3 omega_Id(0,{0u,0u});
  u3::U3 omega0(0,{0u,0u});
  u3::U3 omega1(0,{1u,1u});
  u3::U3 omega2(0,{2u,2u});
  // check these labels
  u3::U3 omega_A(2,{2u,0u});
  u3::U3 omega_B(-2,{0u,2u});
  u3::U3 omega_C(0,{1u,1u});

  HalfInt S_space = HalfInt(4,2);
  HalfInt J_space = HalfInt(4,2);

  HalfInt J_space_alt = HalfInt(6,2);

  // Sp3RSpace
  sp3r::Sp3RSpace sp3r_space(sigma,Nn_max);
  // Sp3RSJ Spaces
  sp3r::Sp3RSJSpace bra_space(sigma,Nn_max,S_space,J_space);
  sp3r::Sp3RSJSpace ket_space(sigma,Nn_max,S_space,J_space);

  std::cout<<ket_space.DebugStr()<<std::endl;



  HalfInt S_space_l1 = HalfInt(0,2);
  HalfInt J_space_l1 = HalfInt(2,2);

  HalfInt S_space_l2 = HalfInt(0,2);
  HalfInt J_space_l2 = HalfInt(4,2);

  // Sp3RSJ Spaces
  sp3r::Sp3RSJSpace bra_space_l1(sigma,Nn_max,S_space_l1,J_space_l1);
  sp3r::Sp3RSJSpace ket_space_l1(sigma,Nn_max,S_space_l1,J_space_l1);
  sp3r::Sp3RSJSpace bra_space_l2(sigma,Nn_max,S_space_l2,J_space_l2);
  sp3r::Sp3RSJSpace ket_space_l2(sigma,Nn_max,S_space_l2,J_space_l2);

  // Operator Spaces

  // J = 0 Operator
  HalfInt S_op = HalfInt(0,2);
  HalfInt J_op = HalfInt(0,2);
  sp3r::OperatorSpaceSp3RSJ op_space_J0(S_op,J_op,QxQ_su3_list);
  std::cout<<op_space_J0.DebugStr()<<std::endl;


  // List of U(3) operator labels
  std::vector<u3::U3>
  op_list_CxC = {
      u3::U3(0,{0u,0u}),
      u3::U3(0,{2u,2u}),
    };
  
  sp3r::OperatorSpaceSp3RSJ op_space_CxC(S_op,J_op,op_list_CxC);


  // List of U(3) operator labels
  std::vector<u3::U3>
  op_list_CxC00 = {
      u3::U3(0,{0u,0u}),
    };
  
  sp3r::OperatorSpaceSp3RSJ op_space_CxC00(S_op,J_op,op_list_CxC00);

  std::vector<u3::U3>
  op_list_CxC22 = {
      u3::U3(0,{2u,2u}),
    };
  
  sp3r::OperatorSpaceSp3RSJ op_space_CxC22(S_op,J_op,op_list_CxC22);


  std::vector<u3::U3>
  op_list_BA = {
      u3::U3(0,{0u,0u}),
    };
  
  sp3r::OperatorSpaceSp3RSJ op_space_BA(S_op,J_op,op_list_BA);

  // J = 2 Operator
  HalfInt J_op_2 = HalfInt(4,2);
  sp3r::OperatorSpaceSp3RSJ op_space_J2(S_op,J_op_2,QxQ_su3_list);

  // Class for observables desired

  // Identity
  sp3r::OperatorInfoSp3RSJ id_ob_op(op_space_J0,&sp3r_operator::IdentitySp3R);
  sp3r::ObservableSp3RSJ id_obs = {
    std::make_pair(id_ob_op,1.0)
  };

  // Harmonic Oscillator Hamiltonian
  sp3r::OperatorInfoSp3RSJ harmonic_oscillator_op(op_space_J0,&sp3r_operator::HarmonicOscillatorSp3R);
  sp3r::ObservableSp3RSJ harmonic_oscillator_obs = {
    std::make_pair(harmonic_oscillator_op,1.0)
  };

  // Mass Quadrupole Operator
  sp3r::OperatorInfoSp3RSJ Q_op(op_space_J2,&sp3r_operator::MassQuadrupoleSp3R);
  sp3r::ObservableSp3RSJ Q_obs = {
    std::make_pair(Q_op,1.0)
  };

  // [CxC]_00
  sp3r::OperatorInfoSp3RSJ CdotC_op(op_space_J0,&sp3r_operator::CtimesCSp3R);
  sp3r::ObservableSp3RSJ CdotC_obs = {
    std::make_pair(CdotC_op,1.0)
  };

  // [AxB]_00
  sp3r::OperatorInfoSp3RSJ AdotB_op(op_space_BA,&sp3r_operator::AtimesBSp3R);
  sp3r::ObservableSp3RSJ AdotB_obs = {
    std::make_pair(AdotB_op,1.0)
  };

  // [CxC]_00
  sp3r::OperatorInfoSp3RSJ CxC_op(op_space_CxC,&sp3r_operator::CtimesCSp3R);

  // [BxA]_00
  sp3r::OperatorInfoSp3RSJ BdotA_op(op_space_BA,&sp3r_operator::BtimesASp3R);

  // L^2
  sp3r::OperatorInfoSp3RSJ LdotL_op(op_space_CxC,&sp3r_operator::LdotLSp3R);
  sp3r::ObservableSp3RSJ LdotL_obs = {
    std::make_pair(LdotL_op,1.0)
  };

  // Elliott Quadrupole
  sp3r::OperatorInfoSp3RSJ elliott_op(op_space_CxC,&sp3r_operator::QdotQElliottSp3R);
  sp3r::ObservableSp3RSJ elliott_obs = {
    std::make_pair(elliott_op,1.0)
  };

  // SU(3) Generator
  sp3r::OperatorInfoSp3RSJ su3_gen(op_space_J0,&sp3r_operator::SU3GeneratorSp3R);


  // J^2
  sp3r::OperatorInfoSp3RSJ J2_op(op_space_J0,&sp3r_operator::JdotJSp3R);
  sp3r::ObservableSp3RSJ J2_obs = {
    std::make_pair(J2_op,1.0)
  };

  // S^2
  sp3r::OperatorInfoSp3RSJ S2_op(op_space_J0,&sp3r_operator::SdotSSp3R);
  sp3r::ObservableSp3RSJ S2_obs = {
    std::make_pair(S2_op,1.0)
  };

  // LdotS as single operator
  sp3r::OperatorInfoSp3RSJ LS_op(op_space_CxC,&sp3r_operator::LdotSSp3R);
  sp3r::ObservableSp3RSJ LS_obs = {
    std::make_pair(LS_op,1.0)
  };

  // LdotS as linear combination of operators
  sp3r::ObservableSp3RSJ LS_obs_2 = {
    std::make_pair(J2_op,0.5),
    std::make_pair(LdotL_op,-0.5),
    std::make_pair(S2_op,-0.5),
  };

  // Mass Quadrupole
  sp3r::OperatorInfoSp3RSJ QdotQ_op(op_space_CxC,&sp3r_operator::QdotQSp3R);
  sp3r::ObservableSp3RSJ QdotQ_obs = {
    std::make_pair(QdotQ_op,1.0)
  };

  // Construct observable
  basis::OperatorBlock<double> id_matrix;
  sp3r::ConstructSp3ROperatorMatrix(sp3r_space,bra_space,ket_space,id_obs,id_matrix);

  std::cout<<bra_space.GetBranchedDimension()<<std::endl;
  std::cout<<ket_space.GetBranchedDimension()<<std::endl;
  
  // Comparison with true identity matrix
  basis::OperatorBlock<double> id_test = basis::OperatorBlock<double>::Identity(
      bra_space.GetBranchedDimension(),
      ket_space.GetBranchedDimension()
    );
  
  if (!mcutils::IsZero(id_matrix-id_test,1e-8))
  {
    std::cout<<fmt::format("Incorrect Identity Matrix")<<std::endl;
    std::cout<<id_matrix<<std::endl;
  }
  else
  {
    std::cout<<fmt::format("Correct Identity Matrix")<<std::endl;
    std::cout<<id_test<<std::endl;
  }


  // Construct observable
  basis::OperatorBlock<double> ho_matrix;
  sp3r::ConstructSp3ROperatorMatrix(sp3r_space,bra_space,ket_space,harmonic_oscillator_obs,ho_matrix);

  std::cout<<fmt::format("Harmonic Oscillator Hamiltonian")<<std::endl;
  std::cout<<ho_matrix<<std::endl;


  // Construct observable
  basis::OperatorBlock<double> Q_matrix;
  sp3r::ConstructSp3ROperatorMatrix(sp3r_space,bra_space,ket_space,Q_obs,Q_matrix);

  std::cout<<fmt::format("Mass Quadrupole Operator")<<std::endl;
  std::cout<<Q_matrix<<std::endl;


  // Construct observable
  basis::OperatorBlock<double> CxC_matrix;
  sp3r::ConstructSp3ROperatorMatrix(sp3r_space,bra_space,ket_space,CdotC_obs,CxC_matrix);

  std::cout<<fmt::format("QdotQ Elliott = [CxC]_00 Operator")<<std::endl;
  std::cout<<CxC_matrix<<std::endl;


  // Construct observable
  basis::OperatorBlock<double> LxL_matrix;
  sp3r::ConstructSp3ROperatorMatrix(sp3r_space,bra_space,ket_space,LdotL_obs,LxL_matrix);

  std::cout<<fmt::format("LdotL Operator")<<std::endl;
  std::cout<<LxL_matrix<<std::endl;


  // Construct observable
  basis::OperatorBlock<double> elliott_matrix;
  sp3r::ConstructSp3ROperatorMatrix(sp3r_space,bra_space,ket_space,elliott_obs,elliott_matrix);

  std::cout<<fmt::format("QdotQ Elliott = [C_2xC_2]_00 Operator")<<std::endl;
  std::cout<<elliott_matrix<<std::endl;

  std::cout<<ket_space.DebugStr()<<std::endl;

  // Construct observable
  basis::OperatorBlock<double> ab_matrix;
  sp3r::ConstructSp3ROperatorMatrix(sp3r_space,bra_space,ket_space,AdotB_obs,ab_matrix);

  std::cout<<fmt::format("[AxB]^(0,0) Operator")<<std::endl;
  std::cout<<ab_matrix<<std::endl;


  // Construct observable
  basis::OperatorBlock<double> J2_matrix;
  sp3r::ConstructSp3ROperatorMatrix(sp3r_space,bra_space,ket_space,J2_obs,J2_matrix);

  std::cout<<fmt::format("J^2 Operator")<<std::endl;
  std::cout<<J2_matrix<<std::endl;


  // Construct observable
  basis::OperatorBlock<double> LS_mat_1;
  sp3r::ConstructSp3ROperatorMatrix(sp3r_space,bra_space,ket_space,LS_obs,LS_mat_1);

  std::cout<<fmt::format("LxS Operator")<<std::endl;
  std::cout<<LS_mat_1<<std::endl;

  // Construct observable
  basis::OperatorBlock<double> LS_mat_2;
  sp3r::ConstructSp3ROperatorMatrix(sp3r_space,bra_space,ket_space,LS_obs_2,LS_mat_2);

  std::cout<<fmt::format("LxS Operator")<<std::endl;
  std::cout<<LS_mat_2<<std::endl;


  if (!mcutils::IsZero(LS_mat_1-LS_mat_2,1e-8))
  {
    std::cout<<fmt::format("Incorrect L-S coupling")<<std::endl;
  }
  else
  {
    std::cout<<fmt::format("Correct L-S coupling")<<std::endl;
  }

  // Construct observable
  basis::OperatorBlock<double> QxQ_matrix;
  sp3r::ConstructSp3ROperatorMatrix(sp3r_space,bra_space,ket_space,QdotQ_obs,QxQ_matrix);

  std::cout<<fmt::format("QdotQ = [Q_2xQ_2]_00 Operator")<<std::endl;
  std::cout<<QxQ_matrix<<std::endl;



  // Test QdotQ for various irreps
  u3::U3 sigma1(16,{2u,1u});
  u3::U3 sigma2(5,{2u,0u});

  std::vector<u3::U3>
  sigma_list_Q = {
      u3::U3(16,{2u,1u}),
      u3::U3({29,2},{2u,1u}),
      u3::U3(5,{2u,0u}),
      u3::U3({3,2},{0u,0u}),
      u3::U3(16,{4u,3u})
    };

  HalfInt S_space_Q = HalfInt(20,2);
  HalfInt J_space_Q = HalfInt(20,2);

  for(const auto& sigma_Q : sigma_list_Q)
  {
    // Sp3RSpace
    sp3r::Sp3RSpace sp3r_space_Q(sigma_Q,Nn_max);
    // Sp3RSJ Spaces
    sp3r::Sp3RSJSpace bra_space_Q(sigma_Q,Nn_max,S_space_Q,J_space_Q);
    sp3r::Sp3RSJSpace ket_space_Q(sigma_Q,Nn_max,S_space_Q,J_space_Q);

    // Construct observable
    basis::OperatorBlock<double> QxQ_matrix_Q;
    //sp3r::ConstructSp3ROperatorMatrix(sp3r_space_Q,bra_space_Q,ket_space_Q,op_space_J0,QdotQ_obs,QxQ_matrix_Q);

    //std::cout<<fmt::format("QdotQ = [QxQ]_00 Operator ")<<std::endl;
    //std::cout<<QxQ_matrix_Q<<std::endl;
    //std::cout<<sp3r_space_Q.DebugStr()<<std::endl;
  }

} //main