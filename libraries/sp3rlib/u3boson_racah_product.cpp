/****************************************************************
  u3boson_racah_product.cpp

  Colin V. Coane
  
  SPDX-License-Identifier: MIT 

  08/18/23 (cvc): Created. Split from racah_product.cpp
  06/03/24 (cvc): Refactored to .h/.cpp file format
****************************************************************/

#include "sp3rlib/u3boson_racah_product.h"

namespace u3boson_operator
{

// Wrapper for Identity
basis::OperatorBlock<double> IdentityU3Boson(
  const u3::U3& omega,
  int rho,
  const u3boson::U3BosonSpace& base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  // Subspace labels
  const auto& omega_bra = bra_subspace.omega();
  const auto& omega_ket = ket_subspace.omega();

  // Return Identity matrix elements for given sector
  if (omega_bra == omega_ket)
  {
    return basis::OperatorBlock<double>::Identity(
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

// Wrapper for Creation and Annihiliation Operators
// FOR a and a^\dagger (ladder 2 U(3) quanta)!!!! 
// NOT b and b^\dagger (Weyl boson operators, 1 quanta)!!!!
// Need better notation
// Creation Operator <bra ||a^\dagger^(2,0)|| ket>_(rho)
basis::OperatorBlock<double> CreationU3Boson(
  const u3::U3& omega,
  int rho,
  const u3boson::U3BosonSpace& base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  const u3::U3 sigma = base_space.sigma();
  // Ask about doing check vs just returning regardless of omega
  u3::U3 omega_adag(2,{2u,0u});
  if (omega == omega_adag)
  {
    // Return creation operator for given sector
    return u3boson::U3BosonCreationOperator(sigma,bra_subspace,ket_subspace,u_coef_cache);
  }
  else
  {
    return basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
    );
  }
}

// Annhiliation Operator <bra ||a^(0,2)|| ket>_(rho)
basis::OperatorBlock<double> AnnihiliationU3Boson(
  const u3::U3& omega,
  int rho,
  const u3boson::U3BosonSpace& base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  const u3::U3 sigma = base_space.sigma();
  u3::U3 omega_a(-2,{0u,2u});
  // Ask about doing check vs just returning regardless of omega
  if (omega == omega_a)
  {
    // Return annihiliation operator for given sector
    return u3boson::U3BosonAnnihilationOperator(sigma,bra_subspace,ket_subspace,u_coef_cache);
  }
  else
  {
    return basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
    );
  }
}


// Wrapper for SU(3) Casimir Operator
basis::OperatorBlock<double> SU3CasimirU3Boson(
  const u3::U3& omega,
  int rho,
  const u3boson::U3BosonSpace& base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  const u3::U3 sigma = base_space.sigma();
  u3::U3 omega_casimir(0,{0u,0u});
  // Subspace labels
  const auto& omega_bra = bra_subspace.omega();
  const auto& omega_ket = ket_subspace.omega();
  
  // Return matrix elements for given sector
  if ((omega_bra == omega_ket) && (omega.SU3() == omega_casimir.SU3()))
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


// Couple U(3) operators [T^wT x S^wS]^w0
// Matrix elements of the form <bra || [T^wT x S^wS]^rho0,w0 || ket>_(outer_mult)
// 
// Function inputs currently must have the form below:
// Tfunc(wT, rhoT, base_space, bra_space, ket_space, u_coef_cache) = <bra ||T^wT|| ket>_(rhoT)
// Sfunc(wS, rhoS, base_space, bra_space, ket_space, u_coef_cache) = <bra ||S^wS|| ket>_(rhoS)
//
// Note: currently U(3)-boson<Sp(3,R) operators within a single Sp(3,R) irrep
basis::OperatorBlock<double> CoupledU3BosonOperators(
  std::function<basis::OperatorBlock<double>(
    const u3::U3&,
    int,
    const u3boson::U3BosonSpace&,
    const u3boson::U3Subspace&,
    const u3boson::U3Subspace&,
    u3::UCoefCache&
    )> Tfunc,
  std::function<basis::OperatorBlock<double>(
    const u3::U3&,
    int,
    const u3boson::U3BosonSpace&,
    const u3boson::U3Subspace&,
    const u3boson::U3Subspace&,
    u3::UCoefCache&
    )> Sfunc,
  const u3::U3& omega_T,
  const u3::U3& omega_S,
  const u3::U3& omega0,
  int rho0,
  int outer_multiplicity,
  const u3boson::U3BosonSpace& base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
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

  // Create matrix
  basis::OperatorBlock<double> TxS = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );

  //auto [Twrap, Swrap] = wrap_functions(Tfunc,Sfunc,omega_bar_vec,omega_T,omega_S,
  //                        omega_bra,base_space,bra_subspace,ket_subspace,u_coef_cache);

  // Calculate matrix elements
  // Iterate over omega = omega_ket, omega_prime = omega_bra
  // Only loop if w_ket x w0 -> w_bra AND wT x wS -> rho0 w0 is allowed
  if ((u3::OuterMultiplicity(omega_T,omega_S,omega0)>=rho0))
    if ((u3::OuterMultiplicity(omega_ket,omega0,omega_bra)>0))
    {
      for (int bra_state_index = 0; bra_state_index < bra_subspace.size(); bra_state_index++)
      {
        for (int ket_state_index = 0; ket_state_index < ket_subspace.size(); ket_state_index++)
          {
            const auto& bra_state = bra_subspace.GetState(bra_state_index);
            const auto& ket_state = ket_subspace.GetState(ket_state_index);
            const u3::U3& n_bra = bra_state.n();
            const u3::U3& n_ket = ket_state.n();
            const int rho_bra_max = bra_state.rho_max();
            const int rho_ket_max = ket_state.rho_max();

            for (int rho_bra = 1; rho_bra <= rho_bra_max; rho_bra++)
              for (int rho_ket = 1; rho_ket <= rho_ket_max; rho_ket++)
              {
                int row = bra_subspace.GetStateOffset(bra_state_index,rho_bra);
                int col = ket_subspace.GetStateOffset(ket_state_index,rho_ket);

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
    }
  // return computed sector
  return TxS;
}

// commutator [a,a^\dagger]^omega_0
// no restrictions on labels here, can use to test failure if incorrectly assigned
basis::OperatorBlock<double> CommutatorBosonQuanta(
  const u3::U3& omega_T,
  const u3::U3& omega_S,
  const u3::U3& omega0,
  int rho0,
  int outer_multiplicity,
  const u3boson::U3BosonSpace& base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
  )
{
  const u3::U3 sigma = base_space.sigma();

  // commutator matrix
  basis::OperatorBlock<double> commutator_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );

  commutator_matrix += CoupledU3BosonOperators(&AnnihiliationU3Boson,&CreationU3Boson,
                          omega_T,omega_S,omega0,rho0,outer_multiplicity,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  commutator_matrix -= CoupledU3BosonOperators(&CreationU3Boson,&AnnihiliationU3Boson,
                          omega_S,omega_T,omega0,rho0,outer_multiplicity,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  
  return commutator_matrix;
}


// commutator [a^\dagger,a^\dagger]^omega
// can be passed into coupling function CoupledU3BosonOperators()
basis::OperatorBlock<double> CommutatorCreationCreation(
  const u3::U3& omega,
  int rho,
  const u3boson::U3BosonSpace& base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
  )
{
  const u3::U3 sigma = base_space.sigma();

  u3::U3 omega_adag(2,{2u,0u});
  u3::U3 omega_a(-2,{0u,2u});
  int rho0 = 1;

  // commutator matrix
  basis::OperatorBlock<double> commutator_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );

  commutator_matrix += CoupledU3BosonOperators(&CreationU3Boson,&CreationU3Boson,
                          omega_adag,omega_adag,omega,rho0,rho,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  commutator_matrix -= CoupledU3BosonOperators(&CreationU3Boson,&CreationU3Boson,
                          omega_adag,omega_adag,omega,rho0,rho,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  
  return commutator_matrix;
}


// commutator [a,a]^omega
// can be passed into coupling function CoupledU3BosonOperators()
basis::OperatorBlock<double> CommutatorAnnihilationAnnihilation(
  const u3::U3& omega,
  int rho,
  const u3boson::U3BosonSpace& base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
  )
{
  const u3::U3 sigma = base_space.sigma();

  u3::U3 omega_adag(2,{2u,0u});
  u3::U3 omega_a(-2,{0u,2u});
  int rho0 = 1;

  // commutator matrix
  basis::OperatorBlock<double> commutator_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );

  commutator_matrix += CoupledU3BosonOperators(&AnnihiliationU3Boson,&AnnihiliationU3Boson,
                          omega_a,omega_a,omega,rho0,rho,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  commutator_matrix -= CoupledU3BosonOperators(&AnnihiliationU3Boson,&AnnihiliationU3Boson,
                          omega_a,omega_a,omega,rho0,rho,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  
  return commutator_matrix;
}


// commutator [a,a^\dagger]^omega
// can be passed into coupling function CoupledU3BosonOperators()
basis::OperatorBlock<double> CommutatorAnnihilationCreation(
  const u3::U3& omega,
  int rho,
  const u3boson::U3BosonSpace& base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
  )
{
  const u3::U3 sigma = base_space.sigma();

  u3::U3 omega_adag(2,{2u,0u});
  u3::U3 omega_a(-2,{0u,2u});
  int rho0 = 1;

  // commutator matrix
  basis::OperatorBlock<double> commutator_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );

  commutator_matrix += CoupledU3BosonOperators(&AnnihiliationU3Boson,&CreationU3Boson,
                          omega_a,omega_adag,omega,rho0,rho,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  commutator_matrix -= CoupledU3BosonOperators(&CreationU3Boson,&AnnihiliationU3Boson,
                          omega_adag,omega_a,omega,rho0,rho,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  
  return commutator_matrix;
}


/////////////////////////////////////////////
/////////////////////////////////////////////
/////////////////  TESTING  /////////////////
/////////////////////////////////////////////
/////////////////////////////////////////////

// number operator [a^\dagger x a]^omega
basis::OperatorBlock<double> CreationAnnihilation(
  const u3::U3& omega,
  int rho,
  const u3boson::U3BosonSpace& base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
  )
{
  const u3::U3 sigma = base_space.sigma();

  u3::U3 omega_adag(2,{2u,0u});
  u3::U3 omega_a(-2,{0u,2u});
  int rho0 = 1;

  // commutator matrix
  basis::OperatorBlock<double> adaga_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );

  adaga_matrix += CoupledU3BosonOperators(&CreationU3Boson,&AnnihiliationU3Boson,
                          omega_adag,omega_a,omega,rho0,rho,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  
  return adaga_matrix;
}


// commutator [[a^\dagger,a]^0(0,0), a^\dagger]^omega
basis::OperatorBlock<double> CommutatorNumberCreation(
  const u3::U3& omega,
  int rho,
  const u3boson::U3BosonSpace& base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
  )
{
  const u3::U3 sigma = base_space.sigma();

  u3::U3 omega0(0,{0u,0u});
  u3::U3 omega_adag(2,{2u,0u});
  u3::U3 omega_a(-2,{0u,2u});
  int rho0 = 1;
  int outer_multiplicity = 1;

  // Subspace labels
  const auto& omega_bra = bra_subspace.omega();
  const auto& omega_ket = ket_subspace.omega();

  // commutator matrix
  basis::OperatorBlock<double> commutator_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );

  commutator_matrix += CoupledU3BosonOperators(&CreationAnnihilation,&CreationU3Boson,
                          omega0,omega_adag,omega,rho0,rho,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  commutator_matrix -= CoupledU3BosonOperators(&CreationU3Boson,&CreationAnnihilation,
              omega_adag,omega0,omega,rho0,rho,
              base_space,bra_subspace,ket_subspace,u_coef_cache);
  
  return commutator_matrix;
}


// commutator [[a^\dagger,a]^0(0,0), a]^omega
basis::OperatorBlock<double> CommutatorNumberAnnihilation(
  const u3::U3& omega,
  int rho,
  const u3boson::U3BosonSpace& base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
  )
{
  const u3::U3 sigma = base_space.sigma();

  u3::U3 omega0(0,{0u,0u});
  u3::U3 omega_adag(2,{2u,0u});
  u3::U3 omega_a(-2,{0u,2u});
  int rho0 = 1;
  int outer_multiplicity = 1;

  // Subspace labels
  const auto& omega_bra = bra_subspace.omega();
  const auto& omega_ket = ket_subspace.omega();

  // commutator matrix
  basis::OperatorBlock<double> commutator_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );

  commutator_matrix += CoupledU3BosonOperators(&CreationAnnihilation,&AnnihiliationU3Boson,
                          omega0,omega_a,omega,rho0,rho,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  commutator_matrix -= CoupledU3BosonOperators(&AnnihiliationU3Boson,&CreationAnnihilation,
              omega_a,omega0,omega,rho0,rho,
              base_space,bra_subspace,ket_subspace,u_coef_cache);
  
  return commutator_matrix;
}


// Generators
basis::OperatorBlock<double> SU3GeneratorC(
  const u3::U3& omega,
  int rho,
  const u3boson::U3BosonSpace& base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
  )
{
  const u3::U3 sigma = base_space.sigma();
  u3::U3 omega_T(2,{2u,0u});
  u3::U3 omega_S(-2,{0u,2u});
  //u3::U3 omega0(0,{0u,0u});
  // check these labels
  u3::U3 omega0_Cas(0,{1u,1u});
  int rho0 = 1;
  int outer_multiplicity = 1;

  // Subspace labels
  const auto& omega_bra = bra_subspace.omega();
  const auto& omega_ket = ket_subspace.omega();

  // commutator matrix
  basis::OperatorBlock<double> generator_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );

  generator_matrix += std::sqrt(2)*CoupledU3BosonOperators(&CreationU3Boson,&AnnihiliationU3Boson,
                          omega_T,omega_S,omega0_Cas,rho0,outer_multiplicity,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  
  return generator_matrix;
}

// SU3 Casimir
basis::OperatorBlock<double> SU3CasimirCoupling(
  const u3::U3& omega,
  int rho,
  const u3boson::U3BosonSpace& base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
  )
{
  const u3::U3 sigma = base_space.sigma();

  u3::U3 omega_Cas(0,{1u,1u});
  u3::U3 omega0(0,{0u,0u});
  int rho0 = 1;
  int outer_multiplicity = 1;

  // Subspace labels
  const auto& omega_bra = bra_subspace.omega();
  const auto& omega_ket = ket_subspace.omega();

  // commutator matrix
  basis::OperatorBlock<double> casimir_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );

  casimir_matrix += (0.5)*CoupledU3BosonOperators(&SU3GeneratorC,&SU3GeneratorC,
                          omega_Cas,omega_Cas,omega0,rho0,outer_multiplicity,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  
  return casimir_matrix;
}

}