/****************************************************************
  racah_product.cpp

  Colin V. Coane
  
  SPDX-License-Identifier: MIT 

  04/02/22 (cvc): Created.

****************************************************************/

#include "sp3rlib/sp3r_operator.h"

#include <iostream>
#include <functional>

#include <fstream>
#include <iterator>
#include <map>
#include <tuple>
#include <vector>
#include <utility>

#include "am/halfint.h"
#include "am/wigner_gsl.h"
#include "fmt/format.h"

#include "sp3rlib/u3coef.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3boson.h"
#include "fmt/format.h"
#include "mcutils/eigen.h"
#include "sp3rlib/u3.h"


// Ask about what to do with wrappers
// Wrapper for Tfunc and Sfunc, returns blocks
// given up a bit, need to fix and templatize!!!!
//template<typename F1, typename F2>
auto wrap_functions(
  //F1&& f1, 
  //F2&& f2, 
  std::function<basis::OperatorBlock<double>(
  const u3::U3&,
  int,
  const u3::U3&,
  const u3boson::U3Subspace&,
  const u3boson::U3Subspace&,
  u3::UCoefCache
  )> f1,
std::function<basis::OperatorBlock<double>(
  const u3::U3&,
  int,
  const u3::U3&,
  const u3boson::U3Subspace&,
  const u3boson::U3Subspace&,
  u3::UCoefCache
  )> f2,
  MultiplicityTagged<u3::U3>::vector omega_bar_vec, 
  const u3::U3 omega_T,
  const u3::U3 omega_S,
  const u3::U3 omega_bra,
  const u3boson::U3BosonSpace base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache u_coef_cache
  )
{
  //std::cout<<fmt::format("Initialization of wrapper")<<std::endl;
  const u3::U3 sigma = base_space.sigma();
  std::vector<std::vector<basis::OperatorBlock<double>>> Tmat;
  std::vector<std::vector<basis::OperatorBlock<double>>> Smat;
  // Iterate over rho_T, rho_S, omega_bar first, get max of max values
  //std::cout<<fmt::format("Iterate over omega_bar, rho_S_max")<<std::endl;
  for (auto& [omega_bar, rho_S_max]: omega_bar_vec)
  {
    std::cout<<omega_bar.Str()<<std::endl;
    std::cout<<rho_S_max<<std::endl;
  }
  for (auto& [omega_bar, rho_S_max]: omega_bar_vec)
  {
    // Get rho_T for summing from w_bar x wT -> wprime
    int rho_T_max = u3::OuterMultiplicity(omega_bar,omega_T,omega_bra);
    // get subspace labels for omega_bar
    //std::cout<<fmt::format("get barred subspace")<<std::endl;
    //std::cout<<base_space.DebugStr()<<std::endl;
    //std::cout<<omega_bar.Str()<<std::endl;
    //if (base_space.ContainsSubspace({omega_bar}))
    u3boson::U3Subspace barred_subspace = base_space.LookUpSubspace({omega_bar});
    // temporary vectors
    std::vector<basis::OperatorBlock<double>> T_temp;
    std::vector<basis::OperatorBlock<double>> S_temp;
    //std::cout<<fmt::format("Iterate rho_T")<<std::endl;
    for (int rho_T = 1; rho_T <= rho_T_max; rho_T++)
    {
      // Get T matrix
      //basis::OperatorBlock<double> f1mat = f1(omega_T,rho_T,sigma,bra_subspace,barred_subspace,u_coef_cache);
      //T_temp.push_back(f1mat);
      T_temp.push_back(f1(omega_T,rho_T,sigma,bra_subspace,barred_subspace,u_coef_cache));
    }
    //std::cout<<fmt::format("Iterate rho_S")<<std::endl;
    for (int rho_S = 1; rho_S <= rho_S_max; rho_S++)
    {
      // Get S matrix
      //basis::OperatorBlock<double> f2mat = f2(omega_S,rho_S,sigma,barred_subspace,ket_subspace,u_coef_cache);
      //S_temp.push_back(f2mat);
      S_temp.push_back(f2(omega_S,rho_S,sigma,barred_subspace,ket_subspace,u_coef_cache));
    }
    //std::cout<<fmt::format("Push back vectors")<<std::endl;
    Tmat.push_back(T_temp);
    Smat.push_back(S_temp);
  }
  //std::cout<<fmt::format("Return functions")<<std::endl;
  return std::make_pair(
    [Twrap = std::move(Tmat)](int u, int rho_T) {return Twrap[u][rho_T-1]; },
    [Swrap = std::move(Smat)](int v, int rho_S) {return Swrap[v][rho_S-1]; }
  );
}

// Wrapper for Identity
basis::OperatorBlock<double> IdentityU3Boson(
  const u3::U3 omega,
  int rho,
  const u3::U3& sigma,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache u_coef_cache
)
{
  // Subspace labels
  const auto& omega_bra = bra_subspace.omega();
  const auto& omega_ket = ket_subspace.omega();

  // Return Identity elements for given block
  if (omega_bra == omega_ket)
  {
    return basis::OperatorBlock<double>::Identity(
      bra_subspace.dimension(),
      ket_subspace.dimension()
    );
  }
  else
  {
    return  basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
    );
  }
}

// Wrapper for Creation and Annihiliation Operators
// Creation Operator
basis::OperatorBlock<double> CreationU3Boson(
  const u3::U3 omega,
  int rho,
  const u3::U3& sigma,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache u_coef_cache
)
{
  // Subspace labels
  const auto& omega_bra = bra_subspace.omega();
  const auto& omega_ket = ket_subspace.omega();

  // Return creation operator for given sectors
  return u3boson::U3BosonCreationOperator(sigma,bra_subspace,ket_subspace,u_coef_cache);
}

// Annhiliation Operator
basis::OperatorBlock<double> AnnihiliationU3Boson(
  const u3::U3 omega,
  int rho,
  const u3::U3& sigma,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache u_coef_cache
)
{
  // Subspace labels
  const auto& omega_bra = bra_subspace.omega();
  const auto& omega_ket = ket_subspace.omega();

  // Return annihiliation operator for given sectors
  return u3boson::U3BosonAnnihilationOperator(sigma,bra_subspace,ket_subspace,u_coef_cache);
}



// Couple U(3) operators
// TODO write example functions for Tfunc and Sfunc
//
// Function inputs currently must have the form below:
// Tfunc(wT, rhoT, sigma, bra_space, ket_space) = <bra ||T^wT|| ket>_(rhoT)
// Sfunc(wS, rhoS, sigma, bra_space, ket_space) = <bra ||S^wS|| ket>_(rhoS)
//
basis::OperatorBlock<double> CoupledU3Operators(
  std::function<basis::OperatorBlock<double>(
    const u3::U3&,
    int,
    const u3::U3&,
    const u3boson::U3Subspace&,
    const u3boson::U3Subspace&,
    u3::UCoefCache
    )> Tfunc,
  std::function<basis::OperatorBlock<double>(
    const u3::U3&,
    int,
    const u3::U3&,
    const u3boson::U3Subspace&,
    const u3boson::U3Subspace&,
    u3::UCoefCache
    )> Sfunc,
  const u3::U3& omega_T,
  const u3::U3& omega_S,
  const u3::U3& omega0,
  int rho0,
  int outer_multiplicity,
  const u3boson::U3BosonSpace base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache u_coef_cache
  )
{
  //std::cout<<fmt::format("Initializations")<<std::endl;
  // Get sigma
  const u3::U3 sigma = base_space.sigma();

  // Subspace labels
  const auto& omega_bra = bra_subspace.omega();
  const auto& omega_ket = ket_subspace.omega();

  // Get labels to sum over
  // rho0 prime for interchanging wS x wT -> w0
  int rho0prime_max = u3::OuterMultiplicity(omega_S,omega_T,omega0);
  // omega_ bar from w x wS -> w_bar
  MultiplicityTagged<u3::U3>::vector omega_bar_vec = u3::KroneckerProduct(omega_ket,omega_S);
  //for (auto& [omega_bar, rho_S_max]: omega_bar_vec) { std::cout<<omega_bar.Str()<<std::endl; }
  std::vector<int> removal_idx;
  //for (int idx_bar = 0; idx_bar<omega_bar_vec.size();idx_bar++) // doesn't seem to iterate at all
  for (auto& [omega_bar, rho_S_max]: omega_bar_vec)
  {
    //auto& [omega_bar, rho_S_max] = omega_bar_vec[idx_bar];
    if (base_space.LookUpSubspaceIndex({omega_bar})>base_space.size())
    {
      std::cout<<fmt::format("\n{} not in base_space",omega_bar.Str())<<std::endl;
      //removal_idx.push_back(idx_bar);
      //omega_bar_vec.erase(std::remove(omega_bar_vec.begin(),omega_bar_vec.end(),[omega_bar, rho_S_max]),omega_bar_vec.end());
      // need to remove omega_bar that don't belong to base_space
    }
  }



  // Matrix elements of coupled operators
  //std::cout<<fmt::format("Create TxS operator")<<std::endl;
  basis::OperatorBlock<double> TxS = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );

  /*
  std::vector<std::vector<basis::OperatorBlock<double>>> Tmat;
  std::vector<std::vector<basis::OperatorBlock<double>>> Smat;
  // Iterate over rho_T, rho_S, omega_bar first, get max of max values
  for (auto& [omega_bar, rho_S_max]: omega_bar_vec)
  {
    // Get rho_T for summing from w_bar x wT -> wprime
    int rho_T_max = u3::OuterMultiplicity(omega_bar,omega_T,omega_bra);
    // get subspace labels for omega_bar
    u3boson::U3Subspace barred_subspace = base_space.LookUpSubspace({omega_bar});
    // temporary vectors
    std::vector<basis::OperatorBlock<double>> T_temp;
    std::vector<basis::OperatorBlock<double>> S_temp;
    for (int rho_T = 1; rho_T <= rho_T_max; rho_T++)
    {
      // Get T matrix
      T_temp.push_back(Tfunc(omega_T,rho_T,sigma,bra_subspace,barred_subspace));
    }
    for (int rho_S = 1; rho_S <= rho_S_max; rho_S++)
    {
      // Get S matrix
      S_temp.push_back(Sfunc(omega_S,rho_S,sigma,barred_subspace,ket_subspace));
    }
    Tmat.push_back(T_temp);
    Smat.push_back(S_temp);
  }
  */

  // Wrapper to look up Tmat and Smat blocks
  //auto [Twrap, Swrap] = wrap_functions(&Tfunc, &Sfunc,omega_bar_vec,omega_T,omega_S,
  //                        omega_bra,base_space,bra_subspace,ket_subspace,u_coef_cache);
  //std::cout<<fmt::format("Wrap functions")<<std::endl;
  auto [Twrap, Swrap] = wrap_functions(Tfunc,Sfunc,omega_bar_vec,omega_T,omega_S,
                          omega_bra,base_space,bra_subspace,ket_subspace,u_coef_cache);

  // Calculate matrix elements of coupled product of operators
  // Iterate over omega = omega_ket, omega_prime = omega_bra
  //std::cout<<fmt::format("Begin loop")<<std::endl;
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

              // TODO move rho0prime sum inside rho_S,rho_T sums to speed up computation
              for (int rho0prime = 1; rho0prime <= rho0prime_max; rho0prime++)
              {
                for (int u = 0; u < omega_bar_vec.size(); u++)
                //for (auto& [omega_bar, rho_S_max]: omega_bar_vec)
                {
                  // extract omega_bar
                  u3::U3 omega_bar = omega_bar_vec[u].irrep;
                  int rho_S_max = omega_bar_vec[u].tag;
                  // get subspace labels for omega_bar
                  u3boson::U3Subspace barred_subspace = base_space.LookUpSubspace({omega_bar});
                  // Get rho_T for summing from w_bar x wT -> wprime
                  int rho_T_max = u3::OuterMultiplicity(omega_bar,omega_T,omega_bra);
                  for (int rho_S = 1; rho_S <= rho_S_max; rho_S++)
                  {
                    for (int rho_T = 1; rho_T <= rho_T_max; rho_T++)
                    {
                      // call Tfunc, Sfunc here
                      // extract and sum matrix elements of uncoupled operators
                      //std::cout<<fmt::format("Calling Tfunc")<<std::endl;
                      basis::OperatorBlock<double> Tmat = Tfunc(omega_T,rho_T,sigma,bra_subspace,barred_subspace,u_coef_cache);
                      //std::cout<<fmt::format("Calling Sfunc")<<std::endl;
                      basis::OperatorBlock<double> Smat = Sfunc(omega_S,rho_S,sigma,barred_subspace,ket_subspace,u_coef_cache);
                      //basis::OperatorBlock<double> TS_dot = Twrap.row(row)*Smat.col(col);
                      //std::cout<<fmt::format("T dot S")<<std::endl;
                      basis::OperatorBlock<double> TS_dot = Twrap(u,rho_T).row(row)*Swrap(u,rho_S).col(col);
                      //basis::OperatorBlock<double> TS_dot = Tmat[u][rho_T-1].row(row)*Smat[u][rho_S-1].col(col);
                      // matrix product
                      // note the u index takes care of extracting omega_bar sectors
                      // sum ends up only being over nbar, rhobar, for a given omega_bar
                      // sum[<n' rho' w'||T^wT||nbar rhobar wbar> <nbar rhobar wbar||S^wS||n rho w>]
                      // index and extract correct product of matrix elements
                      // given known omegaprime, omegabar, omega
                      // omega bar taken care of by taking dot product of sectors
                      // compute sum
                      //std::cout<<fmt::format("Adding to sum")<<std::endl;
                      sum += u3::Phi(omega_T.SU3(),omega_S.SU3(),omega0.SU3(),rho0,rho0prime)
                        * u3::U(omega_ket.SU3(),omega_S.SU3(),omega_bra.SU3(),
                        omega_T.SU3(),omega_bar.SU3(),rho_S,rho_T,
                        omega0.SU3(),rho0prime,outer_multiplicity)
                        * (TS_dot(0,0));
                    }
                  }
                }
              }
              //std::cout<<fmt::format("Add element to TxS")<<std::endl;
              TxS(row,col) = sum;
            }
      }
  }

  return TxS;
}

// commutator [b^dagger,b]^omega_0
basis::OperatorBlock<double> CommutatorBB(
  const u3::U3& omega_T,
  const u3::U3& omega_S,
  const u3::U3& omega0,
  int rho0,
  int outer_multiplicity,
  const u3boson::U3BosonSpace base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache u_coef_cache
  )
{
  // Get sigma
  const u3::U3 sigma = base_space.sigma();

  // Subspace labels
  const auto& omega_bra = bra_subspace.omega();
  const auto& omega_ket = ket_subspace.omega();

  // commutator matrix
  basis::OperatorBlock<double> commutator_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );

  commutator_matrix += CoupledU3Operators(&CreationU3Boson,&AnnihiliationU3Boson,
                          omega_T,omega_S,omega0,rho0,outer_multiplicity,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  commutator_matrix -= CoupledU3Operators(&AnnihiliationU3Boson,&CreationU3Boson,
              omega_S,omega_T,omega0,rho0,outer_multiplicity,
              base_space,bra_subspace,ket_subspace,u_coef_cache);
  
  return commutator_matrix;
}




// Test coupling identity to itself
int main(int argc, char **argv)
{
  //std::cout<<fmt::format("Beginning initialization...")<<std::endl;
  //std::cout<<fmt::format("U3CoefInit...")<<std::endl;
  u3::U3CoefInit(100); // what size should max_lambda_plus_mu be???
  //std::cout<<fmt::format("U3 sigma...")<<std::endl;
  
  u3::U3 sigma(0,{0,0});
  int Nn_max=6;
  //std::cout<<fmt::format("U3BosonSpace...")<<std::endl;
  u3boson::U3BosonSpace u3boson_space(sigma,Nn_max);
  //std::cout<<fmt::format("U3CoefCache...")<<std::endl;
  u3::UCoefCache u_coef_cache;
  int rho0 = 1;
  int outer_multiplicity = 1;
  //std::cout<<fmt::format("U3's...")<<std::endl;
  u3::U3 omega_T(0,{0,0});
  u3::U3 omega_S(0,{0,0});
  u3::U3 omega0(0,{0,0});
  // check these labels
  u3::U3 omega_T_b(2,{2,0});
  u3::U3 omega_S_b(-2,{0,2});
  u3::U3 omega0_b(0,{0,0});
  //std::cout<<fmt::format("U3BosonSpace:")<<std::endl;
  std::cout<<u3boson_space.DebugStr()<<std::endl;
  //std::cout<<u3boson_space.size()<<std::endl;
  //for (int idx = 0; idx<u3boson_space.size();idx++) // doesn't seem to iterate at all
  for (const auto& bra_subspace : u3boson_space)
  {
    //auto& bra_subspace = u3boson_space.GetSubspace(idx);
    //std::cout<<bra_subspace.DebugStr()<<std::endl;
    for (const auto& ket_subspace : u3boson_space)
    {
      //std::cout<<bra_subspace.DebugStr()<<std::endl;
      //std::cout<<ket_subspace.DebugStr()<<std::endl;
      // compute rme of Identity coupled to itself
      std::cout<<fmt::format("Testing I x I")<<std::endl;
      basis::OperatorBlock<double> IxI = CoupledU3Operators(
        &IdentityU3Boson,
        &IdentityU3Boson,
        omega_T,
        omega_S,
        omega0,
        rho0,
        outer_multiplicity,
        u3boson_space,
        bra_subspace,
        ket_subspace,
        u_coef_cache
      );
      //mcutils::IsZero(operator_block-validation_matrix,1e-8) // use this to check
      if (IxI != IdentityU3Boson(omega0,rho0,sigma,bra_subspace,ket_subspace,u_coef_cache))
      {
        std::cout<<fmt::format("I x I incorrect rme")<<std::endl;
        //std::cout<<IxI<<std::endl;
        //std::cout<<IdentityU3Boson(omega0,rho0,sigma,bra_subspace,ket_subspace,u_coef_cache)<<std::endl;
      }
      else
      {
        //std::cout<<fmt::format("I x I correct rme")<<std::endl;
        //std::cout<<IxI<<std::endl;
        //std::cout<<IdentityU3Boson(omega0,rho0,sigma,bra_subspace,ket_subspace,u_coef_cache)<<std::endl;
      }
      // compute rme of [b^dagger,b], check if proportional to identity
      std::cout<<fmt::format("Testing [b^dagger,b]")<<std::endl;
      basis::OperatorBlock<double> bDaggerb = CommutatorBB(
        omega_T_b,
        omega_S_b,
        omega0_b,
        rho0,
        outer_multiplicity,
        u3boson_space,
        bra_subspace,
        ket_subspace,
        u_coef_cache
      );
      std::cout<<bDaggerb(0,0)<<std::endl;
      if (bDaggerb != (-1.0)*std::sqrt(6)*IdentityU3Boson(
            omega0,rho0,sigma,bra_subspace,ket_subspace,u_coef_cache))
      {
        std::cout<<fmt::format("[b^dagger,b] incorrect rme")<<std::endl;
        //std::cout<<bDaggerb<<std::endl;
        //std::cout<<bDaggerb(0,0)*IdentityU3Boson(
        //    omega0,rho0,sigma,bra_subspace,ket_subspace,u_coef_cache)<<std::endl;
      }
      else
      {
        //std::cout<<fmt::format("[b^dagger,b] correct rme")<<std::endl;
        //std::cout<<bDaggerb<<std::endl;
        //std::cout<<bDaggerb(0,0)*IdentityU3Boson(
        //    omega0,rho0,sigma,bra_subspace,ket_subspace,u_coef_cache)<<std::endl;
      }
    }
  }
}



















/*
#include <vector>
#include <utility>
#include <fmt/format.h>

template<typename F1, typename F2>
auto foo(F1&& f1, F2&& f2, int imax)
{
    std::vector<int> f1_vec;
    std::vector<int> f2_vec;

    for (int i = 0; i < imax; ++i)
    {
        f1_vec.push_back(f1(i));
        f2_vec.push_back(f2(i));
    }

    return std::make_pair(
        [f1_vec_ = std::move(f1_vec)](int i) { return f1_vec_[i]; }, 
        [f2_vec_ = std::move(f2_vec)](int i) { return f2_vec_[i]; }
        );
}

int bar1(int i)
{
    fmt::print("in bar1({})\n", i);
    return 2*i;
}

int bar2(int i)
{
    fmt::print("in bar2({})\n", i);
    return i*i;
}

int main()
{
    auto [bar1_wrapped, bar2_wrapped] = foo(&bar1, &bar2, 10);
    fmt::print("{}\n", typeid(bar1_wrapped).name());

    for (int i = 0; i<10; ++i)
    {
        fmt::print("{} {}\n", bar1_wrapped(i), bar2_wrapped(i));
    }
    return 0;
}
*/






/*
// SO3Operator class
class SO3Operator
{
  public:
  constexpr SO3Operator(HalfInt L0_, 
          std::function<double(std::size_t,HalfInt,std::size_t,Halfint)> f,
          std::function<std::pair<std::size_t,std::size_t>(HalfInt,HalfInt)> bounds_f)
            : L0{std::move(L0_)}, func{std::move(f)}, bounds_func{std::move(bounds_f)} {}
  
  double inline operator() (
    std::size_t gamma_bra, HalfInt L_bra, std::size_t gamma_ket, HalfInt L_ket
    )
  {
    return func(gamma_bra, L_bra, gamma_ket, L_ket);
  }

  std::pair<std::size_t,std::size_t> inline gamma_max(HalfInt L_bra, HalfInt L_ket)
  {
    return bounds_func(L_bra, L_ket);
  }

  private:
  HalfInt L0;
  std::function<std::pair<std::size_t,std::size_t>(HalfInt,HalfInt)> bounds_func;
  std::function<double(std::size_t,HalfInt,std::size_t,HalfInt)> func;
};



// Coupling of SO(3) Tensors
basis::OperatorBlock<double> CoupledSO3Operators(
  SO3Operator& T,
  SO3Operator& S,
  HalfInt L0,
  const u3::U3& sigma,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  )
{
  // Angular momentum
  HalfInt LT = T.L0;
  HalfInt LS = S.L0;

  // Upstream labels
  const auto& omega_bra = bra_subspace.omega();
  const auto& omega_ket = ket_subspace.omega();
  if (!u3::OuterMultiplicity(omega_ket,{2,{2,0}},omega_bra))
      return basis::OperatorBlock<double>::Zero(bra_subspace.upsilon_max(),ket_subspace.upsilon_max());

  // Matrix elements of coupled operators
  basis::OperatorBlock<double> TxS
    = basis::OperatorBlock<double>::Zero(
      bra_subspace.size(),
      ket_subspace.size()
      );

  // Calculate matrix elements for each state
  for (int bra_state_index=0; bra_state_index<bra_subspace.size(); ++bra_state_index)
    for (int ket_state_index=0; ket_state_index<ket_subspace.size(); ++ket_state_index)
      {
        const auto& bra_state = bra_subspace.GetState(bra_state_index);
        const auto& ket_state = ket_subspace.GetState(ket_state_index);
        const int rho_bra_max = bra_state.rho_max();
        const int rho_ket_max = ket_state.rho_max();
        const int L_bra = bra_state.L;
        const int L_ket = ket_state.L;

        for (int rho_bra=1; rho_bra<=rho_bra_max; ++rho_bra)
          for( int rho_ket=1; rho_ket<=rho_ket_max; ++rho_ket)
            {
              // State offset in each subspace
              int row = bra_subspace.GetStateOffset(bra_state_index,rho_bra);
              int col = ket_subspace.GetStateOffset(ket_state_index,rho_ket);
              
              // Sum over matrix elements of uncoupled operators
              double sum = 0;

              for (
                // sum over | (l_bar,m_bar) L_bar > < (l_bar,m_bar) L_bar | 
                )
              {
                HalfInt Lbar; // need to compute allowed values
                // sum over barred terms
                sum += pow(-1, LT+LS-L0)*pow(-1, L_ket + LS + LT + L_bra)
                  *Hat(Lbar)*Hat(L0)*am::Wigner6J(L_ket,LS,Lbar,LT,L_bra,L0)
                  *T(bra_state, ket_bar)*S(bra_bar, ket_state);
              }

              // Matrix element of coupled operators
              TxS(row,col) = sum;
            }
      }

  return TxS;
}





void main()
{
    basis::OperatorBlocks<double> Q, L, QxL;
    u3shell::Sectors Q_indexing, L_indexing, QxL_sectors;

    {Q,Q_indexing} = generateQ();
    {L,L_indexing} = generateL();

    SO3Operator Qop{{1,1}, [Q, Q_indexing](std::size_t gamma_bra, HalfInt x_bra, std::size_t gamma_ket, HalfInt x_ket) { 
        // std::size_t bra_subspace_index = Q_indexing.GetSubspace(x_bra);
        // std::size_t ket_subspace_index = Q_

        // return Q[sector_index](gamma_bra, gamma_ket);
        }};
    
    for (auto&& [sector_index,sector] : enumerate(sectors))
        for (std::size_t i=0; i<sector.rows(); ++i)
            for (std::size_t j=0; j<sector.cols(); ++j)
                QxL[sector_index](i,j) = RacahProduct({1,1}, Qop, Lop);
}
*/
