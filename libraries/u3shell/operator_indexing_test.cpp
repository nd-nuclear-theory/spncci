/****************************************************************
  recurrence_indexing_relative_operator_test.cpp

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  3/23/22 (aem): Created.
****************************************************************/
#include "u3shell/operator_indexing_spatial.h"
#include "u3shell/operator_indexing_spin.h"
#include "u3shell/operator_indexing_sectors.h"

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////



int main(int argc, char **argv)
{
  using relative_operator_sectors
  = u3shell::relative::OperatorU3SpinSectors<
    u3shell::spatial::onecoord::OperatorSpace,
    u3shell::spatial::onecoord::OperatorL0Space,
    u3shell::spin::twobody::OperatorSpace,
    u3shell::spin::twobody::OperatorSubspace
    >;


  if(true)
  {
    // Check indexing for Identity
    int Nmax=4;
    int N1v=1;

    // Parameters are N1v, Nmax, J0, Allowed_w0_values,
    //  Allowed_L0_values, Allowed_S0_values, Allowed T0_values
    u3shell::relative::OperatorParameters Identity_parameters(N1v,Nmax,0,{{0,{0u,0u}}},{0},{0},{0});

    std::cout<<"Identity operator"<<std::endl;

    auto spatial_ptr
    = std::make_shared<const u3shell::spatial::onecoord::OperatorSpace>(Identity_parameters);
    std::cout<<spatial_ptr->DebugStr()<<std::endl;

    auto spin_ptr
    = std::make_shared<const u3shell::spin::twobody::OperatorSpace>(Identity_parameters);
    std::cout<<spin_ptr->DebugStr()<<std::endl;

    relative_operator_sectors operator_sectors(spatial_ptr,spin_ptr,Identity_parameters.J0);
    std::cout<<operator_sectors.DebugStr()<<std::endl;

  }

  if(false)
  {
    // Check indexing for Sp(3,R) raising operator
    int Nmax=4;
    int N1v=1;

    // Parameters are N1v, Nmax, J0, Allowed_w0_values,
    //  Allowed_L0_values, Allowed_S0_values, Allowed T0_values
    u3shell::relative::OperatorParameters Arel_parameters(N1v,Nmax,0,{{2,{2u,0u}}},{0,2},{0},{0});

    std::cout<<"Relative Sp(3,R) raising operator"<<std::endl;

    auto spatial_ptr
    = std::make_shared<const u3shell::spatial::onecoord::OperatorSpace>(Arel_parameters);
    std::cout<<spatial_ptr->DebugStr()<<std::endl;

    auto spin_ptr
    = std::make_shared<const u3shell::spin::twobody::OperatorSpace>(Arel_parameters);
    std::cout<<spin_ptr->DebugStr()<<std::endl;

    relative_operator_sectors operator_sectors(spatial_ptr,spin_ptr,Arel_parameters.J0);
    std::cout<<operator_sectors.DebugStr()<<std::endl;

  }


  if(false)
  {

    // Check indexing for Trel
    int Nmax=2;
    int N1v=1;
    u3shell::relative::OperatorParameters
      Trel_parameters(N1v,Nmax,0,{{2,{2u,0u}},{-2,{0u,2u}},{0,{0u,0u}}},{0},{0},{0});

    std::cout<<"Relative kinetic energy"<<std::endl;
    auto spatial_ptr
    = std::make_shared<const u3shell::spatial::onecoord::OperatorSpace>(Trel_parameters);
    std::cout<<spatial_ptr->DebugStr()<<std::endl;

    auto spin_ptr
    = std::make_shared<const u3shell::spin::twobody::OperatorSpace>(Trel_parameters);
    std::cout<<spin_ptr->DebugStr()<<std::endl;

    relative_operator_sectors operator_sectors(spatial_ptr,spin_ptr,Trel_parameters.J0);
    std::cout<<operator_sectors.DebugStr()<<std::endl;

  }


  if(false)
  {
    // No contraints other than J0.
    std::cout<<"No constraints"<<std::endl;
    int Nmax=20;
    int N1v=1;
    unsigned int J0=0;
    std::unordered_set<u3::U3> Allowed_w0_values;
    std::set<unsigned int> Allowed_L0_values;
    std::set<uint8_t> Allowed_S0_values,Allowed_T0_values;
    u3shell::relative::OperatorParameters
      hamiltonian_parameters(
        N1v,Nmax,J0,
        Allowed_w0_values,
        Allowed_L0_values,
        Allowed_S0_values,
        Allowed_T0_values
        );

    auto spatial_ptr
    = std::make_shared<const u3shell::spatial::onecoord::OperatorSpace>(hamiltonian_parameters);
    // std::cout<<spatial_ptr->DebugStr()<<std::endl;

    auto spin_ptr
    = std::make_shared<const u3shell::spin::twobody::OperatorSpace>(hamiltonian_parameters);
    // std::cout<<spin_ptr->DebugStr()<<std::endl;

    relative_operator_sectors operator_sectors(spatial_ptr,spin_ptr,hamiltonian_parameters.J0);
    std::cout<<operator_sectors.DebugStr()<<std::endl;


  }
  return 0;
}
