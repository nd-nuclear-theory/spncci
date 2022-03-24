/****************************************************************
  recurrence_indexing_relative_operator_test.cpp

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  3/23/22 (aem): Created.
****************************************************************/
#include "spncci_basis/recurrence_indexing_operator.h"

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////



int main(int argc, char **argv)
{

  // u3::U3CoefInit(100);
  // test_relative();

if(true)
  {
    int Nmax=2;
    int N1v=1;
    unsigned int J0=0;
    std::vector<unsigned int> allowed_L0_values={0};

    std::unordered_set<u3::SU3> Allowed_x0_values={{0,0},{2,2}};
    std::set<unsigned int> Allowed_L0_values={0,2};
    std::set<unsigned int> Allowed_S0_values={};
    std::set<unsigned int> Allowed_T0_values={};
    relative::OperatorParameters
      operator_parameters(
        N1v,Nmax,J0,
        Allowed_x0_values,
        Allowed_L0_values,
        Allowed_S0_values,
        Allowed_T0_values
        );

    relative::spatial::OperatorSpace
      spatial_operator_space(operator_parameters);

    std::cout<<spatial_operator_space.DebugStr()<<std::endl;

    relative::spin::OperatorSpace
      spin_operator_space(operator_parameters);
    std::cout<<spin_operator_space.DebugStr()<<std::endl;

    relative::OperatorSectors operator_sectors(spatial_operator_space,spin_operator_space);
    std::cout<<operator_sectors.DebugStr()<<std::endl;
  }
  // termination
  return 0;
}
