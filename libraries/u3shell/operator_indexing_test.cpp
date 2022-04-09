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

  if(true)
  {
    // Check indexing for Sp(3,R) raising operator
    int Nmax=2;
    int N1v=1;
    relative::OperatorParameters Arel_parameters(N1v,Nmax,0,{{2,{2,0}}},{0,2},{0},{0});

    std::cout<<"Relative Sp(3,R) raising operator"<<std::endl;
    relative::spatial::OperatorSpace
      spatial_operator_space(Arel_parameters);

    std::cout<<spatial_operator_space.DebugStr()<<std::endl;

    relative::spin::OperatorSpace
      spin_operator_space(Arel_parameters);
    std::cout<<spin_operator_space.DebugStr()<<std::endl;

    relative::OperatorSectors operator_sectors(spatial_operator_space,spin_operator_space);
    std::cout<<operator_sectors.DebugStr()<<std::endl;

  }


  if(true)
  {

    // Check indexing for Trel
    int Nmax=2;
    int N1v=1;
    relative::OperatorParameters Trel_parameters(N1v,Nmax,0,{{2,{2,0}},{-2,{0,2}},{0,{0,0}}},{0},{0},{0});

    std::cout<<"Relative kinetic energy"<<std::endl;
    relative::spatial::OperatorSpace
      spatial_operator_space(Trel_parameters);

    std::cout<<spatial_operator_space.DebugStr()<<std::endl;

    relative::spin::OperatorSpace
      spin_operator_space(Trel_parameters);
    std::cout<<spin_operator_space.DebugStr()<<std::endl;

    relative::OperatorSectors operator_sectors(spatial_operator_space,spin_operator_space);
    std::cout<<operator_sectors.DebugStr()<<std::endl;

  }


  if(false)
  {
    int Nmax=2;
    int N1v=1;
    unsigned int J0=0;
    std::vector<unsigned int> allowed_L0_values={0};

    std::unordered_set<u3::U3> Allowed_w0_values={{0,{0,0}},{0,{2,2}}};
    std::set<unsigned int> Allowed_L0_values={0,2};
    std::set<unsigned int> Allowed_S0_values={};
    std::set<unsigned int> Allowed_T0_values={};
    relative::OperatorParameters
      operator_parameters(
        N1v,Nmax,J0,
        Allowed_w0_values,
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

    for(int i=0; i<operator_sectors.size(); ++i)
      {
        const auto& sector = operator_sectors.GetSector(i);
        std::cout<<fmt::format("bra_degeneracy: {}  ket_degeneracy: {}",
          sector.bra_subspace_degeneracy(),
          sector.ket_subspace_degeneracy()
        )<<std::endl;

        std::cout<<fmt::format("num_elements: {}  num_block_elements: {}",
          sector.num_elements(), sector.num_block_elements()
        )<<std::endl;

        for(int bra_degeneracy_index=1; bra_degeneracy_index<=sector.bra_subspace_degeneracy(); bra_degeneracy_index++)
          for(int ket_degeneracy_index=1; ket_degeneracy_index<=sector.ket_subspace_degeneracy(); ket_degeneracy_index++)
            {
              std::cout<<fmt::format("bra_degeneracy_index: {}  ket_degeneracy_index: {}  block offset: {}",
                bra_degeneracy_index,ket_degeneracy_index,sector.GetBlockOffset(bra_degeneracy_index,ket_degeneracy_index)
              )<<std::endl;
            }
      }
  }
  // termination
  return 0;
}
