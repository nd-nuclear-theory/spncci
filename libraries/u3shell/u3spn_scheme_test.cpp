/****************************************************************
  u3spn_scheme_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT

  9/6/16 (aem): Created.
****************************************************************/

#include "u3shell/u3spn_scheme.h"

#include "fmt/format.h"


////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{


  std::map<u3shell::U3SPN,int> subspace_dimensions;
  subspace_dimensions[u3shell::U3SPN(u3::U3S(u3::U3(16,u3::SU3(4,0)),1),0,1)] = 42;
  subspace_dimensions[u3shell::U3SPN(u3::U3S(u3::U3(11,u3::SU3(0,1)),0),1,1)] = 2;

  std::cout << "Building space"
            << std::endl;
  // build space
  u3shell::SpaceU3SPN space(subspace_dimensions);

  // dump subspace contents
      
  for (int subspace_index=0; subspace_index<space.size(); ++subspace_index)
    {
      const u3shell::SubspaceU3SPN& subspace = space.GetSubspace(subspace_index);
      
      std::cout << fmt::format("subspace {} dimension {}",subspace.LabelStr(),subspace.size()) << std::endl;
         
    }

  std::cout << "Sector enumeration -- constrained"
            << std::endl;

  // (4,0) x (x) -> (0,1); find usable operator labels by coupling (4,0) x (1,0) -> (5,0), ...
    
  u3shell::OperatorLabelsU3S op(-5,u3::SU3(0,5),1,1);  // N0, x0, S0, g0

  u3shell::SectorsU3SPN sectors(space,op,false);
      
  for (int sector_index=0; sector_index<sectors.size(); ++sector_index)
    {
      auto sector=sectors.GetSector(sector_index);

      std::cout << fmt::format(
          "{:3} ({:2},{:2},{:1}): {} {}",
          sector_index,
          sector.bra_subspace_index(),sector.ket_subspace_index(),sector.multiplicity_index(),
          sector.bra_subspace().LabelStr(),sector.ket_subspace().LabelStr()
        )
                << std::endl;
    }

  // termination
  return 0;
}
