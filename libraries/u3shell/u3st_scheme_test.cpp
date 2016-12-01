/****************************************************************
  u3st_scheme_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "u3shell/u3st_scheme.h"

#include "cppformat/format.h"


////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  ////////////////////////////////////////////////////////////////
  // Relative state tests
  ////////////////////////////////////////////////////////////////

  if (true)
    {
      // build space
      int Nmax=8;
      u3shell::RelativeSpaceU3ST space(Nmax);
      std::cout << space.Str();

      // dump subspace contents
      
      for (int subspace_index=0; subspace_index<space.size(); ++subspace_index)
        {
          const u3shell::RelativeSubspaceU3ST& subspace = space.GetSubspace(subspace_index);
          std::cout << fmt::format("subspace {}",subspace.Str()) << std::endl;          
        }

    }
  ////////////////////////////////////////////////////////////////
  // two-body state tests
  ////////////////////////////////////////////////////////////////

  if (true)
    {
      // build space
      int Nmax=8;
      u3shell::TwoBodySpaceU3ST space(Nmax);
      std::cout << space.Str();

      // dump subspace contents
      
      for (int subspace_index=0; subspace_index<space.size(); ++subspace_index)
        {
          const u3shell::TwoBodySubspaceU3ST& subspace = space.GetSubspace(subspace_index);
      
          std::cout << fmt::format("subspace {}",subspace.Str()) << std::endl;

          for (int state_index=0; state_index<subspace.size(); ++state_index)
            {
          
              const u3shell::TwoBodyStateU3ST state(subspace,state_index);
              std::string state_string = fmt::format("  | {:2} {:2} >",state.N1(),state.N2());
              std::cout << state_string << std::endl;
            }
          
        }
      
    }
  
  ////////////////////////////////////////////////////////////////
  // two-body sectors
  ////////////////////////////////////////////////////////////////

  if (true)
    {
      std::cout << "Sector enumeration"
                << std::endl;

      // build space
      int Nmax=4;
      u3shell::TwoBodySpaceU3ST space(Nmax);
      std::cout << space.Str();
      
      u3shell::OperatorLabelsU3ST op(0,u3::SU3(0,0),0,0,0);

      u3shell::TwoBodySectorsU3ST sectors(space,op);
      
      for (int sector_index=0; sector_index<sectors.size(); ++sector_index)
        {
          auto sector=sectors.GetSector(sector_index);

          std::cout << fmt::format(
                                   "{:3} ({:2},{:2},{:1}): {} {}",
                                   sector_index,
                                   sector.bra_subspace_index(),sector.ket_subspace_index(),sector.multiplicity_index(),
                                   sector.bra_subspace().Str(),sector.ket_subspace().Str()
                                   )
                    << std::endl;
        }
    }


  if (true)
    {
      std::cout << "Sector enumeration -- all-to-all"
                << std::endl;

      // build space
      int Nmax=4;
      u3shell::TwoBodySpaceU3ST space(Nmax);
      std::cout << space.Str();
      
      u3shell::TwoBodySectorsU3ST sectors(space);
      
      for (int sector_index=0; sector_index<sectors.size(); ++sector_index)
        {
          auto sector=sectors.GetSector(sector_index);

          std::cout << fmt::format(
                                   "{:3} ({:2},{:2},{:1}): {} {}",
                                   sector_index,
                                   sector.bra_subspace_index(),sector.ket_subspace_index(),sector.multiplicity_index(),
                                   sector.bra_subspace().Str(),sector.ket_subspace().Str()
                                   )
                    << std::endl;
        }
    }

  if (true)
    {
      std::cout << "Sector enumeration -- constrained"
                << std::endl;

      // build space
      int Nmax=4;
      u3shell::TwoBodySpaceU3ST space(Nmax);
      std::cout << space.Str();
      
      u3shell::OperatorLabelsU3ST op(2,u3::SU3(2,0),0,0,0);

      u3shell::TwoBodySectorsU3ST sectors(space,op);
      
      for (int sector_index=0; sector_index<sectors.size(); ++sector_index)
        {
          auto sector=sectors.GetSector(sector_index);

          std::cout << fmt::format(
                                   "{:3} ({:2},{:2},{:1}): {} {}",
                                   sector_index,
                                   sector.bra_subspace_index(),sector.ket_subspace_index(),sector.multiplicity_index(),
                                   sector.bra_subspace().Str(),sector.ket_subspace().Str()
                                   )
                    << std::endl;
        }
    }

  // termination
  return 0;
}
