/****************************************************************
  indexing_u3st_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "cppformat/format.h"

#include "u3shell/indexing_u3st.h"

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{


  ////////////////////////////////////////////////////////////////
  // two-body state tests
  ////////////////////////////////////////////////////////////////

  if (true)
    {
      // build space
      int Nmax=4;
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
  


  // termination
  return 0;
}
