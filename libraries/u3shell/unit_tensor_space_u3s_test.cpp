/****************************************************************
  u3st_scheme_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "u3shell/unit_tensor_space_u3s.h"
#include "cppformat/format.h"

#include "u3shell/relative_operator.h"

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
      int Nmax=0;
      int N1v=1;
      std::vector<u3shell::RelativeUnitTensorLabelsU3ST> unit_tensor_labels;
      u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax,N1v,unit_tensor_labels,-1,-1,false);

      for(auto& unit_tensor : unit_tensor_labels)
        std::cout<<unit_tensor.Str()<<std::endl;

      u3shell::RelativeUnitTensorSpaceU3S space(Nmax,N1v,unit_tensor_labels);
      std::cout << space.Str();

      // dump subspace contents
      
      for (int subspace_index=0; subspace_index<space.size(); ++subspace_index)
        {
          const u3shell::RelativeUnitTensorSubspaceU3S& subspace = space.GetSubspace(subspace_index);
          std::cout << fmt::format("subspace {}",subspace.Str()) << std::endl;          
        

          for (int state_index=0; state_index<subspace.size(); ++state_index)
            {
          
              const u3shell::RelativeUnitTensorStateU3S state(subspace,state_index);
              std::string state_string = fmt::format("  [{:2} ({:2} {:2}, {:2} {:2})]",
                state.T0(),state.Sp(),state.Tp(),state.S(),state.T());
              std::cout << state_string << std::endl;
            }
        }
    }

  // termination
  return 0;
}
