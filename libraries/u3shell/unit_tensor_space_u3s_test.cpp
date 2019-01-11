/****************************************************************
  u3st_scheme_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "u3shell/unit_tensor_space_u3s.h"
#include "fmt/format.h"

#include "u3shell/relative_operator.h"

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  ////////////////////////////////////////////////////////////////
  // Relative state tests
  ////////////////////////////////////////////////////////////////

  if (false)
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

  if(true)
    {
      // build list of observable labels 
      int Nmax=0;
      int N1v=1;
      std::vector<u3shell::RelativeUnitTensorLabelsU3ST> unit_tensor_labels;
      u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax,N1v,unit_tensor_labels,-1,-1,false);
      std::set<u3shell::IndexedOperatorLabelsU3S> observable_labels_set;
      for(auto& unit_tensor : unit_tensor_labels)
        {
          int N0=unit_tensor.N0();
          u3::SU3 x0=unit_tensor.x0();
          HalfInt S0=unit_tensor.S0();
          int kappa0_max=u3::BranchingMultiplicitySO3(x0,0);
          if(kappa0_max>0)
            {
              u3shell::OperatorLabelsU3S labels_u3s(N0,x0,S0);
              
              // kappa0=1 and L0=0
              observable_labels_set.emplace(labels_u3s,1,0);
            }
        }

      std::vector<u3shell::IndexedOperatorLabelsU3S>observable_labels;
      for(auto& labels : observable_labels_set)
        observable_labels.push_back(labels);
      // Construct observable space


    u3shell::ObservableSpaceU3S observable_space(observable_labels);
    std::cout<<observable_space.Str()<<std::endl;
 

    }

  // termination
  return 0;
}
