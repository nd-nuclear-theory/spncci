/****************************************************************
  interaction_truncation_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  11/2/17 (aem): Created.
****************************************************************/
#include <iostream>
#include <fstream>
#include <unordered_set>
#include "u3shell/relative_branching.h"
#include "am/halfint.h"
#include "am/halfint_fmt.h"
#include "fmt/format.h"
#include "u3shell/interaction_truncation.h"

double zero_threshold=1e-6;

int main(int argc, char **argv)
{
  u3::U3CoefInit();

  int Nmax=10;
  int Jmax=4;
  double truncation_threshold=1e-2;

  int J0=0;
  int g0=0;
  int T0min=0;
  int T0max=0;
  basis::SymmetryPhaseMode symmetry_phase_mode=basis::SymmetryPhaseMode::kHermitian;
  basis::OperatorLabelsJT operator_labels2(J0,g0,T0min,T0max,symmetry_phase_mode);

  // std::string interaction_filename;
  std::string interaction_filename="../../data/relative_interactions/jisp16_Nmax20_hw20.0_rel.dat";
  u3shell::RelativeRMEsU3ST interaction_u3st;

  u3shell::GetRelativeRMEsU3ST(interaction_filename,Nmax,Jmax,interaction_u3st);
  // u3shell::PrintRelativeRMEsU3ST(interaction_u3st)

  u3shell::RelativeRMEsU3ST truncated_interaction_u3st;
  u3shell::TruncateInteractionU3ST(interaction_u3st,truncated_interaction_u3st,truncation_threshold);
  // u3shell::PrintRelativeRMEsU3ST(truncated_interaction_u3st);
  std::cout<<"number of rmes u3st truncated "<<truncated_interaction_u3st.size()<<std::endl;

   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // get distribution of interaction over U(3) tensors
  if(false)
    {
      std::unordered_map<u3shell::SpacialSymmetries,double,boost::hash<u3shell::SpacialSymmetries>> probability_by_u3;
      u3shell::AccumulateInteractionByU3(interaction_u3st,probability_by_u3);

      double total;
      for(auto it=probability_by_u3.begin(); it!=probability_by_u3.end(); ++it)
        {
          u3::SU3 x0;
          HalfInt N0;
          std::tie(N0,x0)=it->first;
          if(N0>=0)
            total+=it->second;
        }

      for(auto it=probability_by_u3.begin(); it!=probability_by_u3.end(); ++it)
        {
          u3::SU3 x0;
          HalfInt N0;
          std::tie(N0,x0)=it->first;
          if(N0>=0)
            std::cout<<fmt::format("{:3} {:8}    {:10f}",N0,x0.Str(),it->second/total)<<std::endl;
        }
    }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    u3shell::RelativeRMEsU3ST truncated_interaction_u3;
    u3shell::TruncateInteractionByU3(interaction_u3st,truncated_interaction_u3,truncation_threshold);
      std::cout<<"number of rmes u3 truncated "<<truncated_interaction_u3.size()<<std::endl;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Write out branched relative rmes
  if(false)
    {
      basis::RelativeSpaceLSJT relative_space_lsjt2(Nmax, Jmax);
      std::array<basis::RelativeSectorsLSJT,3> isospin_component_sectors_lsjt2;
      std::array<basis::MatrixVector,3> isospin_component_blocks_lsjt2;

      std::string relative_filename=fmt::format("test_{:2d}_{:e}",Nmax,truncation_threshold);
      bool verbose=true;
      basis::WriteRelativeOperatorLSJT(relative_filename,
          relative_space_lsjt2,
          operator_labels2,
          isospin_component_sectors_lsjt2,
          isospin_component_blocks_lsjt2,
          verbose
        );
    }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
