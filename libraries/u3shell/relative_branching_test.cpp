/****************************************************************
  relative_branching_test.cpp

  Anna E. McCoy
  TRIUMF

  SPDX-License-Identifier: MIT

  11/2/17 (aem): Created.
  5/14/19 (aem): Updated   basis::OperatorLabelsJT to 
    basis::RelativeOperatorParametersLSJT for reading in operators
****************************************************************/
#include <iostream>
#include <fstream>
#include "u3shell/relative_branching.h"

double zero_threshold=10e-6;
namespace u3shell
{
} // end namespace


int main(int argc, char **argv)
{
  u3::U3CoefInit();

  std::string interaction_filename="../../data/relative_interactions/jisp16_Nmax20_hw20.0_rel.dat";
  int Nmax=4;
  int Jmax=4;

  int J0=0;
  int g0=0;
  int T0min=0;
  int T0max=0;
  
  basis::SymmetryPhaseMode symmetry_phase_mode=basis::SymmetryPhaseMode::kHermitian;
  basis::OperatorLabelsJT operator_labels2(J0,g0,T0min,T0max,symmetry_phase_mode);

  // Read in the interaction from file
  basis::RelativeSpaceLSJT relative_space_lsjt;
  std::array<basis::RelativeSectorsLSJT,3> isospin_component_sectors_lsjt;
  std::array<basis::MatrixVector,3> isospin_component_blocks_lsjt;
  basis::RelativeOperatorParametersLSJT operator_labels;

  // Reads in relative operator and fills out isospin_component_sectors_lsjt, 
  // isospin_component_blocks_lsjt, operator_labels and relative_space_lsjt.
  basis::ReadRelativeOperatorLSJT(
    interaction_filename,relative_space_lsjt,operator_labels,
    isospin_component_sectors_lsjt, isospin_component_blocks_lsjt, true
    );

  // upcouple interaction
  u3shell::RelativeRMEsU3ST interaction_u3st;
  std::cout<<"upcoupling "<<std::endl; 
  for (int T0=T0min; T0<=T0max; ++T0)
    {
      u3shell::Upcoupling(
        relative_space_lsjt,
        isospin_component_sectors_lsjt,
        isospin_component_blocks_lsjt,
        J0, g0, T0,Nmax, interaction_u3st);
    }


  basis::RelativeSpaceLSJT relative_space_lsjt2(Nmax, Jmax);
  std::array<basis::RelativeSectorsLSJT,3> isospin_component_sectors_lsjt2;
  std::array<basis::MatrixVector,3> isospin_component_blocks_lsjt2;
  u3shell::BranchRelativeRMEs(
    operator_labels2,Nmax, Jmax, 
    interaction_u3st,relative_space_lsjt2,
    isospin_component_sectors_lsjt2,
    isospin_component_blocks_lsjt2
    );


  // Upcouple branched rmes
  u3shell::RelativeRMEsU3ST interaction_u3st2;
  std::cout<<"upcoupling "<<std::endl; 
  for (int T0=T0min; T0<=T0max; ++T0)
    {
      u3shell::Upcoupling(
        relative_space_lsjt2,
        isospin_component_sectors_lsjt2,
        isospin_component_blocks_lsjt2,
        J0, g0, T0,Nmax, interaction_u3st2);
    }

  // Compare upcoupled rmes from file with upcoupled rmes from branched rmes
  for(auto it=interaction_u3st.begin(); it!=interaction_u3st.end(); ++it)
    std::cout<<it->second<<"  "<< interaction_u3st2[it->first]<<"  "<<it->second/interaction_u3st2[it->first]<< std::endl;

  // Write branched rmes to file
  std::string relative_filename="test";
  bool verbose=true;
  basis::WriteRelativeOperatorLSJT(relative_filename,
      relative_space_lsjt2,
      operator_labels2,
      isospin_component_sectors_lsjt2,
      isospin_component_blocks_lsjt2,
      verbose
    );

}



