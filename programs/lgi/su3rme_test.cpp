/****************************************************************
  su3rme_test.cpp

  Anna E. McCoy 
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  10/24/21 (aem): Created.
****************************************************************/
#include "u3ncsm/su3rme.h"

#include <fstream>
#include <iostream>
#include "fmt/format.h"
#include "utilities/utilities.h"
#include "utilities/nuclide.h"

#include "SU3ME/CInteractionPN.h"
#include "SU3ME/InteractionPPNN.h"



int main(int argc, char **argv)
{
  // su3::init();
  nuclide::NuclideType nuclide({3,3});
  int Nmax=2;

  std::string spncci_root_dir=utils::get_spncci_project_root_dir();
  std::string operator_base_name = fmt::format("{}/spncci/data/relative_operators/Brel",spncci_root_dir);
  std::cout<<operator_base_name<<std::endl;

  const auto&[Z,N]=nuclide;
  CBaseSU3Irreps baseSU3Irreps(Z, N, Nmax);
  
  std::ofstream interaction_log_file("/dev/null");
  bool log_is_on = false;
  bool generate_missing_rme = true;
  CInteractionPPNN interactionPPNN(baseSU3Irreps, log_is_on, interaction_log_file);
  std::string ppnn_file_name(operator_base_name+".PPNN");
  
  CInteractionPN interactionPN(baseSU3Irreps, generate_missing_rme, log_is_on, interaction_log_file);
  std::string pn_file_name(operator_base_name+".PN");
  interactionPPNN.LoadTwoBodyOperator(ppnn_file_name);
  interactionPN.AddOperator(pn_file_name);


}//end main
