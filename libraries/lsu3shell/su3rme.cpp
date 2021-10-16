/****************************************************************
  su3rme.cpp

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT
****************************************************************/
#include "lsu3shell/su3rme.h"


namespace lsu3shell
{

  InteractionTerms 
  GetOperatorFromFile(
    const nuclide::NuclideType& nuclide,
    const int Nmax,
    const std::string& operator_filename
  )
  {  
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Load operators for Brel and Nrel 
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    const auto&[Z,N]=nuclide;
    CBaseSU3Irreps baseSU3Irreps(Z, N, Nmax);
    std::ofstream interaction_log_file("/dev/null");
    bool log_is_on = false;

    CInteractionPPNN interactionPPNN(baseSU3Irreps, log_is_on, interaction_log_file);

    bool generate_missing_rme = true;
    CInteractionPN interactionPN(baseSU3Irreps, generate_missing_rme, log_is_on, interaction_log_file);

    std::string ppnn_file_name(operator_filename+".PPNN");
    std::string pn_file_name(operator_filename+".PN");
    interactionPPNN.LoadTwoBodyOperator(ppnn_file_name);
    interactionPN.AddOperator(pn_file_name);

    interactionPPNN.TransformTensorStrengthsIntoPP_NN_structure();

    return InteractionTerms(interactionPPNN,interactionPN);
  }
}// end namespace
