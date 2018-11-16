/****************************************************************
  truncate_interaction.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  11/2/17 (aem): Created.
****************************************************************/
#include <iostream>
#include <fstream>
#include <unordered_set>
#include "u3shell/relative_branching.h"
#include "cppformat/format.h"
#include "u3shell/interaction_truncation.h"
namespace u3shell
{} // end namespace

// Best result using 1e-2 to 1e-4
double zero_threshold=1e-4;

int main(int argc, char **argv)
{
  u3::U3CoefInit();

  if(argc<4)
  {
    std::cout<<"Nmax Jmax truncation_type (zero_threshold, N_truncation)"<<std::endl;
    std::cout<<"zero_threshold: truncate rmes contributing less than zero_threshold percentage"<<std::endl;
    std::cout<<"N_truncation: maximum allowed value of Nrel for Nrel truncation or N0 for N0 truncation"<<std::endl;
  }

  int N0_max;
  int Nrel_max;
  int x0_max;
  double truncation_threshold;

  int Nmax=std::stoi(argv[1]);
  int Jmax=std::stoi(argv[2]);
  std::string truncation_type=argv[3];

  // depending on decomposition type, initilize on of the truncation parameters 
  if(truncation_type=="Nrel")
    Nrel_max=std::stoi(argv[4]);

  else if(truncation_type=="N0")
    N0_max=std::stoi(argv[4]);

  else if(truncation_type=="x0")
    x0_max=std::stoi(argv[4]);

  else
    truncation_threshold=atof(argv[4]);
  

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Get interaction
  int J0=0;
  int g0=0;
  int T0min=0;
  int T0max=0;
  basis::SymmetryPhaseMode symmetry_phase_mode=basis::SymmetryPhaseMode::kHermitian;
  basis::OperatorLabelsJT operator_labels(J0,g0,T0min,T0max,symmetry_phase_mode);

  // std::string interaction_filename;
  std::string interaction_directory="../../data/relative_interactions";
  std::string interaction_basename="jisp16_Nmax20_hw20.0_rel";
  //std::string interaction_basename="ksqr_Nmax16_rel";
  std::string interaction_filename=fmt::format("{}/{}.dat",interaction_directory, interaction_basename);
  u3shell::RelativeRMEsU3ST interaction_u3st;

  u3shell::GetRelativeRMEsU3ST(interaction_filename,Nmax,Jmax,interaction_u3st);
  u3shell::PrintRelativeRMEsU3ST(interaction_u3st);
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Truncate interaction 
  u3shell::RelativeRMEsU3ST truncated_interaction_u3st;

  if(truncation_type=="u3st")
    u3shell::TruncateInteractionU3ST(interaction_u3st,truncated_interaction_u3st,truncation_threshold);
  
  else if (truncation_type=="u3")
    u3shell::TruncateInteractionByU3(interaction_u3st,truncated_interaction_u3st,truncation_threshold);
 
  else if (truncation_type=="Nrel")
    TruncateInteractionNrel(interaction_u3st,truncated_interaction_u3st,Nrel_max);

  else if (truncation_type=="N0")
    TruncateInteractionN0(interaction_u3st,truncated_interaction_u3st,N0_max); 

  else if (truncation_type=="x0")
    TruncateInteractionx0(interaction_u3st,truncated_interaction_u3st,x0_max);   

  else 
    {
      std::cout<<fmt::format("{} is an invalid truncation type. ", truncation_type)<<std::endl;
      std::cout<<"valid truncation types are u3st, u3, Nrel and N0"<<std::endl;
      std::exit(EXIT_FAILURE);
    }

  // u3shell::PrintRelativeRMEsU3ST(truncated_interaction_u3st);
  std::string rme_filename;
  if(truncation_type=="N0" )
    rme_filename=fmt::format("u3st_{}_{:02d}_double_{:.1e}.dat",truncation_type, N0_max, zero_threshold);  
  else if(truncation_type=="Nrel")
    rme_filename=fmt::format("u3st_{}_{:02d}_double_{:.1e}.dat",truncation_type, Nrel_max, zero_threshold);
  else if(truncation_type=="x0")
    rme_filename=fmt::format("u3st_{}_{:02d}_double_{:.1e}.dat",truncation_type, x0_max, zero_threshold);
  else
    rme_filename=fmt::format("u3st_{}_{:.1e}_double_{:.1e}.dat",truncation_type, truncation_threshold, zero_threshold);
  //u3shell::WriteRelativeRMEsU3ST(rme_filename, truncated_interaction_u3st);
  u3shell::WriteRelativeOperatorU3ST(rme_filename, truncated_interaction_u3st, true);
  std::cout<<"number of rmes u3st truncated "<<truncated_interaction_u3st.size()<<std::endl;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Write out branched relative rmes 
  basis::RelativeSpaceLSJT relative_space_lsjt(Nmax, Jmax);
  std::array<basis::RelativeSectorsLSJT,3> isospin_component_sectors_lsjt;
  std::array<basis::MatrixVector,3> isospin_component_blocks_lsjt;

  u3shell::BranchRelativeRMEs(operator_labels,Nmax,Jmax, 
    truncated_interaction_u3st,
    relative_space_lsjt,
    isospin_component_sectors_lsjt,
    isospin_component_blocks_lsjt
    );

  std::string relative_filename;
  if(truncation_type=="N0" )
    relative_filename=fmt::format("{}_Nmax{:02d}_{}_{:02d}.dat",interaction_basename,Nmax,truncation_type,N0_max);  
  else if(truncation_type=="Nrel")
    relative_filename=fmt::format("{}_Nmax{:02d}_{}_{:02d}.dat",interaction_basename,Nmax,truncation_type,Nrel_max);
  else if(truncation_type=="x0")
    relative_filename=fmt::format("{}_Nmax{:02d}_{}_{:02d}.dat",interaction_basename,Nmax,truncation_type,x0_max);    
  else
    relative_filename=fmt::format("{}_Nmax{:02d}_{}_{:0.1e}.dat",interaction_basename,Nmax,truncation_type,truncation_threshold);
 
  bool verbose=true;
  basis::WriteRelativeOperatorLSJT(
    relative_filename,
    relative_space_lsjt,
    operator_labels,
    isospin_component_sectors_lsjt,
    isospin_component_blocks_lsjt,
    verbose
  );
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}



